"""Build a display-ready view of the bystander analysis (read-only; does not change the pipeline).

Reads the existing enumeration output, the locked Step 0 anchor table, and the rubric's functional
bystander annotations, then writes a self-contained JSON for a display page. A bystander is an
accidental edit NEAR the intended target, inside the editing window; this is NOT genome-wide
off-target analysis. Every predicted rate is labelled model-predicted, not HSC-measured. It uses
only real values from the enumeration; a missing value becomes null and is never made up. It
produces no clonal or lifetime-risk value (G3).
"""
from __future__ import annotations

import hashlib
import json

from . import outputs
from . import step0_anchor as anchor
from ..scoring.rubric import bystander_function

_ENUMERATION_FILE = outputs.enumeration_path()
_DISPLAY_FILE = outputs.output_path("bystander_display.json")

EDITOR_NOTE = ("bystander = accidental edits NEAR the target (within the editing window); "
               "this is NOT off-target (genome-wide) analysis")

# Primary engineering anchors (docs/step0-coordinates.md): strongest exact-allele HSPC evidence.
PRIMARY_ANCHORS = {175, 113, 123, 124}

# Display groupings of the in-scope positions (the -123/-124 synthetic pair shown together).
TARGET_GROUPS = [
    ("-198", [198]),
    ("-175", [175]),
    ("-113", [113]),
    ("-123/-124", [123, 124]),
    ("-117", [117]),
    ("-114", [114]),
]

# Per-target editorial notes, grounded in the locked anchor table / rubric annotations only.
TARGET_NOTES = {
    "-123/-124": ("Synthetic double-mutant: KLF1 binding was shown for the -123/-124 PAIR "
                  "(Ravi 2022), not either single position; both are intended edits here."),
    "-175": ("Dominant bystander is functionally adverse (TAL1 E-box). Editor choice is a tradeoff: "
             "a higher-bystander editor (ABE8e) can still be the better therapeutic (Mayuranathan "
             "2023); evaluate empirically per editor (only ABE8e was run here)."),
    "-113": ("Adjacent bystander -112 is functionally adverse (destroys the created GATA1 motif); "
             "-116 is co-productive. Function, not count, drives the bystander concern."),
    "-114": ("Adjacent C edits (-114/-115) are part of the intended distal-BCL11A-motif disruption "
             "(productive bystander)."),
}


def _rate(value, model_name):
    """Label a predicted rate as model-predicted (not measured); return null if the value is absent."""
    return {
        "value": value,
        "model_predicted": True,
        "measured": False,
        "source": model_name,
        "note": "BE-DICT prediction in cell-line/library context; not HSC-measured",
    }


def _anchor_by_offset() -> dict:
    return {p.promoter_offset: p for p in anchor.IN_SCOPE}


def _functional_source(promoter_offset: int) -> str | None:
    cls, _ = bystander_function(promoter_offset)
    if cls == "unknown":
        return None
    return ("rubric FUNCTIONAL_BYSTANDER (EMERGING; Mayuranathan 2023 PMID 37400614 / "
            "Martyn 2019 PMID 30617196 / Wienert 2015 PMID 25971621)")


def build_display(records: list[dict]) -> dict:
    anchors = _anchor_by_offset()

    # Index HBG1 records by guide_id (HBG2 mirrors via homology; rates are per-protospacer).
    by_guide: dict[str, list[dict]] = {}
    for r in records:
        if r["gene"] != "HBG1":
            continue
        by_guide.setdefault(r["provenance"]["guide_id"], []).append(r)

    def representative_guide(offsets: list[int]) -> str | None:
        best, best_rate = None, None
        for gid, recs in by_guide.items():
            ot = [x for x in recs if x["on_target_or_bystander"] == "on_target"
                  and x["promoter_offset"] in offsets]
            if not ot:
                continue
            rate = max((x["predicted_rate"] or 0.0) for x in ot)
            if best is None or rate > best_rate:
                best, best_rate = gid, rate
        return best

    targets = []
    for label, offsets in TARGET_GROUPS:
        ap = anchors[offsets[0]]
        gid = representative_guide(offsets)
        recs = by_guide.get(gid, []) if gid else []
        model_name = recs[0]["model_name"] if recs else None

        intended = []
        for off in offsets:
            a = anchors[off]
            ot = [x for x in recs if x["on_target_or_bystander"] == "on_target"
                  and x["promoter_offset"] == off]
            rate_val = ot[0]["predicted_rate"] if ot else None
            intended.append({
                "position": f"-{off}",
                "genomic": {"chrom": anchor.CHROM,
                            "hbg1_pos": a.hbg1_pos, "hbg2_pos": a.hbg2_pos,
                            "ref_alt": f"{a.genomic_plus_ref}>{a.genomic_plus_alt}"},
                "sense": f"{a.promoter_sense_ref}>{a.promoter_sense_alt} (promoter/sense strand)",
                "predicted_rate": _rate(rate_val, model_name),
            })

        bystanders = []
        for x in recs:
            if x["on_target_or_bystander"] != "bystander":
                continue
            off = x["promoter_offset"]
            cls, mech = bystander_function(off)
            bystanders.append({
                "position": f"-{off}",
                "genomic": {"chrom": x["chrom"], "hbg1_pos": x["pos_1based"],
                            "ref_alt": f"{x['genomic_plus_ref']}>{x['genomic_plus_alt']}"},
                "sense": f"{x['promoter_sense_ref']}>{x['promoter_sense_alt']} (promoter/sense strand)",
                "predicted_rate": _rate(x["predicted_rate"], x["model_name"]),
                "functional_class": cls,
                "mechanism": mech,
                "functional_source": _functional_source(off),
            })
        bystanders.sort(key=lambda b: -(b["predicted_rate"]["value"] or 0.0))

        targets.append({
            "label": label,
            "gene_loci": ["HBG1", "HBG2"],
            "editor_class": ap.editor_class,
            "intended_mechanism": ap.mechanism,
            "is_primary_anchor": any(o in PRIMARY_ANCHORS for o in offsets),
            "variant_db_status": {"HBG1": anchors[offsets[0]].hbg1_status,
                                  "HBG2": anchors[offsets[0]].hbg2_status},
            "representative_guide": gid,
            "intended_edits": intended,
            "bystanders": bystanders,
            "note": TARGET_NOTES.get(label),
        })

    return {
        "meta": {
            "generated_from": {
                "enumeration_file": "data/outputs/enumeration_HBG1_HBG2.json",
                "run_ref_sha256": _enumeration_sha256(),
                "pipeline_step": "view export of Steps 2/4 results (no pipeline change)",
            },
            "claim_disclaim": outputs.CLAIM_DISCLAIM,
            "editor_note": EDITOR_NOTE,
        },
        "targets": targets,
    }


def _enumeration_sha256() -> str:
    with open(_ENUMERATION_FILE, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()[:16]


def write_display() -> str:
    with open(_ENUMERATION_FILE) as handle:
        records = json.load(handle)["records"]
    display = build_display(records)
    with open(_DISPLAY_FILE, "w") as handle:
        json.dump(display, handle, indent=2)
    return _DISPLAY_FILE


if __name__ == "__main__":
    print(write_display())
