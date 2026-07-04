"""Build one design record per candidate guide.

Each record collects what we know about a candidate: the predicted on-target and bystander edits, the
off-target classes, the HBG paralog hazard, the CH-driver annotation, how it matches the HBG evidence
registry, its evidence tier, and the lab tests needed to confirm it. It is meant to be read, not
reduced to a single number. The safety score is included too, but only as a secondary, provisional
tiebreaker.

A design that installs more than one edit (for example the -123/-124 pair) keeps every intended
on-target edit, and the registry match runs over all of them, not just the highest-rate one.
"""
from __future__ import annotations

from . import annotation, evidence, hbg_registry, offtarget_taxonomy, paralog_hazard
from .io import step0_anchor as anchor
from .pipeline import environment
from .scoring.rubric import RubricWeights, score_guides

SCHEMA_VERSION = "1"
_MECHANISM_BY_OFFSET = {p.promoter_offset: p.mechanism for p in anchor.IN_SCOPE}
_REL_ORDER = (hbg_registry.NO_MATCH, hbg_registry.NEARBY,
              hbg_registry.SAME_MOTIF, hbg_registry.EXACT_EDIT)


def _edit_str(record: dict) -> dict:
    return {
        "genomic_plus": f"{record['genomic_plus_ref']}>{record['genomic_plus_alt']}",
        "promoter_sense": f"{record['promoter_sense_ref']}>{record['promoter_sense_alt']}",
    }


def _strongest(values, order) -> str:
    present = [v for v in values if v in order]
    return max(present, key=order.index) if present else order[0]


def _validation_checklist(bystanders, off_triggers, hazard, ch_annot) -> list[dict]:
    """Assemble the mandatory empirical-validation checklist for a candidate."""
    items = [{"validates": "on_target_edit_rate",
              "assay": "amplicon NGS in the target cell type (CD34+ HSPCs)", "status": "required"}]
    if bystanders:
        items.append({"validates": "bystander_spectrum",
                      "assay": "amplicon NGS quantifying bystander edits across the editing window",
                      "status": "required"})
    items.extend(off_triggers)
    items.extend(hazard.get("required_confirmation", []))
    if ch_annot:
        items.append({"validates": "ch_driver_offtarget",
                      "assay": "targeted off-target sequencing + functional confirmation at the CH-driver locus",
                      "status": "required"})
    return items


def _limitations() -> list[str]:
    return [
        "On-target and bystander rates are model-predicted (BE-DICT), not measured in HSCs.",
        "Off-target nominations are Cas9-style reference-genome homology only; base-editor-specific, "
        "guide-independent deaminase, and RNA off-target classes are not assessed here.",
        "CH annotation is gene-membership plus coordinate/mechanism context for documentation only; it "
        "produces no clonal-expansion, leukemia, incidence, time-to-event, or clinical estimate (G3).",
        "Registry relationships to clinical programs are motif/mechanism analogs; exact clinical target "
        "positions may be non-public, and published literature never promotes the evidence tier.",
        "Scope is the HBG1/HBG2 proximal promoter only.",
    ]


def build_design_records(records: list[dict], weights: RubricWeights | None = None,
                         empirical: dict | None = None) -> dict:
    """Build the HBG design-history artifact set from annotated enumeration records (Steps 2/4).

    records:   enumeration EditRecords-as-dicts (on-target + bystander) with provenance.
    empirical: optional {guide_id: [measurement, ...]} that elevates a candidate to
               empirically_measured (measured in this run/context). Absent empirical evidence never
               becomes a model value, and published literature never promotes the tier.
    """
    weights = weights or RubricWeights()
    empirical = empirical or {}
    scores = {s.guide_id: s for s in score_guides(records, weights)}

    by_guide: dict[str, list[dict]] = {}
    for r in records:
        by_guide.setdefault(r["provenance"]["guide_id"], []).append(r)

    out = []
    for guide_id, recs in by_guide.items():
        anchors = sorted((r for r in recs if r["on_target_or_bystander"] == "on_target"),
                         key=lambda r: abs(r["promoter_offset"]))
        if not anchors:
            continue  # only guides that install an in-scope anchor become design records
        bystanders = [r for r in recs if r["on_target_or_bystander"] == "bystander"]
        prov = recs[0]["provenance"]
        ch_hits = prov.get("ch_driver_offtargets", [])
        score = scores.get(guide_id)
        base = anchors[0]

        # Registry match over the SET of on-target anchors (biology of e.g. -123/-124 preserved).
        by_edit = []
        for a in anchors:
            mech = _MECHANISM_BY_OFFSET.get(abs(a["promoter_offset"]))
            m = hbg_registry.relationship(a["gene"], a["promoter_offset"], a["editor_class"], mech)
            by_edit.append({"promoter_offset": a["promoter_offset"], "relationship": m["relationship"],
                            "evidence_strength": m["evidence_strength"],
                            "matched_entries": m["matched_entries"]})
        best_rel = _strongest([m["relationship"] for m in by_edit], _REL_ORDER)
        best_strength = _strongest([m["evidence_strength"] for m in by_edit],
                                   list(hbg_registry.STRENGTH_ORDER))
        empirical_support = empirical.get(guide_id, [])
        tier = evidence.evidence_tier(best_rel, empirical_support)
        measured_in = "this_run" if tier == evidence.EMPIRICALLY_MEASURED else "literature"

        taxonomy = offtarget_taxonomy.classify_offtargets(prov.get("offtarget", {}))
        off_triggers = offtarget_taxonomy.assay_triggers(taxonomy)
        hazard = paralog_hazard.paralog_hazard(base)
        ch_annot = [
            {**{k: h.get(k) for k in ("ch_driver_gene", "chrom", "start", "mismatch_count", "feature")},
             "consequence_class": annotation.ch_consequence_class(h.get("ch_driver_gene"), h.get("feature"))}
            for h in ch_hits
        ]

        out.append({
            "record_type": "hbg_design_history",
            "schema_version": SCHEMA_VERSION,
            "candidate": {
                "guide_id": guide_id, "gene": base["gene"], "editor_class": base["editor_class"],
                "editor_name": prov.get("editor"),
                "target": {"locus": "HBG1/HBG2 proximal promoter", "chrom": base["chrom"],
                           "intended_edits": [{"promoter_offset": a["promoter_offset"],
                                               "pos_1based": a["pos_1based"],
                                               "mechanism": _MECHANISM_BY_OFFSET.get(abs(a["promoter_offset"])),
                                               "edit": _edit_str(a)} for a in anchors]},
            },
            "evidence_tier": tier,
            "registry_match": {
                "relationship": best_rel,
                "evidence_strength": best_strength,
                "measured_in": measured_in,
                "by_edit": by_edit,
                "note": ("Relationship and evidence strength are separate. Published literature "
                         "(including exact published HSPC edits) is analog_supported at most; "
                         "empirically_measured means measured in this run/context."),
            },
            "predicted_edit_spectrum": {
                "on_target": [{"promoter_offset": a["promoter_offset"], "predicted_rate": a.get("predicted_rate"),
                               "model": a.get("model_name"), "confidence": a.get("confidence")} for a in anchors],
                "bystanders": [{"promoter_offset": b["promoter_offset"], "predicted_rate": b.get("predicted_rate"),
                                "genomic_plus": _edit_str(b)["genomic_plus"]} for b in bystanders],
                "model_disagreement": {"status": "single_predictor",
                                       "note": "one predictor (BE-DICT) wired; no cross-model comparison available"},
            },
            "off_target": {"taxonomy": taxonomy, "assay_triggers": off_triggers},
            "hbg_paralog_hazard": hazard,
            "ch_annotation": {"driver_offtargets": ch_annot,
                              "note": "documentation/monitoring only; no clonal-expansion or clinical estimate (G3)"},
            "empirical_support": empirical_support,
            "required_validation": _validation_checklist(bystanders, off_triggers, hazard, ch_annot),
            "provisional_scoring": {
                "note": "SECONDARY and PROVISIONAL; the evidence vectors above are the primary read, not this score.",
                **(score.as_dict() if score else {}),
            },
            "provenance": {"guide_id": guide_id, "model_name": base.get("model_name"),
                           "predictor_runs": prov.get("predictor_runs"),
                           "protospacer": prov.get("protospacer"), "pam": prov.get("pam"),
                           "strand": prov.get("strand")},
            "limitations": _limitations(),
        })

    return {
        "record_type": "hbg_design_history_set",
        "schema_version": SCHEMA_VERSION,
        "count": len(out),
        "records": out,
        "reproducibility": {
            "generated_at": environment.run_stamp(),
            "repo": environment.repo_state(),
            "genome_build": anchor.GENOME_BUILD,
            "models": sorted({r.get("model_name") for r in records if r.get("model_name")}),
            "config": {"rubric_weights": weights.as_dict()},
        },
    }
