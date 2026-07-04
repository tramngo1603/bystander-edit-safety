"""Step 3 activation: fill in the consequence layer with the external regulatory model.

This reads the enumeration output, scores every edit's regulatory consequence with AlphaGenome
(cached, deduped), scores each CH-driver off-target locus to sharpen the CH gate (functional vs
benign), and writes the enumeration output back with consequence filled in (it replaces the deferred
passthrough). It reads the API key from the environment only and never saves it. Coding off-target
sites would use AlphaMissense (out of scope in this build; source noted in annotation.py) instead of the regulatory model; none come up for the
promoter edits here.
"""
from __future__ import annotations

import json

from ..adapters.consequence import AlphaGenomeConsequenceAnnotator
from ..io import outputs

_ENUMERATION_FILE = outputs.enumeration_path()


def run_consequence(annotator: AlphaGenomeConsequenceAnnotator | None = None,
                    write: bool = True, score_ch_offtargets: bool = True) -> dict:
    """Annotate the enumeration records' consequence (and CH off-target consequence) in place."""
    annotator = annotator or AlphaGenomeConsequenceAnnotator()
    with open(_ENUMERATION_FILE) as handle:
        doc = json.load(handle)
    records = doc["records"]

    # 1) regulatory consequence per edit (genomic-plus allele; noncoding/regulatory).
    n_scored = 0
    for r in records:
        r["consequence"] = annotator.score(
            r["chrom"], r["pos_1based"], r["genomic_plus_ref"], r["genomic_plus_alt"])
        r["provenance"] = {
            **r["provenance"],
            "consequence_source": annotator.name,
            "consequence_note": "regulatory consequence MODEL-PREDICTED (AlphaGenome); requires wet-lab confirmation",
        }
        n_scored += 1

    # 2) CH-driver off-target functional consequence: an APPROXIMATE, clearly-labelled proxy.
    # The true off-target base at a multi-mismatch site is not modelled, so we score a fixed A>G proxy
    # edit at the locus. This is a locus-level regulatory-sensitivity signal, not the real off-target
    # consequence; it is marked is_proxy=True and, by default, does NOT harden the CH gate on its own
    # (see RubricWeights.ch_gate_trust_proxy_consequence). Pass score_ch_offtargets=False to skip it.
    n_ch = 0
    if score_ch_offtargets:
        for r in records:
            for hit in r["provenance"].get("ch_driver_offtargets", []):
                if hit.get("ch_functional_consequence"):
                    continue
                cons = annotator.score(hit["chrom"], int(hit["start"]), "A", "G")
                hit["ch_functional_consequence"] = {
                    **cons,
                    "is_proxy": True,
                    "proxy_edit": "A>G",
                    "note": "APPROXIMATE locus-level proxy edit (A>G) at the off-target site; the true "
                            "off-target base is not modelled. Regulatory severity only.",
                }
                n_ch += 1

    doc["meta"] = {**doc.get("meta", {}),
                   "consequence_layer": "ACTIVE (AlphaGenome, non-commercial; outputs not used to train models)",
                   "consequence_scored_edits": n_scored,
                   "ch_offtarget_consequence_scored": n_ch}
    if write:
        with open(_ENUMERATION_FILE, "w") as handle:
            json.dump(doc, handle, indent=2)
    return {"scored_edits": n_scored, "ch_offtargets_scored": n_ch}


if __name__ == "__main__":
    print(run_consequence())
