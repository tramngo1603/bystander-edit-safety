"""Steps 5 & 6: safety scoring + ranking (Step 5), and the credibility layer (Step 6).

Both work off the same Steps 2/4 enumeration output and the same RubricWeights. Step 5 ranks the
therapeutic designs by the safety rubric. Step 6 shows the ranking diverges from naive sorts and
attaches conservative per-guide uncertainty. WEIGHTS ARE PROVISIONAL (see scoring/rubric.py) and are
read from a RubricWeights config object, so you can retune them without touching pipeline code.
"""
from __future__ import annotations

import json

from ..io import outputs
from ..scoring.credibility import annotate_uncertainty, divergence_report
from ..scoring.rubric import EFFICACY_CAVEAT, RubricWeights, WEIGHTS_STATUS, score_guides


def _load_enumeration() -> list[dict]:
    with open(outputs.enumeration_path()) as handle:
        return json.load(handle)["records"]


def run_ranking(record_dicts: list[dict] | None = None,
                weights: RubricWeights | None = None, write: bool = True) -> list[dict]:
    """Step 5: rank designs by safety, return their decomposed score dicts (best first)."""
    weights = weights or RubricWeights()
    if record_dicts is None:
        record_dicts = _load_enumeration()

    ranked = [s.as_dict() for s in score_guides(record_dicts, weights)]
    if write:
        outputs.write_ranking(
            ranked,
            filename="ranking_HBG1_HBG2.json",
            meta={
                "step": "5 (safety scoring + ranking)",
                "weights_status": WEIGHTS_STATUS,
                "weights": weights.as_dict(),
                "efficacy_caveat": EFFICACY_CAVEAT,
                "bystander_scoring": "FUNCTIONAL (motif consequence), not bystander count",
                "offtarget_scoring": "mismatch-stratified (4-mm down-weighted); 0-mm paralog excluded; "
                                     "no hard count/VAF gate (would gate on noise)",
                "ch_gate": "predicted CH-driver off-target vetoes (worst tier) ONLY if low-mismatch AND "
                           "exonic/UTR; high-mismatch or intronic CH hits are downgraded to a noted "
                           "flag. Nomination requiring empirical confirmation, not a verdict.",
                "consequence_slot": "active (w_consequence=2.0): penalizes the predicted regulatory "
                                    "severity of non-productive bystanders once Step 3 "
                                    "(annotate_consequence.py) fills it in. The default enumeration leaves "
                                    "consequence 'deferred', so the term is zero until that pass runs; "
                                    "the CH functional-consequence field is sharpened by the same pass.",
                "ranking_unit": "design = guide+editor that installs an in-scope anchor",
                "score_interpretation": "Read the tier and the decomposed `terms` first. safety_score is "
                                        "a PROVISIONAL within-tier ordering signal (weights are provisional "
                                        "calibration anchors), not an absolute safety measure and not a clearance.",
            },
        )
    return ranked


def run_credibility(record_dicts: list[dict] | None = None,
                    weights: RubricWeights | None = None, write: bool = True) -> dict:
    """Step 6: divergence vs naive baselines + conservative per-guide uncertainty."""
    weights = weights or RubricWeights()
    if record_dicts is None:
        record_dicts = _load_enumeration()

    scores = [s.as_dict() for s in score_guides(record_dicts, weights)]
    divergence = divergence_report(record_dicts, weights)
    annotated = [{**s, "uncertainty": annotate_uncertainty(s)} for s in scores]

    body = {
        "divergence": divergence,
        "ranked_designs_with_uncertainty": annotated,
    }
    if write:
        outputs.write_credibility(
            body,
            filename="credibility_HBG1_HBG2.json",
            meta={
                "step": "6 (credibility core)",
                "weights_status": WEIGHTS_STATUS,
                "naive_baselines": ["on-target efficiency alone", "flat off-target count alone"],
                "uncertainty": "conservative; model-predicted/nominated, not measured; testable predictions",
                "human_gated": "results presented for review; not finalized unilaterally",
            },
        )
    return body


if __name__ == "__main__":
    # Step 5 (ranking) then Step 6 (credibility) off the same enumeration output.
    ranked = run_ranking()
    print(f"WEIGHTS: {WEIGHTS_STATUS}")
    print(f"ranked designs: {len(ranked)}\n")
    hdr = "{:>3} {:<26} {:<18} {:<13} {:>6} {:>8} {:>9} {:>10} {:>8}".format(
        "#", "tier", "guide", "anchors", "eff", "bystPen", "offtPen", "CHflag", "score")
    print(hdr); print("-" * len(hdr))
    for s in ranked:
        chf = (",".join(s["ch_genes"]) + ("*GATE" if s["ch_gate"] else ("*noted" if s["ch_flag_noted"] else ""))) if s["ch_genes"] else "-"
        print("{:>3} {:<26} {:<18} {:<13} {:>6} {:>8} {:>9} {:>10} {:>8}".format(
            s["rank"], s["tier_label"], s["guide_id"], str(s["anchors"]),
            s["on_target_efficiency"], s["terms"]["bystander_penalty"],
            s["terms"]["offtarget_penalty"], chf, s["safety_score"]))

    body = run_credibility()
    d = body["divergence"]
    print(f"\n=== DIVERGENCE: {d['n_designs']} designs ===")
    print(f"changed vs naive-efficiency: {d['n_changed_vs_naive_efficiency']} | "
          f"changed vs naive-flat-count: {d['n_changed_vs_naive_offtarget_count']}\n")
    print("biggest movers vs naive-efficiency (rubric_rank <- naive_eff_rank):")
    for m in d["biggest_movers_vs_naive_efficiency"]:
        print(f"  {m['guide_id']:<20} rubric#{m['rubric_rank']:<3} naiveEff#{m['naive_efficiency_rank']:<3} "
              f"naiveCnt#{m['naive_offtarget_count_rank']:<3} | {m['rationale']}")
