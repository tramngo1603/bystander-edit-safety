"""Step 6 credibility layer: divergence vs naive baselines, and conservative per-guide uncertainty.

Two views over the same ranked designs. divergence_report shows the rubric ranks differently from
naive single-dimension sorts, so the rubric's judgment is doing real work. annotate_uncertainty
attaches a conservative, testable uncertainty block per guide: everything is model-predicted or
nominated, never measured, and each item names the assay that would confirm or refute it. No
clonal/lifetime-risk value is produced (G3).
"""
from __future__ import annotations

from .rubric import GuideScore, RubricWeights, score_guides


# --- Conservative per-guide uncertainty -----------------------------------------------------------

# Systemic caveats: true for every design (limits of the prediction sources themselves).
SOURCE_CAVEATS = {
    "on_target_edit_rate": (
        "MODEL-PREDICTED (BE-DICT bystander model, trained on cell-line/library data), NOT measured "
        "in hematopoietic stem cells. Confirm by amplicon NGS in the target cell type (e.g. CD34+ HSPCs)."
    ),
    "off_target_nominations": (
        "NOMINATED by reference-only, ungapped, <=4-mismatch CRISPOR search (no bulges, no population "
        "variants; for population/ancestry-aware nomination see the companion ancestry-aware-offtarget "
        "project). Not empirically confirmed; counts are dominated by low-likelihood 4-mismatch sites "
        "and are likelihood-weighted in scoring. Confirm by an empirical assay (CHANGE-seq-BE / GUIDE-seq)."
    ),
    "regulatory_consequence": (
        "MODEL-PREDICTED by AlphaGenome (non-commercial; regulatory-perturbation magnitude only, "
        "not a direction and not an outcome estimate), NOT measured. Confirm by a functional/reporter assay."
    ),
    "ch_feature_call": (
        "Exon/intron/intergenic is COORDINATE annotation; the CH off-target's functional consequence "
        "is now MODEL-PREDICTED by AlphaGenome (regulatory), still not measured. A noted CH off-target "
        "requires targeted off-target sequencing and functional confirmation."
    ),
    "efficacy_floor": (
        "Tier anchors are ALLELIC-editing values in mixed peripheral-blood/marrow cells; allelic % "
        "does NOT convert to edited-HSC fraction. Re-estimate with the final guide/editor/process and "
        "a real potency assay."
    ),
}

OVERALL = "MODEL-PREDICTED / NOMINATED: requires wet-lab confirmation; not a measured or clinical result"


def annotate_uncertainty(score: dict) -> dict:
    """Build the conservative uncertainty block for one ranked design (a GuideScore.as_dict())."""
    testable = []

    eff = score["on_target_efficiency"]
    band = score["efficacy_band"]
    testable.append(
        f"Predicted on-target editing ~{round(eff * 100, 1)}% at anchor(s) {score['anchors']} "
        f"(band: {band}); amplicon NGS in target cells would confirm/refute."
    )

    hist = score.get("offtarget_mm_histogram", {}) or {}
    low_mm = {m: n for m, n in hist.items() if int(m) <= 2}
    if low_mm:
        testable.append(
            f"Highest-likelihood off-targets: {sum(low_mm.values())} site(s) at <=2 mismatches "
            f"({low_mm}); CHANGE-seq-BE / GUIDE-seq would confirm/refute these specifically."
        )
    testable.append(
        f"Off-target burden is likelihood-weighted; {hist.get('4', 0)} of the nominations are 4-mm "
        f"(low-likelihood); empirical assay needed before treating any as real."
    )

    cons_sev = score.get("terms", {}).get("consequence_severity_nonproductive_bystanders")
    if cons_sev:
        testable.append(
            f"Predicted regulatory consequence of non-productive bystanders (AlphaGenome "
            f"severity-weighted {round(cons_sev, 3)}); a functional/reporter assay would confirm/refute."
        )

    if score.get("ch_gate"):
        testable.append(
            f"CH-driver off-target in {score['ch_genes']} GATES this design (predicted, low-mismatch, "
            f"exonic); requires empirical off-target sequencing + functional confirmation."
        )
    elif score.get("ch_flag_noted"):
        testable.append(
            f"CH-driver off-target noted ({score['ch_genes']}) but NOT gating "
            f"({score.get('ch_flag_reason')}); confirm coordinate, mismatch, and functional effect empirically."
        )

    # Confidence is qualitative and conservative, not a fabricated probability.
    return {
        "overall": OVERALL,
        "source_caveats": SOURCE_CAVEATS,
        "testable_predictions": testable,
        "confidence_basis": (
            "All terms are model predictions or reference-only nominations: not a measured rate, "
            "not a confirmed off-target, and not a per-individual outcome estimate."
        ),
    }


# --- Divergence vs naive baselines ----------------------------------------------------------------

def _ranked_ids(scores: list[GuideScore], key, reverse: bool) -> list[str]:
    ordered = sorted(scores, key=lambda s: (key(s), s.guide_id), reverse=reverse)
    # stable, deterministic ranks (1-based); ties broken by guide_id
    return [s.guide_id for s in ordered]


def _rationale(s: GuideScore, eff_delta: int) -> str:
    """One-line reason the rubric rank differs from the naive-efficiency rank."""
    bits = []
    if s.ch_gate:
        bits.append(f"CH-driver off-target GATES (likely+exonic, {','.join(s.ch_genes)})")
    elif s.ch_flag_noted:
        bits.append(f"CH-driver off-target noted ({','.join(s.ch_genes)}, not gating)")
    if s.tier_label == "BELOW_EVIDENCED_ACTIVITY":
        bits.append("below clinical efficacy anchor")
    elif s.tier_label == "PRECLINICAL_ACTIVITY_ONLY":
        bits.append("preclinical-only activity")
    hist = s.offtarget_mm_histogram or {}
    low_mm = sum(n for m, n in hist.items() if int(m) <= 2)
    if low_mm:
        bits.append(f"{low_mm} off-target(s) at <=2 mismatches (likelihood-weighted)")
    if s.bystander_penalty > 0:
        adv = s.bystander_breakdown.get("adverse", 0)
        if adv:
            bits.append("adverse (motif-destroying) bystander")
    direction = "down-ranked for safety" if eff_delta > 0 else ("up-ranked" if eff_delta < 0 else "unchanged")
    return f"{direction}: " + ("; ".join(bits) if bits else "safety terms reorder vs efficiency")


def divergence_report(records: list[dict], weights: RubricWeights | None = None) -> dict:
    """Compare rubric ranking to naive-efficiency and naive-flat-off-target-count rankings."""
    scores = score_guides(records, weights or RubricWeights())
    rubric_rank = {s.guide_id: s.rank for s in scores}

    eff_order = _ranked_ids(scores, key=lambda s: s.on_target_efficiency, reverse=True)
    eff_rank = {gid: i + 1 for i, gid in enumerate(eff_order)}

    # naive flat count: fewer raw off-targets ranks better (ascending)
    cnt_order = _ranked_ids(scores, key=lambda s: s.offtarget_count, reverse=False)
    cnt_rank = {gid: i + 1 for i, gid in enumerate(cnt_order)}

    rows = []
    for s in scores:
        # positive delta = rubric ranks it WORSE (down-ranked) than the naive method would
        eff_delta = rubric_rank[s.guide_id] - eff_rank[s.guide_id]
        cnt_delta = rubric_rank[s.guide_id] - cnt_rank[s.guide_id]
        rows.append({
            "guide_id": s.guide_id,
            "rubric_rank": rubric_rank[s.guide_id],
            "naive_efficiency_rank": eff_rank[s.guide_id],
            "naive_offtarget_count_rank": cnt_rank[s.guide_id],
            "delta_vs_naive_efficiency": eff_delta,
            "delta_vs_naive_count": cnt_delta,
            "tier": s.tier_label,
            "on_target_efficiency": s.on_target_efficiency,
            "offtarget_count": s.offtarget_count,
            "offtarget_penalty": s.terms["offtarget_penalty"],
            "rationale": _rationale(s, eff_delta),
        })

    changed_eff = [r for r in rows if r["delta_vs_naive_efficiency"] != 0]
    changed_cnt = [r for r in rows if r["delta_vs_naive_count"] != 0]
    movers = sorted(rows, key=lambda r: abs(r["delta_vs_naive_efficiency"]), reverse=True)[:6]

    return {
        "n_designs": len(rows),
        "n_changed_vs_naive_efficiency": len(changed_eff),
        "n_changed_vs_naive_offtarget_count": len(changed_cnt),
        "biggest_movers_vs_naive_efficiency": movers,
        "per_design": rows,
        "note": ("Naive baselines rank by a single dimension; the rubric encodes safety-dominant "
                 "tiers, functional bystanders, likelihood-weighted off-targets, and a conditioned "
                 "CH gate. Divergence demonstrates the rubric adds judgment a naive sort does not."),
    }
