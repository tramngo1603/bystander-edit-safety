"""Step 6: divergence vs naive baselines, conservative uncertainty, claim/disclaim framing."""
from base_edit_safety.scoring.credibility import annotate_uncertainty, divergence_report
from base_edit_safety.scoring.rubric import score_guides
from base_edit_safety.io.outputs import CLAIM_DISCLAIM


def _rec(guide_id, offset, role, rate, hist=None, ch=None):
    hist = hist or {}
    return {
        "gene": "HBG1", "editor_class": "ABE", "promoter_offset": offset,
        "on_target_or_bystander": role, "predicted_rate": rate,
        "provenance": {
            "guide_id": guide_id,
            "offtarget": {"offtarget_count": str(sum(hist.values()) + 1), "mit_specificity": "40",
                          "mm_histogram": hist, "cfd_weighted_burden": 0.0, "paralog_count": 1},
            "ch_driver_offtargets": ch or [],
        },
    }


def _ch(gene, mm, feature):
    return {"ch_driver_gene": gene, "mismatch_count": str(mm), "feature": feature,
            "chrom": "chr2", "start": "25264110"}


def test_divergence_down_ranks_high_efficiency_unsafe_guide():
    # naive efficiency would rank the CH-gating guide #1; the rubric must down-rank it.
    recs = ([_rec("clean", 113, "on_target", 0.4, {"2": 1})]
            + [_rec("dirty_ch", 113, "on_target", 0.5, {}, ch=[_ch("DNMT3A", 1, "exon")])])
    rep = divergence_report(recs)
    row = {r["guide_id"]: r for r in rep["per_design"]}
    assert row["dirty_ch"]["naive_efficiency_rank"] == 1     # most efficient
    assert row["dirty_ch"]["rubric_rank"] > row["clean"]["rubric_rank"]  # rubric down-ranks it
    assert row["dirty_ch"]["delta_vs_naive_efficiency"] > 0  # positive = down-ranked for safety
    assert "CH-driver" in row["dirty_ch"]["rationale"]
    assert rep["n_changed_vs_naive_efficiency"] >= 1


def test_divergence_reports_count_and_movers():
    recs = [_rec("a", 113, "on_target", 0.4, {"4": 500}),
            _rec("b", 175, "on_target", 0.4, {"2": 3})]
    rep = divergence_report(recs)
    assert rep["n_designs"] == 2
    assert len(rep["biggest_movers_vs_naive_efficiency"]) >= 1


def test_uncertainty_is_conservative_and_testable():
    s = score_guides([_rec("g", 113, "on_target", 0.4, {"2": 1, "4": 100})])[0].as_dict()
    u = annotate_uncertainty(s)
    assert set(u["source_caveats"]) == {
        "on_target_edit_rate", "off_target_nominations", "regulatory_consequence",
        "ch_feature_call", "efficacy_floor"}
    assert "MODEL-PREDICTED" in u["source_caveats"]["on_target_edit_rate"]
    assert any("amplicon NGS" in t for t in u["testable_predictions"])
    assert any("<=2 mismatches" in t for t in u["testable_predictions"])
    assert "not a measured rate" in u["confidence_basis"]


def test_uncertainty_flags_ch_noted_guide():
    s = score_guides([_rec("g", 123, "on_target", 0.19, {}, ch=[_ch("DNMT3A", 4, "intron")])])[0].as_dict()
    u = annotate_uncertainty(s)
    assert any("CH-driver off-target noted" in t for t in u["testable_predictions"])


def test_uncertainty_has_no_risk_value_token():
    s = score_guides([_rec("g", 113, "on_target", 0.4, {"2": 1})])[0].as_dict()
    blob = str(annotate_uncertainty(s)).lower().replace("clonal-hematopoiesis", "")
    for bad in ("clonal", "expansion", "leukemia", "lifetime", "risk", "trajectory"):
        assert bad not in blob


def test_claim_disclaim_framing_complete():
    assert "LIKELIHOOD-WEIGHTED SAFETY SHORTLIST" in CLAIM_DISCLAIM["what_this_is"]
    joined = " ".join(CLAIM_DISCLAIM["what_this_is_not"]).lower()
    assert "not a safety clearance" in joined
    assert "not an efficacy ranking" in joined
    assert "not an individual-outcome prediction" in joined
    assert "potency assay" in CLAIM_DISCLAIM["thresholds"]
