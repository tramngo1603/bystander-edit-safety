"""Step 5: functional bystander, efficacy tiers, mismatch-weighted off-targets, conditioned CH gate."""
from base_edit_safety.scoring.rubric import (
    RubricWeights, WEIGHTS_STATUS, bystander_function, score_guides,
)


def _rec(guide_id, offset, role, rate, hist=None, ch=None):
    hist = hist or {}
    return {
        "gene": "HBG1", "editor_class": "ABE", "promoter_offset": offset,
        "on_target_or_bystander": role, "predicted_rate": rate,
        "provenance": {
            "guide_id": guide_id,
            "offtarget": {
                "offtarget_count": str(sum(hist.values()) + 1),  # +1 paralog
                "mit_specificity": "40",
                "mm_histogram": hist,
                "cfd_weighted_burden": 0.0,
                "paralog_count": 1,
            },
            "ch_driver_offtargets": ch or [],
        },
    }


def _ch(gene, mm, feature):
    return {"ch_driver_gene": gene, "mismatch_count": str(mm), "feature": feature,
            "chrom": "chr2", "start": "25264110"}


def test_weights_status_provisional_and_knobs_flagged():
    assert "PROVISIONAL" in WEIGHTS_STATUS and RubricWeights().provisional is True
    knobs = RubricWeights().as_dict()["unsettled_policy_knobs"]
    assert any("off_w_mm" in k for k in knobs)
    assert any("ch_gate_mm_threshold" in k for k in knobs)


def test_functional_bystander_lookup():
    assert bystander_function(112)[0] == "adverse"
    assert bystander_function(116)[0] == "productive"
    assert bystander_function(181)[0] == "adverse"
    assert bystander_function(999)[0] == "unknown"


def test_function_not_count_drives_bystander_penalty():
    adv = score_guides([_rec("adv", 113, "on_target", 0.4, {}),
                        _rec("adv", 112, "bystander", 0.2, {})])[0]
    prod = score_guides([_rec("prod", 113, "on_target", 0.4, {}),
                         _rec("prod", 116, "bystander", 0.2, {})])[0]
    assert adv.bystander_penalty > 0 and prod.bystander_penalty == 0.0


def test_offtarget_penalty_is_mismatch_weighted_not_flat():
    # 200 four-mismatch sites should penalize far less than 2 two-mismatch sites (default weights)
    many_4mm = score_guides([_rec("far", 113, "on_target", 0.4, {"4": 200})])[0]
    few_2mm = score_guides([_rec("near", 113, "on_target", 0.4, {"2": 2})])[0]
    assert few_2mm.terms["offtarget_penalty"] > many_4mm.terms["offtarget_penalty"]
    # the decomposition keeps the raw histogram; the paralog hit is excluded
    assert many_4mm.terms["offtarget_mm_histogram"] == {"4": 200}
    assert many_4mm.terms["offtarget_paralog_excluded"] is True


def test_offtarget_weights_configurable():
    recs = [_rec("g", 113, "on_target", 0.4, {"4": 100})]
    base = score_guides(recs)[0].terms["offtarget_penalty"]
    heavier = score_guides(recs, RubricWeights(off_w_mm4=0.1))[0].terms["offtarget_penalty"]
    assert heavier > base


def test_ch_gate_vetoes_only_lowmm_exonic():
    recs = [_rec("clean", 113, "on_target", 0.4, {}),
            _rec("exonic", 113, "on_target", 0.49, {}, ch=[_ch("DNMT3A", 1, "exon")])]
    by = {s.guide_id: s for s in score_guides(recs)}
    assert by["exonic"].ch_gate is True
    assert "REQUIRES_CONFIRMATION" in by["exonic"].tier_label
    assert by["clean"].rank < by["exonic"].rank


def test_ch_intronic_or_highmm_downgraded_to_noted_flag():
    # 4-mm intronic CH hit (the 61forw/DNMT3A case): noted, NOT gating
    s = score_guides([_rec("g", 123, "on_target", 0.19, {}, ch=[_ch("DNMT3A", 4, "intron")])])[0]
    assert s.ch_gate is False
    assert s.ch_flag_noted is True
    assert "NOT gating" in s.ch_flag_reason
    assert s.tier_label != "CH_DRIVER_FLAGGED__REQUIRES_CONFIRMATION"
    assert s.efficacy_band == "PRECLINICAL_ACTIVITY_ONLY"


def test_ch_gate_threshold_is_tunable():
    rec = [_rec("g", 113, "on_target", 0.4, {}, ch=[_ch("DNMT3A", 3, "exon")])]
    assert score_guides(rec)[0].ch_gate is True                              # default threshold 3
    assert score_guides(rec, RubricWeights(ch_gate_mm_threshold=2))[0].ch_gate is False


def test_efficacy_tiers_and_reward():
    def band(e):
        return score_guides([_rec("g", 113, "on_target", e, {})])[0].efficacy_band
    assert band(0.7) == "STRONG_CLINICAL_RANGE"
    assert band(0.4) == "MEETS_CLINICAL_PRECEDENT"
    assert band(0.2) == "PRECLINICAL_ACTIVITY_ONLY"
    assert band(0.1) == "BELOW_EVIDENCED_ACTIVITY"
    hi = score_guides([_rec("hi", 113, "on_target", 0.55, {}),
                       _rec("lo", 175, "on_target", 0.32, {})])
    assert {s.guide_id: s.rank for s in hi}["hi"] == 1


def test_reproducible_and_no_risk_field():
    recs = [_rec("a", 113, "on_target", 0.4, {"4": 10}),
            _rec("a", 112, "bystander", 0.1, {"4": 10})]
    assert [s.safety_score for s in score_guides(recs)] == [s.safety_score for s in score_guides(recs)]
    d = score_guides(recs)[0].as_dict()
    forbidden = ("clonal", "expansion", "leukemia", "lifetime", "risk", "trajectory")
    for key in d:
        assert not any(bad in key.lower() for bad in forbidden)


def test_guides_without_anchor_excluded():
    assert score_guides([_rec("byst_only", 116, "bystander", 0.2, {})]) == []


def _rec_cons(gid, off, role, rate, sev):
    r = _rec(gid, off, role, rate, {})
    r["consequence"] = {"severity": sev, "classification": "functionally-consequential",
                        "source": "AlphaGenome"}
    return r


def test_consequence_penalizes_nonproductive_bystanders_only():
    # adverse bystander (-112) with regulatory severity -> penalized; productive (-116) -> not.
    adv = score_guides([_rec_cons("a", 113, "on_target", 0.4, 0.7),
                        _rec_cons("a", 112, "bystander", 0.2, 0.7)])[0]
    prod = score_guides([_rec_cons("p", 113, "on_target", 0.4, 0.7),
                         _rec_cons("p", 116, "bystander", 0.2, 0.7)])[0]
    assert adv.terms["consequence_penalty"] > 0
    assert prod.terms["consequence_penalty"] == 0.0


def _ch_fc(mm, classification, feature="intron"):
    return {"ch_driver_gene": "DNMT3A", "mismatch_count": str(mm), "feature": feature,
            "chrom": "chr2", "start": "25264110",
            "ch_functional_consequence": {"classification": classification, "severity": 0.6}}


def test_ch_gate_uses_alphagenome_functional_consequence():
    # low-mm + functionally-consequential -> GATES (even though intronic by coordinate)
    g1 = score_guides([_rec("g1", 113, "on_target", 0.4, {}, ch=[_ch_fc(1, "functionally-consequential")])])[0]
    assert g1.ch_gate is True
    # low-mm + low-consequence -> noted, NOT gating (functional call overrides the gate)
    g2 = score_guides([_rec("g2", 113, "on_target", 0.4, {}, ch=[_ch_fc(1, "low-consequence")])])[0]
    assert g2.ch_gate is False and g2.ch_flag_noted is True
    # high-mm + functionally-consequential -> noted (mismatch threshold still excludes)
    g3 = score_guides([_rec("g3", 113, "on_target", 0.4, {}, ch=[_ch_fc(4, "functionally-consequential")])])[0]
    assert g3.ch_gate is False and g3.ch_flag_noted is True
    assert g3.ch_functional_consequence["classification"] == "functionally-consequential"


def test_proxy_ch_consequence_does_not_silently_gate():
    # A low-mm INTRONIC CH hit whose only functional signal is the APPROXIMATE A>G proxy must NOT gate
    # by default (it falls back to the coordinate feature = intron); it gates only if the proxy is trusted.
    hit = {"ch_driver_gene": "DNMT3A", "mismatch_count": "1", "feature": "intron",
           "chrom": "chr2", "start": "25264110",
           "ch_functional_consequence": {"classification": "functionally-consequential",
                                         "severity": 0.7, "is_proxy": True}}
    rec = [_rec("g", 113, "on_target", 0.4, {}, ch=[hit])]
    assert score_guides(rec)[0].ch_gate is False
    assert score_guides(rec, RubricWeights(ch_gate_trust_proxy_consequence=True))[0].ch_gate is True
