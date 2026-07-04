"""Tests for the HBG design-history artifact: evidence tiers, schema, off-target taxonomy, paralog
hazard, assay checklist, CH consequence class, no clinical-risk fields, reproducibility, score demotion.
"""
from base_edit_safety import annotation, evidence, hbg_registry
from base_edit_safety.design_record import build_design_records

_FORBIDDEN = ("clonal", "expansion", "leukemia", "leukaemia", "mds", "lifetime", "risk",
              "trajectory", "incidence")


def _rec(guide_id, gene, offset, role, rate, *, ref="T", alt="C", sref="A", salt="G",
         pos=5_249_970, editor="ABE", offtarget=None, ch=None, co_target=None):
    prov = {"guide_id": guide_id, "editor": "ABE8e", "protospacer": "CTTGACCAATAGCCTTGACA",
            "pam": "AGG", "strand": "-", "predictor_runs": 5}
    if offtarget is not None:
        prov["offtarget"] = offtarget
    if ch is not None:
        prov["ch_driver_offtargets"] = ch
    if co_target is not None:
        prov["co_targets_paralog"] = True
        prov["co_target_gene"] = co_target
    return {"genome_build": "GRCh38", "chrom": "chr11", "pos_1based": pos,
            "genomic_plus_ref": ref, "genomic_plus_alt": alt,
            "promoter_sense_ref": sref, "promoter_sense_alt": salt,
            "gene": gene, "gene_strand": "-", "promoter_offset": offset,
            "editor_class": editor, "on_target_or_bystander": role, "predicted_rate": rate,
            "model_name": "BE-DICT-bystander:ABE8e", "consequence": "deferred", "provenance": prov}


def _records():
    offt = {"source": "CRISPOR", "offtarget_count": "10", "mit_specificity": "40",
            "mm_histogram": {"2": 1, "4": 5}, "cfd_weighted_burden": 0.3, "paralog_count": 1}
    ch = [{"ch_driver_gene": "TET2", "chrom": "chr4", "start": "105275000",
           "mismatch_count": "1", "feature": "exon"}]
    return [
        _rec("HBG1_m113_ABE", "HBG1", 113, "on_target", 0.45, offtarget=offt, ch=ch, co_target="HBG2"),
        _rec("HBG1_m113_ABE", "HBG1", 116, "bystander", 0.10, pos=5_249_973, offtarget=offt,
             ch=ch, co_target="HBG2"),
    ]


def test_evidence_tier_enforcement():
    rec = build_design_records(_records())["records"][0]
    rm = rec["registry_match"]
    assert rm["relationship"] == hbg_registry.EXACT_EDIT                      # -113 ABE is a known HPFH edit
    assert rm["evidence_strength"] == "published_hspc"                        # not flattened to a weaker match
    assert rm["measured_in"] == "literature"                                 # published, not measured here
    assert rec["evidence_tier"] == evidence.ANALOG_SUPPORTED                 # registry, but no in-run empirical
    rs2 = build_design_records(_records(),
                               empirical={"HBG1_m113_ABE": [{"assay": "amplicon", "edit_rate": 0.4}]})
    assert rs2["records"][0]["evidence_tier"] == evidence.EMPIRICALLY_MEASURED
    assert rs2["records"][0]["registry_match"]["measured_in"] == "this_run"


def test_missing_empirical_never_becomes_model_value():
    rec = build_design_records(_records(), empirical={})["records"][0]
    assert rec["evidence_tier"] != evidence.EMPIRICALLY_MEASURED


def test_predicted_only_when_no_registry_and_no_empirical():
    rec = build_design_records([_rec("g", "HBG1", 999, "on_target", 0.4)])["records"][0]
    assert rec["evidence_tier"] == evidence.PREDICTED_ONLY


def test_design_record_schema():
    rec = build_design_records(_records())["records"][0]
    required = {"record_type", "schema_version", "candidate", "evidence_tier", "registry_match",
                "predicted_edit_spectrum", "off_target", "hbg_paralog_hazard", "ch_annotation",
                "empirical_support", "required_validation", "provisional_scoring", "provenance",
                "limitations"}
    assert required <= set(rec)
    assert [e["promoter_offset"] for e in rec["candidate"]["target"]["intended_edits"]] == [113]


def test_offtarget_taxonomy_separation():
    tax = build_design_records(_records())["records"][0]["off_target"]["taxonomy"]
    assert set(tax) == {"spacer_homology_cas_style", "base_editor_context",
                        "guide_independent_deaminase", "rna_off_target"}
    assert tax["spacer_homology_cas_style"]["mm_histogram"] == {"2": 1, "4": 5}
    assert tax["base_editor_context"]["status"] == "not_assessed_here"
    assert tax["rna_off_target"]["status"] == "not_assessed_here"


def test_hbg_paralog_hazard_flagged():
    hz = build_design_records(_records())["records"][0]["hbg_paralog_hazard"]
    assert hz["paralog_hazard"] is True
    assert any("co-targets" in f for f in hz["flags"])
    assert hz["required_confirmation"]


def test_required_validation_checklist():
    validates = {i["validates"] for i in
                 build_design_records(_records())["records"][0]["required_validation"]}
    assert {"on_target_edit_rate", "bystander_spectrum", "base_editor_context", "rna_off_target",
            "hbg_paralog_rearrangement", "ch_driver_offtarget"} <= validates


def test_ch_consequence_class_lof_vs_gof():
    ch = build_design_records(_records())["records"][0]["ch_annotation"]["driver_offtargets"]
    cc = ch[0]["consequence_class"]
    assert cc["driver_mechanism"] == "lof_tumor_suppressor"           # TET2
    assert cc["consequence_relevance"] == "coding_disruption_relevant"  # exon


def test_no_forbidden_clinical_risk_fields():
    def keys(o):
        if isinstance(o, dict):
            for k, v in o.items():
                yield k
                yield from keys(v)
        elif isinstance(o, list):
            for x in o:
                yield from keys(x)
    for key in keys(build_design_records(_records())):
        assert not any(bad in key.lower() for bad in _FORBIDDEN), key


def test_reproducibility_block_present():
    repro = build_design_records(_records())["reproducibility"]
    assert "generated_at" in repro
    assert "commit" in repro["repo"] and "dirty" in repro["repo"]
    assert "rubric_weights" in repro["config"]


def test_safety_score_demoted_to_secondary():
    rec = build_design_records(_records())["records"][0]
    assert "safety_score" not in rec                       # not a top-level primary field
    ps = rec["provisional_scoring"]
    assert "PROVISIONAL" in ps["note"] and "provisional_score_label" in ps


def test_multi_anchor_design_carries_all_edits():
    # the -123/-124 synthetic pair: both intended edits carried, both registry-matched (not collapsed).
    recs = [_rec("HBG_pair", "HBG1", 123, "on_target", 0.30, pos=5_249_980),
            _rec("HBG_pair", "HBG1", 124, "on_target", 0.28, pos=5_249_981)]
    rec = build_design_records(recs)["records"][0]
    assert {e["promoter_offset"] for e in rec["candidate"]["target"]["intended_edits"]} == {123, 124}
    assert len(rec["registry_match"]["by_edit"]) == 2
    assert {a["relationship"] for a in rec["registry_match"]["by_edit"]} == {hbg_registry.EXACT_EDIT}
    assert len(rec["predicted_edit_spectrum"]["on_target"]) == 2


def test_ch_mechanism_refinements():
    assert annotation.ch_consequence_class("TP53", "exon")["driver_mechanism"] == \
        "tumor_suppressor_mixed_lof_dominant_negative"
    assert annotation.ch_consequence_class("DNMT3A", "exon")["driver_mechanism"] == \
        "loss_or_hypomorphic_with_R882_hotspot"
    mdm4 = annotation.ch_consequence_class("MDM4", "exon")
    assert mdm4["driver_mechanism"] == "copy_number_or_overexpression_driver"
    assert mdm4["consequence_relevance"] == "snv_consequence_poorly_represented"


def test_no_unverified_nejm_dois_in_registry():
    blob = " ".join(c for e in hbg_registry.REGISTRY for c in e.citations)
    assert "10.1056" not in blob   # unverified NEJM DOIs must not ship in an evidence artifact
