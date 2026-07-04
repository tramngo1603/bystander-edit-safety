"""Fixture-driven parser, validation, and ranking tests.

Exercises the CRISPOR and BE-DICT parsing/validation paths and the ranking rubric against SYNTHETIC,
hand-authored representative outputs (see tests/fixtures/README.md; these are NOT captured real-tool
regression fixtures), with NO external tool installed and no subprocess or network.
"""
import csv
import json
import pathlib

import pytest

from base_edit_safety.adapters.bedict import BedictBystanderPredictor, _validated_haplotypes
from base_edit_safety.adapters.crispor import AmbiguousTargetMapping, CrisporGuideDesigner
from base_edit_safety.config import BedictConfig, CrisporConfig
from base_edit_safety.io.reference import ReferenceWindow
from base_edit_safety.pipeline.model import Guide, spectrum_to_records
from base_edit_safety.scoring.rubric import score_guides

FIXTURES = pathlib.Path(__file__).parent / "fixtures"

# HBG1 proximal-promoter window, GRCh38 plus strand (same window used in test_crispor_mapping).
HBG1_WIN = ReferenceWindow(
    genome_build="GRCh38", chrom="chr11", start_1based=5_249_937, end_1based=5_250_067,
    plus_sequence=(
        "GACTATTGGTCAAGTTTGCCTTGTCAAGGCTATTGGTCAAGGCAAGGCTGGCCAACCCATGG"
        "GTGGAGTTTAGCCAGGGACCGTTTCAGACAGATATTTGCATTGAGATAGTGTGGGGAAGGGGCCCCCAA"
    ),
)

GUIDE_M113 = Guide(
    guide_id="HBG1_m113_ABE", protospacer="CTTGACCAATAGCCTTGACA", pam="AGG", strand="-",
    genome_build="GRCh38", chrom="chr11", gene="HBG1", gene_strand="-", tss_1based=5_249_857,
    editor_class="ABE", editor_name="ABE8e", edit_from="A", edit_to="G", proto_5p_genomic=5_249_977,
)


def _predictor():
    # config supplied so __init__ does not try to read env / load a real predictor
    return BedictBystanderPredictor(config=BedictConfig("x", "y"))


def _designer():
    return CrisporGuideDesigner(config=CrisporConfig("x", "y", "hg38", "z"))


def _rows(name):
    with open(FIXTURES / name) as handle:
        return [dict(r) for r in csv.DictReader(handle, delimiter="\t")]


# --- BE-DICT response parsing + validation (recorded fixture) --------------------------------------

def test_bedict_response_parses_to_expected_records():
    response = json.loads((FIXTURES / "bedict_response_HBG1_m113_ABE.json").read_text())
    spectrum = _predictor().spectrum_from_response(GUIDE_M113, "ABE8e", response)
    records = spectrum_to_records(spectrum, on_target_offset=113)
    by_role = {r.on_target_or_bystander: r for r in records}
    assert by_role["on_target"].promoter_offset == 113
    assert abs(by_role["on_target"].predicted_rate - 0.25) < 1e-9   # position 8: 0.20 + 0.05
    assert abs(by_role["bystander"].predicted_rate - 0.20) < 1e-9   # position 5: 0.15 + 0.05
    assert all(r.confidence == "model_predicted" for r in records)  # populated plain label, not None


def test_bedict_rejects_wrong_length_haplotype():
    with pytest.raises(ValueError, match="length"):
        _validated_haplotypes(GUIDE_M113, [{"output_seq": "ACGT", "pred_score": 0.5}])


def test_bedict_rejects_probability_out_of_range():
    with pytest.raises(ValueError, match="outside"):
        _validated_haplotypes(GUIDE_M113, [{"output_seq": "A" * 20, "pred_score": 1.5}])


def test_bedict_rejects_probability_mass_over_one():
    raw = [{"output_seq": "A" * 20, "pred_score": 0.7}, {"output_seq": "C" * 20, "pred_score": 0.7}]
    with pytest.raises(ValueError, match="sum to"):
        _validated_haplotypes(GUIDE_M113, raw)


# --- CRISPOR row parsing + off-target summary (recorded fixture) -----------------------------------

def test_crispor_guide_row_maps_and_summarizes_offtargets():
    guide_rows = _rows("crispor_guides.tsv")
    sites_by_guide = {}
    for s in _rows("crispor_offtargets.tsv"):
        sites_by_guide.setdefault(s["guideId"], []).append(s)
    row = next(r for r in guide_rows if r["guideId"] == "19rev")
    cand = _designer()._row_to_candidate(row, HBG1_WIN, "HBG1", sites_by_guide)
    assert cand.strand == "-" and cand.proto_5p_genomic == 5_249_977
    assert cand.offtarget["mm_histogram"] == {"2": 1, "4": 2}
    assert cand.offtarget["paralog_count"] == 1               # the 0-mismatch site is the paralog
    assert cand.offtarget["offtarget_count_excl_paralog"] == 3


def test_crispor_ambiguous_target_is_rejected():
    # A target that occurs twice in the window must NOT be silently mapped to the first hit.
    target = "ACGTACGTACGTACGTACGTAGG"  # 20 nt protospacer + AGG PAM
    win = ReferenceWindow("GRCh38", "chr11", 1, 200,
                          plus_sequence="TTTTT" + target + "TTTTT" + target + "TTTTT")
    row = {"guideId": "dup", "targetSeq": target, "mitSpecScore": "1", "cfdSpecScore": "1",
           "offtargetCount": "0", "targetGenomeGeneLocus": "x"}
    with pytest.raises(AmbiguousTargetMapping):
        _designer()._row_to_candidate(row, win, "HBG1")


# --- ranking rubric on a representative record set -------------------------------------------------

def _rec(guide_id, offset, role, rate, hist=None, ch=None):
    return {
        "gene": "HBG1", "editor_class": "ABE", "promoter_offset": offset,
        "on_target_or_bystander": role, "predicted_rate": rate,
        "provenance": {
            "guide_id": guide_id,
            "offtarget": {"offtarget_count": str(sum((hist or {}).values()) + 1),
                          "mit_specificity": "40", "mm_histogram": hist or {},
                          "cfd_weighted_burden": 0.0, "paralog_count": 1},
            "ch_driver_offtargets": ch or [],
        },
    }


def test_ranking_is_tier_first_and_decomposed():
    ch_hit = [{"ch_driver_gene": "DNMT3A", "mismatch_count": "1", "feature": "exon",
               "chrom": "chr2", "start": "25264110"}]
    recs = (
        [_rec("clean", 113, "on_target", 0.45, {"2": 1})]
        + [_rec("dirty_ch", 113, "on_target", 0.50, {}, ch=ch_hit)]   # highest efficiency, but CH-gated
        + [_rec("weak", 175, "on_target", 0.10, {"4": 5})]            # below evidenced activity
    )
    scores = {s.guide_id: s for s in score_guides(recs)}
    assert scores["dirty_ch"].ch_gate is True
    assert "REQUIRES_CONFIRMATION" in scores["dirty_ch"].tier_label
    assert scores["clean"].rank < scores["dirty_ch"].rank            # tier dominates raw efficiency
    assert scores["weak"].tier_label == "BELOW_EVIDENCED_ACTIVITY"
    for s in scores.values():                                        # decomposed evidence is present
        assert {"offtarget_penalty", "bystander_penalty", "efficacy_band"} <= set(s.terms)


def test_ranking_rejects_malformed_histogram():
    with pytest.raises(ValueError, match="malformed off-target histogram"):
        score_guides([_rec("bad", 113, "on_target", 0.4, {"two": 1})])
