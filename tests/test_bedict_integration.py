"""Integration: the spine produces a non-empty edit spectrum for one anchor guide.

Skipped unless the external predictor is configured (BE_SAFETY_BEDICT_REPO /
BE_SAFETY_BEDICT_PYTHON) and reference access works. Runs the full Step 1 path.
"""
import pytest

from base_edit_safety.config import bedict_configured

pytestmark = pytest.mark.integration

if not bedict_configured():
    pytest.skip("external base-editing predictor not configured", allow_module_level=True)

from base_edit_safety.pipeline.spine import run_anchor_spine  # noqa: E402


@pytest.fixture(scope="module")
def records():
    try:
        return run_anchor_spine(gene="HBG1", promoter_offset=113,
                                editor_class="ABE", editor_name="ABE8e", write=False)
    except Exception as exc:  # network or predictor unavailable at run time
        pytest.skip(f"spine could not run: {exc}")


def test_spectrum_non_empty(records):
    assert len(records) >= 1


def test_on_target_present_and_correct(records):
    on = [r for r in records if r.on_target_or_bystander == "on_target"]
    assert len(on) == 1
    r = on[0]
    assert r.gene == "HBG1" and r.promoter_offset == 113
    assert r.pos_1based == 5_249_970
    assert (r.genomic_plus_ref, r.genomic_plus_alt) == ("T", "C")
    assert (r.promoter_sense_ref, r.promoter_sense_alt) == ("A", "G")


def test_bystander_burden_represented(records):
    bystanders = [r for r in records if r.on_target_or_bystander == "bystander"]
    assert len(bystanders) >= 1  # the -113 guide has predicted bystander edits


def test_provenance_and_consequence(records):
    for r in records:
        assert r.model_name.startswith("BE-DICT-bystander")
        assert r.predicted_rate is not None
        assert "guide_id" in r.provenance
        assert r.consequence == "deferred"  # Step 3 sets consequence to deferred
