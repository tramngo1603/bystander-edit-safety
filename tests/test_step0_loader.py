"""Step 0 loader: reproduces locked in-scope rows with both orientations and the homology gate."""
import pytest

from base_edit_safety.io import step0_anchor as anchor
from base_edit_safety.io.step0_anchor import load_in_scope_records


def test_homology_invariant_holds():
    anchor.check_homology_invariant()  # must not raise


def test_invariant_failure_halts(monkeypatch):
    # A wrong TSS separation must HALT the build (the G8 gate), not pass silently.
    monkeypatch.setattr(anchor, "HBG2_TSS", anchor.HBG1_TSS + anchor.HOMOLOGY_OFFSET + 1)
    with pytest.raises(AssertionError):
        anchor.check_homology_invariant()


def test_loader_emits_both_paralogs_for_every_position():
    records = load_in_scope_records()
    assert len(records) == 2 * len(anchor.IN_SCOPE)
    genes = {r.gene for r in records}
    assert genes == {"HBG1", "HBG2"}


def test_dual_orientation_populated_everywhere():
    for r in load_in_scope_records():
        assert r.genomic_plus_ref and r.genomic_plus_alt
        assert r.promoter_sense_ref and r.promoter_sense_alt
        # sense base is the complement of the genomic-plus base for these minus-strand genes
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        assert r.promoter_sense_ref == comp[r.genomic_plus_ref]
        assert r.promoter_sense_alt == comp[r.genomic_plus_alt]


def test_known_anchor_coordinates_match_locked_doc():
    by_key = {(r.gene, r.promoter_offset): r for r in load_in_scope_records()}
    # HBG2 -114 = 5,254,895 (the corrected row), genomic G>A / sense C>T
    r = by_key[("HBG2", 114)]
    assert r.pos_1based == 5_254_895
    assert (r.genomic_plus_ref, r.genomic_plus_alt) == ("G", "A")
    assert (r.promoter_sense_ref, r.promoter_sense_alt) == ("C", "T")
    # HBG1 -113 = 5,249,970, genomic T>C / sense A>G
    r = by_key[("HBG1", 113)]
    assert r.pos_1based == 5_249_970
    assert (r.genomic_plus_ref, r.genomic_plus_alt) == ("T", "C")
    assert (r.promoter_sense_ref, r.promoter_sense_alt) == ("A", "G")


def test_variant_db_tiers_not_overclaimed():
    by_key = {(r.gene, r.promoter_offset): r for r in load_in_scope_records()}
    # HBG2 -175 is DB-confirmed (rs63750654); HBG1 -175 must NOT be DB-confirmed.
    assert by_key[("HBG2", 175)].variant_db_status == anchor.DB
    assert by_key[("HBG1", 175)].variant_db_status == anchor.UNVERIFIED
    # rsID never cross-assigned between paralogs.
    assert by_key[("HBG2", 175)].provenance.get("dbsnp") == "rs63750654"
    assert "dbsnp" not in by_key[("HBG1", 175)].provenance
    # Synthetic screen hits flagged as such.
    assert by_key[("HBG1", 124)].variant_db_status == anchor.SYNTHETIC
    assert by_key[("HBG1", 123)].variant_db_status == anchor.SYNTHETIC
