"""Step 4: CH-driver annotation checks gene membership only (no risk), and does not overclaim."""
from base_edit_safety.annotation import CH_DRIVER_GENES, annotate_record, ch_driver_offtargets
from base_edit_safety.pipeline.model import EditRecord


def _record(chrom, pos, gene="HBG1"):
    return EditRecord(
        genome_build="GRCh38", chrom=chrom, pos_1based=pos,
        genomic_plus_ref="T", genomic_plus_alt="C",
        promoter_sense_ref="A", promoter_sense_alt="G",
        gene=gene, gene_strand="-", promoter_offset=113,
        editor_class="ABE", on_target_or_bystander="on_target",
    )


def test_driver_panel_is_the_twelve_genes():
    assert CH_DRIVER_GENES == {
        "ASXL1", "CHEK2", "DNMT3A", "GNAS", "IDH2", "MDM4", "PPM1D", "SF3B1", "SRSF2",
        "TET2", "TP53", "U2AF1",
    }


def test_promoter_edit_is_not_flagged():
    r = annotate_record(_record("chr11", 5_249_970))  # HBG1 -113
    assert r.ch_driver_flag is False
    assert r.ch_driver_gene is None
    assert "ch_driver_set" in r.provenance


def test_edit_inside_ch_gene_is_flagged():
    r = annotate_record(_record("chr2", 25_264_110, gene="DNMT3A"))  # within DNMT3A span
    assert r.ch_driver_flag is True
    assert r.ch_driver_gene == "DNMT3A"


def test_offtarget_in_ch_gene_detected():
    sites = [{"chrom": "chr2", "start": "25264110", "strand": "+", "mismatchCount": "4",
              "cfdOfftargetScore": "0.1", "locusDesc": "exon:DNMT3A"}]
    hits = ch_driver_offtargets(sites)
    assert len(hits) == 1
    assert hits[0]["ch_driver_gene"] == "DNMT3A"
    assert hits[0]["boostdm_ch_score"] is None  # scores require the CC-BY-NC download


def test_intergenic_near_ch_gene_is_not_overclaimed():
    # a site NEAR TP53 (outside its span) named "intergenic:TP53-WRAP53" must NOT be flagged
    sites = [{"chrom": "chr17", "start": "7700000", "strand": "+", "mismatchCount": "3",
              "locusDesc": "intergenic:TP53-WRAP53"}]
    assert ch_driver_offtargets(sites) == []


def test_no_risk_value_anywhere_in_ch_outputs():
    sites = [{"chrom": "chr2", "start": "25264110", "locusDesc": "exon:DNMT3A"}]
    hit = ch_driver_offtargets(sites)[0]
    forbidden = ("clonal", "expansion", "leukemia", "lifetime", "risk", "trajectory")
    for key in hit:
        assert not any(bad in key.lower() for bad in forbidden)
