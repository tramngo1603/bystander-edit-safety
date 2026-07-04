"""Step 4: static clonal-hematopoiesis driver-gene annotation of edits and off-targets.

ANNOTATION ONLY (G3): for each edit, sets a boolean ch_driver_flag and the driver gene symbol by
checking which gene span the edit falls in, and reports off-target sites that fall in CH-driver genes.
Membership is decided from coordinates only. It never produces a clonal-expansion, leukemia/MDS,
lifetime-risk, or trajectory value.

The driver-gene set is the boostDM-CH panel of 12 genes; spans are GRCh38 (Ensembl). Per-mutation
driver scores are out of scope for this flag-level build. If they are ever wired in, the sources are
boostDM-CH (CC-BY-NC, https://www.intogen.org/boostdm-ch) and, for coding sites, AlphaMissense
(CC-BY-4.0, Zenodo record 10813168, file AlphaMissense_hg38.tsv.gz).
"""
from __future__ import annotations

from typing import Optional

from .pipeline.model import EditRecord

# boostDM-CH driver-gene panel (12 genes; upstream scores are CC-BY-NC).
CH_DRIVER_GENES = frozenset(
    {"ASXL1", "CHEK2", "DNMT3A", "GNAS", "IDH2", "MDM4", "PPM1D", "SF3B1", "SRSF2", "TET2",
     "TP53", "U2AF1"}
)

# GRCh38 gene spans (Ensembl lookup-by-symbol), 1-based inclusive.
CH_GENE_SPANS = {
    "ASXL1": ("chr20", 32358266, 32439319),
    "CHEK2": ("chr22", 28687738, 28742422),
    "DNMT3A": ("chr2", 25227855, 25342590),
    "GNAS": ("chr20", 58839669, 58911192),
    "IDH2": ("chr15", 90083045, 90102511),
    "MDM4": ("chr1", 204516369, 204558120),
    "PPM1D": ("chr17", 60599735, 60666293),
    "SF3B1": ("chr2", 197388515, 197435079),
    "SRSF2": ("chr17", 76734115, 76737919),
    "TET2": ("chr4", 105145875, 105279816),
    "TP53": ("chr17", 7661779, 7687546),
    "U2AF1": ("chr21", 43092956, 43107605),
}
CH_DRIVER_SET_PROVENANCE = "boostDM-CH 12-gene panel (CC-BY-NC); GRCh38 spans from Ensembl"


def _norm_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def gene_at(chrom: str, pos_1based: int) -> Optional[str]:
    """Return the CH-driver gene whose span contains this coordinate, else None."""
    chrom = _norm_chrom(chrom)
    for gene, (gchrom, start, end) in CH_GENE_SPANS.items():
        if chrom == gchrom and start <= pos_1based <= end:
            return gene
    return None


def annotate_record(record: EditRecord) -> EditRecord:
    """Set the flag for whether this edit's own position falls in a CH-driver gene (static)."""
    gene = gene_at(record.chrom, record.pos_1based)
    record.ch_driver_flag = gene is not None
    record.ch_driver_gene = gene
    record.provenance = {
        **record.provenance,
        "ch_driver_set": CH_DRIVER_SET_PROVENANCE,
        "ch_annotation": "static driver-gene-membership flag only; not a predictive estimate",
    }
    return record


def annotate_records(records: list[EditRecord]) -> list[EditRecord]:
    """Apply the per-edit CH-driver flag to every record."""
    return [annotate_record(r) for r in records]


def ch_driver_offtargets(offtarget_sites: list[dict]) -> list[dict]:
    """Keep only the off-target sites that fall in a CH-driver gene (by coordinate).

    Each input site should carry chrom/start (and optionally locusDesc). Returns one short record
    per CH-driver off-target with the driver gene and the tool's mismatch/score fields.
    """
    hits: list[dict] = []
    for site in offtarget_sites:
        chrom = site.get("chrom")
        start = site.get("start")
        if chrom is None or start is None:
            continue
        try:
            pos = int(start)
        except (TypeError, ValueError):
            continue
        # Membership is by coordinate only. The tool's locus name (e.g. "intergenic:TP53-WRAP53")
        # can name a CH gene the site only sits NEAR; using it to claim membership would overclaim,
        # so we record it for provenance but do not use it to set the flag.
        gene = gene_at(chrom, pos)
        if gene is not None:
            locus_desc = site.get("locusDesc", "") or ""
            # CRISPOR locus feature prefix: exon / intron / intergenic (utr is reported as exon).
            feature = locus_desc.split(":", 1)[0] if ":" in locus_desc else "unknown"
            hits.append({
                "ch_driver_gene": gene,
                "chrom": _norm_chrom(chrom),
                "start": pos,
                "strand": site.get("strand"),
                "mismatch_count": site.get("mismatchCount"),
                "cfd_offtarget_score": site.get("cfdOfftargetScore"),
                "locus_desc": locus_desc,
                "feature": feature,                # exon|intron|intergenic (annotation-only; not AlphaGenome)
                "boostdm_ch_score": None,          # populated only for coding sites once scores are present
            })
    return hits


# Consequence logic differs by driver mechanism (documentation/monitoring context only, NOT a risk
# score): most CH drivers are loss-of-function tumor suppressors where a coding-disrupting edit is the
# concern; some act through specific gain-of-function hotspots where only particular residues matter.
CH_DRIVER_MECHANISM = {
    "ASXL1": "lof_tumor_suppressor",
    "CHEK2": "lof_tumor_suppressor",
    "TET2": "lof_tumor_suppressor",
    "DNMT3A": "loss_or_hypomorphic_with_R882_hotspot",     # LoF/hypomorph plus the R882 hotspot
    "TP53": "tumor_suppressor_mixed_lof_dominant_negative",
    "PPM1D": "truncating_activating",       # C-terminal truncations stabilize PPM1D
    "SF3B1": "spliceosome_hotspot",
    "SRSF2": "spliceosome_hotspot",
    "U2AF1": "spliceosome_hotspot",
    "IDH2": "gof_hotspot",                  # R140 / R172
    "GNAS": "gof_hotspot",                  # R201
    "MDM4": "copy_number_or_overexpression_driver",   # driven by CNV/overexpression; SNV/base-edit consequence poorly represented
}


def ch_consequence_class(gene: str, feature: str | None) -> dict:
    """Consequence-aware annotation for a CH-driver off-target (documentation/monitoring only).

    Combines the coordinate feature (exon/intron/utr/intergenic) with the driver's mechanism so a
    coding disruption in a loss-of-function tumor suppressor is distinguished from a gain-of-function
    hotspot gene where only specific residues matter. Produces NO risk, incidence, or clinical estimate.
    """
    mechanism = CH_DRIVER_MECHANISM.get(gene, "unknown")
    feat = (feature or "unknown").lower()
    coding = feat in ("exon", "utr")
    coding_relevant = {"lof_tumor_suppressor", "loss_or_hypomorphic_with_R882_hotspot",
                       "tumor_suppressor_mixed_lof_dominant_negative"}
    hotspot = {"gof_hotspot", "spliceosome_hotspot", "truncating_activating"}
    if mechanism in coding_relevant:
        relevance = "coding_disruption_relevant" if coding else "non_coding_lower_concern"
    elif mechanism in hotspot:
        relevance = "hotspot_specific" if coding else "non_coding_lower_concern"
    elif mechanism == "copy_number_or_overexpression_driver":
        relevance = "snv_consequence_poorly_represented"
    else:
        relevance = "unknown"
    return {
        "gene": gene,
        "feature": feat,
        "driver_mechanism": mechanism,
        "consequence_relevance": relevance,
        "scope": ("documentation/monitoring support only; not a clonal-expansion, leukemia, incidence, "
                  "time-to-event, or clinical-risk estimate (G3)"),
    }
