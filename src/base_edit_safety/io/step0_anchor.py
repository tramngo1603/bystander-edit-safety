"""Step 0 anchor: locked HBG1/HBG2 promoter coordinates and in-scope edit positions.

Copied exactly from docs/step0-coordinates.md (LOCKED, verified). Do not re-derive coordinates
here; this module only structures the verified table and emits it as EditRecords. It carries both
allele orientations (genomic plus-strand and promoter/sense-strand). On load, it checks the
4,924 bp HBG1/HBG2 homology invariant and stops the build if it fails (a hard build-stop).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

from ..pipeline.model import EditRecord

GENOME_BUILD = "GRCh38"
CHROM = "chr11"
GENE_STRAND = "-"  # HBG1 and HBG2 are on the minus strand
HBG1_TSS = 5_249_857
HBG2_TSS = 5_254_781
HOMOLOGY_OFFSET = 4924  # HBG2_TSS - HBG1_TSS; hard build-stop invariant


@dataclass(frozen=True)
class AnchorPosition:
    """One in-scope promoter position, with an allele definition shared across both paralogs.

    HBG1 and HBG2 share the same alleles because their proximal promoters are homologous. The
    evidentiary tier (variant_db_status) and dbSNP id differ per gene, so it stores them per gene.
    """

    promoter_offset: int                 # literature numbering vs cap site (positive magnitude)
    editor_class: Literal["ABE", "CBE"]
    genomic_plus_ref: str
    genomic_plus_alt: str
    promoter_sense_ref: str
    promoter_sense_alt: str
    mechanism: str
    hbg1_status: str
    hbg2_status: str
    hbg1_rsid: Optional[str] = None
    hbg2_rsid: Optional[str] = None

    @property
    def hbg1_pos(self) -> int:
        return HBG1_TSS + self.promoter_offset

    @property
    def hbg2_pos(self) -> int:
        return HBG2_TSS + self.promoter_offset


# Evidentiary tiers, per docs/step0-coordinates.md. rsIDs are gene-specific and never
# cross-assigned between HBG1 and HBG2.
DB = "DB-CONFIRMED"
UNVERIFIED = "variant_db_UNVERIFIED"
SYNTHETIC = "SYNTHETIC"

IN_SCOPE: list[AnchorPosition] = [
    AnchorPosition(198, "ABE", "A", "G", "T", "C", "de novo KLF1 motif",
                   DB, UNVERIFIED, hbg1_rsid="rs35710727"),
    AnchorPosition(175, "ABE", "A", "G", "T", "C", "TAL1 E-box (looping)",
                   UNVERIFIED, DB, hbg2_rsid="rs63750654"),
    AnchorPosition(124, "ABE", "A", "G", "T", "C", "synthetic screen hit; KLF1 (double-mutant probe)",
                   SYNTHETIC, SYNTHETIC),
    AnchorPosition(123, "ABE", "A", "G", "T", "C", "synthetic screen hit; KLF1 (double-mutant probe)",
                   SYNTHETIC, SYNTHETIC),
    AnchorPosition(117, "CBE", "C", "T", "G", "A", "distal CCAAT / BCL11A motif (tiled/motif-level)",
                   DB, UNVERIFIED, hbg1_rsid="rs35378915"),
    AnchorPosition(114, "CBE", "G", "A", "C", "T", "distal BCL11A motif (tiled/motif-level)",
                   DB, DB, hbg1_rsid="rs281860601", hbg2_rsid="rs34809449"),
    AnchorPosition(113, "ABE", "T", "C", "A", "G", "de novo GATA1 motif",
                   UNVERIFIED, DB, hbg2_rsid=None),
]
# Note: HBG1 -113 is UNVERIFIED (experimentally demonstrated, no HBG1-specific DB record);
# HBG2 -113 has no listed DB record either and is therefore treated as UNVERIFIED below.


def check_homology_invariant() -> None:
    """Hard build-stop: equivalent HBG1/HBG2 offsets must differ by exactly 4,924 bp."""
    if HBG2_TSS - HBG1_TSS != HOMOLOGY_OFFSET:
        raise AssertionError(
            f"TSS separation {HBG2_TSS - HBG1_TSS} != {HOMOLOGY_OFFSET}; build halted (G8)."
        )
    for pos in IN_SCOPE:
        if pos.hbg2_pos - pos.hbg1_pos != HOMOLOGY_OFFSET:
            raise AssertionError(
                f"offset {pos.promoter_offset}: {pos.hbg2_pos - pos.hbg1_pos} != "
                f"{HOMOLOGY_OFFSET}; build halted (G8)."
            )


def _record_for_gene(pos: AnchorPosition, gene: str) -> EditRecord:
    if gene == "HBG1":
        coord, status, rsid = pos.hbg1_pos, pos.hbg1_status, pos.hbg1_rsid
    elif gene == "HBG2":
        coord = pos.hbg2_pos
        # HBG2 -113 has no listed gene-specific DB record; keep it UNVERIFIED, never DB-confirmed.
        status = pos.hbg2_status if pos.hbg2_status else UNVERIFIED
        rsid = pos.hbg2_rsid
    else:
        raise ValueError(gene)
    provenance = {
        "source": "step0-coordinates.md (locked, verified)",
        "promoter_offset_numbering": "vs transcription cap site (+1 first transcribed base, no zero)",
        "mechanism": pos.mechanism,
    }
    if rsid:
        provenance["dbsnp"] = rsid
    return EditRecord(
        genome_build=GENOME_BUILD,
        chrom=CHROM,
        pos_1based=coord,
        genomic_plus_ref=pos.genomic_plus_ref,
        genomic_plus_alt=pos.genomic_plus_alt,
        promoter_sense_ref=pos.promoter_sense_ref,
        promoter_sense_alt=pos.promoter_sense_alt,
        gene=gene,
        gene_strand=GENE_STRAND,
        promoter_offset=pos.promoter_offset,
        editor_class=pos.editor_class,
        on_target_or_bystander="on_target",
        variant_db_status=status,
        provenance=provenance,
    )


def load_in_scope_records() -> list[EditRecord]:
    """Emit the locked in-scope positions as EditRecords for both paralogs.

    Checks the 4,924 bp homology invariant first; any failure halts the build (G8).
    """
    check_homology_invariant()
    records: list[EditRecord] = []
    for pos in IN_SCOPE:
        records.append(_record_for_gene(pos, "HBG1"))
        records.append(_record_for_gene(pos, "HBG2"))
    return records
