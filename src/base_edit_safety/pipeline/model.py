"""Core data model that flows through every pipeline step.

The types here move in order: a Guide (or an editor-agnostic CandidateGuide) goes to a predictor,
which returns an EditSpectrum, which reduces to one EditRecord per edited position. Both allele
orientations are always carried (see docs/step0-coordinates.md): a single ref_base is unsafe because
the literature uses promoter/sense-strand alleles while genome files use the genomic plus-strand.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Literal

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(base: str) -> str:
    """Watson-Crick complement of a base (case-preserving)."""
    return base.translate(_COMPLEMENT)


def reverse_complement(seq: str) -> str:
    """Reverse complement of a sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


# Canonical base-editing window within the protospacer (1-based, PAM-distal numbering).
EDIT_WINDOW = range(3, 11)
EDIT_SOURCE_BASE = {"ABE": "A", "CBE": "C"}
EDIT_TARGET_BASE = {"ABE": "G", "CBE": "T"}


@dataclass
class EditRecord:
    """One predicted edit at one genomic position (the pipeline's shared output row)."""

    genome_build: str
    chrom: str
    pos_1based: int
    genomic_plus_ref: str
    genomic_plus_alt: str
    promoter_sense_ref: str
    promoter_sense_alt: str
    gene: str                      # "HBG1" | "HBG2" | "HBG1/HBG2_AMBIGUOUS"
    gene_strand: str               # "-"
    promoter_offset: int
    editor_class: Literal["ABE", "CBE"]
    on_target_or_bystander: Literal["on_target", "bystander"]
    predicted_rate: Optional[float] = None
    model_name: Optional[str] = None
    confidence: Optional[str] = None   # qualitative label: model_predicted | reference_nomination | proxy_consequence | deferred
    consequence: Optional[str] = None
    ch_driver_flag: bool = False
    ch_driver_gene: Optional[str] = None
    variant_db_status: str = "unset"   # DB-CONFIRMED | variant_db_UNVERIFIED | SYNTHETIC
    provenance: dict = field(default_factory=dict)
    # NOTE: no clonal-expansion / lifetime-risk field by design (annotation, not prediction).


@dataclass(frozen=True)
class Guide:
    """A protospacer placed on the genome, with strand-aware coordinate mapping.

    Positions are numbered 1-based from the PAM-distal 5' end (the convention base-editing outcome
    models use). Mapping a position to a genomic coordinate needs the strand the protospacer lies on.
    """

    guide_id: str
    protospacer: str                       # 20 nt, 5'->3' on `strand`
    pam: str                               # immediately 3' of protospacer on `strand`
    strand: Literal["+", "-"]
    genome_build: str
    chrom: str
    gene: str                              # explicit; never inferred from homology offset (G8)
    gene_strand: Literal["+", "-"]
    tss_1based: int                        # transcription start site of `gene`
    editor_class: Literal["ABE", "CBE"]
    editor_name: str                       # upstream editor identifier (e.g. ABE8e, BE4max)
    edit_from: str                         # edited base on the protospacer strand (A for ABE, C for CBE)
    edit_to: str
    proto_5p_genomic: int = field(default=0)  # genomic coord of protospacer position 1

    def genomic_coord(self, proto_pos_1based: int) -> int:
        """Genomic coordinate of a 1-based protospacer position."""
        if self.strand == "+":
            return self.proto_5p_genomic + (proto_pos_1based - 1)
        return self.proto_5p_genomic - (proto_pos_1based - 1)

    def promoter_offset(self, proto_pos_1based: int) -> int:
        """Promoter offset (vs cap site) of a protospacer position, given the gene TSS."""
        return self.genomic_coord(proto_pos_1based) - self.tss_1based

    def editable_positions(self) -> list[int]:
        """1-based protospacer positions whose base equals the editor's source base."""
        return [i + 1 for i, b in enumerate(self.protospacer) if b == self.edit_from]

    def alleles_for(self, proto_pos_1based: int) -> dict[str, str]:
        """Return genomic-plus and promoter-sense ref/alt for an edit at a protospacer position.

        The edit is `edit_from`->`edit_to` on the protospacer strand. Genomic-plus alleles follow
        from the protospacer strand; promoter/sense alleles follow from the gene strand.
        """
        if self.strand == "+":
            g_ref, g_alt = self.edit_from, self.edit_to
        else:
            g_ref, g_alt = complement(self.edit_from), complement(self.edit_to)
        # Sense = the gene's own strand. For a minus-strand gene the sense base is the
        # complement of the genomic-plus base.
        if self.gene_strand == "+":
            s_ref, s_alt = g_ref, g_alt
        else:
            s_ref, s_alt = complement(g_ref), complement(g_alt)
        return {
            "genomic_plus_ref": g_ref,
            "genomic_plus_alt": g_alt,
            "promoter_sense_ref": s_ref,
            "promoter_sense_alt": s_alt,
        }


@dataclass(frozen=True)
class CandidateGuide:
    """An enumerated, editor-agnostic guide with its genomic position and off-target summary.

    The guide-design adapter builds it in Step 2. Call for_editor() to make an editor-specific
    Guide from it before predicting an edit spectrum.
    """

    guide_id: str
    protospacer: str
    pam: str
    strand: str
    genome_build: str
    chrom: str
    gene: str
    gene_strand: str
    tss_1based: int
    proto_5p_genomic: int
    offtarget: dict = field(default_factory=dict)
    offtarget_sites: list = field(default_factory=list)  # transient; not all persisted

    def has_editable_base(self, editor_class: str) -> bool:
        """True if the protospacer carries the editor's source base inside the editing window."""
        source = EDIT_SOURCE_BASE[editor_class]
        return any(self.protospacer[p - 1] == source for p in EDIT_WINDOW if p - 1 < len(self.protospacer))

    def for_editor(self, editor_class: str, editor_name: str) -> Guide:
        """Make an editor-specific Guide for outcome prediction."""
        return Guide(
            guide_id=f"{self.guide_id}_{editor_name}",
            protospacer=self.protospacer,
            pam=self.pam,
            strand=self.strand,
            genome_build=self.genome_build,
            chrom=self.chrom,
            gene=self.gene,
            gene_strand=self.gene_strand,
            tss_1based=self.tss_1based,
            editor_class=editor_class,
            editor_name=editor_name,
            edit_from=EDIT_SOURCE_BASE[editor_class],
            edit_to=EDIT_TARGET_BASE[editor_class],
            proto_5p_genomic=self.proto_5p_genomic,
        )


@dataclass(frozen=True)
class OutcomeHaplotype:
    """One predicted output sequence and its probability."""

    output_seq: str
    pred_score: float


@dataclass
class EditSpectrum:
    guide: Guide
    model_name: str
    haplotypes: list[OutcomeHaplotype]
    provenance: dict = field(default_factory=dict)

    def edited_positions(self, hap: OutcomeHaplotype) -> list[int]:
        """1-based protospacer positions that differ from the guide protospacer in this haplotype."""
        ref = self.guide.protospacer
        out = hap.output_seq
        return [i + 1 for i in range(min(len(ref), len(out))) if ref[i] != out[i]]

    def per_position_rate(self) -> dict[int, float]:
        """Marginal probability that each editable position is edited, summed over haplotypes."""
        rates: dict[int, float] = {p: 0.0 for p in self.guide.editable_positions()}
        for hap in self.haplotypes:
            for pos in self.edited_positions(hap):
                if pos in rates:
                    rates[pos] += hap.pred_score
        return rates


def spectrum_to_records(
    spectrum: EditSpectrum,
    on_target_offset: int | None = None,
    on_target_coords: set[int] | None = None,
    allow_overunity: bool = False,
) -> list[EditRecord]:
    """Reduce a spectrum to one EditRecord per editable position.

    A position is flagged on_target when it matches a known therapeutic anchor: either by
    promoter offset (on_target_offset, used for the single-guide spine) or by genomic coordinate
    (on_target_coords, used when enumerating against the locked in-scope anchor set). All other
    edited positions are bystanders. A per-position marginal rate above 1 means the haplotype
    probabilities are not a valid distribution; that raises unless allow_overunity is set.
    """
    guide = spectrum.guide
    rates = spectrum.per_position_rate()
    # Plain qualitative confidence label for these records (model-predicted edit rates); the predictor
    # and ensemble run count live in provenance. Vocabulary is a fixed set of evidence categories.
    confidence = "model_predicted"
    records: list[EditRecord] = []
    for proto_pos, rate in sorted(rates.items()):
        if rate > 1.0 + 1e-6 and not allow_overunity:
            raise ValueError(
                f"per-position edit rate {rate:.4f} > 1 for {guide.guide_id} at protospacer position "
                f"{proto_pos}: haplotype probabilities are not a valid distribution. Fix the predictor "
                f"output or pass allow_overunity=True to keep it as a flagged value."
            )
        offset = guide.promoter_offset(proto_pos)
        coord = guide.genomic_coord(proto_pos)
        if on_target_coords is not None:
            is_on_target = coord in on_target_coords
        else:
            is_on_target = on_target_offset is not None and abs(offset) == abs(on_target_offset)
        # Emit positions the model actually edits, plus the intended on-target position.
        if rate <= 0 and not is_on_target:
            continue
        alleles = guide.alleles_for(proto_pos)
        role = "on_target" if is_on_target else "bystander"
        provenance = {
            "guide_id": guide.guide_id,
            "protospacer": guide.protospacer,
            "pam": guide.pam,
            "strand": guide.strand,
            "protospacer_position": proto_pos,
            "editor": guide.editor_name,
            "predictor_runs": spectrum.provenance.get("ensemble_runs"),
        }
        records.append(
            EditRecord(
                genome_build=guide.genome_build,
                chrom=guide.chrom,
                pos_1based=coord,
                genomic_plus_ref=alleles["genomic_plus_ref"],
                genomic_plus_alt=alleles["genomic_plus_alt"],
                promoter_sense_ref=alleles["promoter_sense_ref"],
                promoter_sense_alt=alleles["promoter_sense_alt"],
                gene=guide.gene,
                gene_strand=guide.gene_strand,
                promoter_offset=offset,
                editor_class=guide.editor_class,
                on_target_or_bystander=role,
                predicted_rate=round(rate, 6),
                model_name=spectrum.model_name,
                confidence=confidence,
                variant_db_status="unset",
                provenance=provenance,
            )
        )
    return records
