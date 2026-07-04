"""Off-target nomination taxonomy for base editing.

Base-editing off-target risk is several distinct things, and a Cas9-style spacer-homology search
answers only one of them. This module keeps the classes separate and, for the ones this repo does not
assess, states that plainly and names the assay that would. A homology-only search is never presented
as base-editor off-target truth.

Classes:
  - spacer_homology_cas_style: protospacer-similarity nominations from a Cas9-style tool (CRISPOR);
    candidate sites by sequence homology on the reference genome, NOT base-editor off-target truth.
  - base_editor_context: base-editor-specific off-targets (editing window + editor context), measured
    empirically by CHANGE-seq-BE / GUIDE-seq2 / equivalent. Not computed here.
  - guide_independent_deaminase: Cas-independent / guide-independent genomic deaminase activity.
    Not computed here.
  - rna_off_target: transcriptome-wide RNA base editing (C>U / A>I) by the deaminase. Not computed here.
"""
from __future__ import annotations

SPACER_HOMOLOGY = "spacer_homology_cas_style"
BASE_EDITOR_CONTEXT = "base_editor_context"
GUIDE_INDEPENDENT_DEAMINASE = "guide_independent_deaminase"
RNA_OFF_TARGET = "rna_off_target"

NOT_ASSESSED = "not_assessed_here"


def classify_offtargets(offtarget_summary: dict | None) -> dict:
    """Wrap a guide's off-target nominations into the base-editing taxonomy.

    Only spacer_homology_cas_style is populated (from the Cas9-style tool). The base-editor-specific
    classes are explicitly marked not-assessed with the assay that would assess them.
    """
    summary = offtarget_summary or {}
    return {
        SPACER_HOMOLOGY: {
            "source": summary.get("source"),
            "status": "reference_genome_homology_nomination",
            "mm_histogram": summary.get("mm_histogram", {}),
            "cfd_weighted_burden": summary.get("cfd_weighted_burden"),
            "mit_specificity": summary.get("mit_specificity"),
            "caveat": ("Cas9-style protospacer-similarity nominations on the reference genome only "
                       "(no bulges, no population variants); NOT base-editor off-target truth."),
        },
        BASE_EDITOR_CONTEXT: {"status": NOT_ASSESSED,
                              "assay": "CHANGE-seq-BE or GUIDE-seq2 (base-editor context)"},
        GUIDE_INDEPENDENT_DEAMINASE: {"status": NOT_ASSESSED,
                                      "assay": "orthogonal genomic deaminase assay (e.g. Detect-seq) "
                                               "plus a guide-independent editing control"},
        RNA_OFF_TARGET: {"status": NOT_ASSESSED,
                         "assay": "RNA-seq editome (transcriptome-wide C>U / A>I)"},
    }


def assay_triggers(taxonomy: dict) -> list[dict]:
    """Assay items triggered by the off-target taxonomy (the not-assessed classes always trigger)."""
    triggers = []
    for cls in (BASE_EDITOR_CONTEXT, GUIDE_INDEPENDENT_DEAMINASE, RNA_OFF_TARGET):
        entry = taxonomy.get(cls, {})
        if entry.get("status") == NOT_ASSESSED:
            triggers.append({"validates": cls, "assay": entry["assay"], "status": "required"})
    homology = taxonomy.get(SPACER_HOMOLOGY, {})
    if homology.get("mm_histogram"):
        triggers.append({"validates": SPACER_HOMOLOGY,
                         "assay": "empirical off-target sequencing (CHANGE-seq-BE / GUIDE-seq2) at the "
                                  "nominated sites", "status": "required"})
    return triggers
