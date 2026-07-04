"""Curated registry of de-risked HBG1/HBG2 fetal-hemoglobin editing concepts, and candidate matching.

A local, HBG-specific registry of edit concepts with clinical or experimental support: the natural
HPFH promoter lesions already in Step 0, published HBG-promoter base-editing / HSPC evidence, and the
clinical HBG1/HBG2 promoter programs. Each candidate is reported by its structural relationship to the
registry (exact_edit / same_motif / nearby / no_match) AND a separate evidence strength, so an exact
published HSPC edit is not flattened into a weaker match.

Scope notes:
  - Clinical-program exact target positions may be non-public, so those entries are recorded at the
    modality/motif level and are only ever reported as analogs, never as an identity match.
  - Citations here are trial-registry / program-level identifiers. Primary-literature DOIs are NOT
    asserted until independently verified: evidence artifacts should not carry unverifiable citations.
"""
from __future__ import annotations

from dataclasses import dataclass

from .io import step0_anchor as anchor

# How a candidate relates to the registry, structurally.
EXACT_EDIT = "exact_edit"
SAME_MOTIF = "same_motif"
NEARBY = "nearby"
NO_MATCH = "no_match"

# Separate evidence strength (weakest -> strongest). Kept apart from the evidence TIER so that
# published literature never launders into `empirically_measured` (which means measured in this run).
STRENGTH_ORDER = ("none", "analog", "motif", "published_hspc", "clinical")


@dataclass(frozen=True)
class RegistryEntry:
    entry_id: str
    label: str
    genes: tuple
    modality: str                 # base_editing | nuclease
    editor_class: str | None      # ABE | CBE | None (nucleases, or unconfirmed base-editor class)
    offsets: tuple                # positive promoter-offset magnitudes, or () if not position-specific
    motif: str                    # BCL11A | GATA1 | KLF1 | TAL1 | ...
    strength: str                 # one of STRENGTH_ORDER
    citations: tuple
    note: str


# Real citations (PMIDs) for the Step 0 HPFH positions (from docs/step0-coordinates.md).
_CITATIONS = {
    198: ("PMID 2430647 (clinical)", "PMID 28659276 (KLF1 mechanism)", "PMID 37400614 (HSPC ABE)"),
    175: ("PMID 2449926 (clinical HBG2)", "PMID 25971621 (TAL1 mechanism)", "PMID 37400614 (HSPC ABE)"),
    124: ("PMID 35147495 (synthetic screen; KLF1 double-mutant)",),
    123: ("PMID 35147495 (synthetic screen; KLF1 double-mutant)",),
    117: ("PMID 2417646 (clinical)", "PMID 3181130 (distal CCAAT / BCL11A)", "PMID 37989316 (motif BE)"),
    114: ("PMID 1698280 (clinical HBG2)", "PMID 1704803 (clinical HBG1)", "PMID 35147495 / 37989316 (BE)"),
    113: ("PMID 23621512 (clinical)", "PMID 30617196 (de novo GATA1)",
          "PMID 36006707 (exact -113 A>G ABE8e)", "PMID 37400614 (HSPC ABE)"),
}

# Evidence strength per HPFH offset (docs/step0-coordinates.md engineering-anchor tiering).
_STRENGTH = {
    198: "published_hspc",   # HSPC ABE (Mayuranathan 2023)
    175: "published_hspc",   # strongest exact-allele HSPC ABE (Mayuranathan 2023)
    124: "published_hspc",   # synthetic -123/-124 pair, HSPC-derived erythroid (Ravi 2022)
    123: "published_hspc",
    117: "motif",            # tiled / motif-level (Han 2023)
    114: "motif",            # tiled / motif-level (Ravi 2022 / Han 2023)
    113: "published_hspc",   # exact -113 A>G ABE8e in patient CD34+ HSCs (Li 2022)
}


def _motif_for(mechanism: str) -> str:
    m = mechanism.upper()
    for motif in ("BCL11A", "GATA1", "KLF1", "TAL1", "CCAAT"):
        if motif in m:
            return "BCL11A" if motif == "CCAAT" else motif   # distal CCAAT box is the BCL11A site
    return "unknown"


def _hpfh_entries() -> list[RegistryEntry]:
    return [
        RegistryEntry(
            entry_id=f"hpfh_-{pos.promoter_offset}_{pos.editor_class}",
            label=f"HBG -{pos.promoter_offset} {pos.editor_class} ({pos.mechanism})",
            genes=("HBG1", "HBG2"), modality="base_editing", editor_class=pos.editor_class,
            offsets=(pos.promoter_offset,), motif=_motif_for(pos.mechanism),
            strength=_STRENGTH.get(pos.promoter_offset, "motif"),
            citations=_CITATIONS.get(pos.promoter_offset, ()),
            note="Natural HPFH lesion and/or published HSPC base-editing evidence (Step 0 in-scope).",
        )
        for pos in anchor.IN_SCOPE
    ]


# Clinical HBG1/HBG2 promoter programs. Recorded at the modality/motif level; positions may be
# non-public. Base-editor CLASS for BEAM-101 is publicly reported but not asserted here (editor_class
# left None), and primary-literature DOIs are omitted pending verification.
_CLINICAL = [
    RegistryEntry(
        "beam_101", "BEAM-101 (Beam Therapeutics, BEACON)", ("HBG1", "HBG2"),
        "base_editing", None, (), "BCL11A", "clinical",
        ("NCT05456880 (BEACON trial)", "Beam Therapeutics BEAM-101 program materials (public)"),
        "Base-edited CD34+ HSPCs in the HBG1/HBG2 promoter, designed to inhibit BCL11A promoter "
        "binding (HPFH-mimicking): same modality (base editing) and same locus/mechanism as this "
        "tool's scope. Base-editor class (reported as ABE publicly) is not asserted here.",
    ),
    RegistryEntry(
        "reni_cel_edit301", "reni-cel / EDIT-301 (Editas, RUBY / EdiThal)", ("HBG1", "HBG2"),
        "nuclease", None, (), "BCL11A", "clinical",
        ("NCT04853576 (RUBY)", "NCT05444894 (EdiThal)", "Editas reni-cel program materials (8-K summary)"),
        "AsCas12a editing of the HBG1/HBG2 promoter BCL11A site: SAME locus and mechanism, DIFFERENT "
        "modality (nuclease, not base editing), so its edits are indels/deletions rather than point edits.",
    ),
    RegistryEntry(
        "bcl11a_enhancer", "BCL11A erythroid-enhancer editing (e.g. exa-cel)", ("BCL11A",),
        "nuclease", None, (), "BCL11A_enhancer", "clinical",
        ("Casgevy / exa-cel (program-level)",),
        "Cas9 disruption of the BCL11A +58 erythroid enhancer: same THERAPEUTIC GOAL (HbF reactivation) "
        "but a DIFFERENT target (BCL11A enhancer, not the HBG promoter).",
    ),
]

REGISTRY = _hpfh_entries() + _CLINICAL
_BY_ID = {e.entry_id: e for e in REGISTRY}


def entry(entry_id: str) -> RegistryEntry | None:
    return _BY_ID.get(entry_id)


def relationship(gene: str, promoter_offset: int, editor_class: str, mechanism: str | None = None) -> dict:
    """Report a candidate's registry match: the strongest structural relationship, a separate evidence
    strength (over all matched entries), and the per-entry matches.

    exact_edit: a base-editing registry entry at the same promoter offset and editor class.
    same_motif: same motif (e.g. BCL11A) at a registry HBG-promoter entry, not an exact edit.
    nearby: within 10 bp of a registry offset.
    """
    offset = abs(promoter_offset)
    motif = _motif_for(mechanism) if mechanism else None
    matches: list[tuple[str, RegistryEntry]] = []
    for e in REGISTRY:
        if offset in e.offsets and e.modality == "base_editing" and e.editor_class == editor_class:
            matches.append((EXACT_EDIT, e))
        elif motif and motif != "unknown" and motif == e.motif and "HBG1" in e.genes:
            matches.append((SAME_MOTIF, e))
        elif e.offsets and any(abs(o - offset) <= 10 for o in e.offsets):
            matches.append((NEARBY, e))

    levels = {lvl for lvl, _ in matches}
    if EXACT_EDIT in levels:
        rel = EXACT_EDIT
    elif SAME_MOTIF in levels:
        rel = SAME_MOTIF
    elif matches:
        rel = NEARBY
    else:
        rel = NO_MATCH

    strengths = [e.strength for _, e in matches] or ["none"]
    strength = max(strengths, key=STRENGTH_ORDER.index)

    return {
        "relationship": rel,
        "evidence_strength": strength,
        "matched_entries": [{"entry_id": e.entry_id, "label": e.label, "relationship": lvl,
                             "strength": e.strength, "citations": list(e.citations)}
                            for lvl, e in matches],
        "note": "Clinical-program exact positions may be non-public; motif matches are analogs.",
    }
