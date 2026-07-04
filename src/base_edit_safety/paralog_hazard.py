"""HBG1/HBG2 paralog-region hazard layer (a locus-specific assay-planning flag, not a risk score).

HBG1 and HBG2 are a recent duplication: their proximal promoters are nearly identical and the genes
sit ~4.9 kb apart on chr11. Editing in this duplicated region can produce events that generic
off-target tools structurally miss: co-editing of both paralogs, and large deletions / rearrangements
driven by the HBG1-HBG2 repeat (notably the ~4.9 kb inter-HBG deletion). This flags candidates in the
high-homology paralog context and names the confirmation assays they need. It is assay planning, not
a proven clinical risk.
"""
from __future__ import annotations

from .io import step0_anchor as anchor


def paralog_hazard(record: dict) -> dict:
    """Assess HBG1/HBG2 paralog hazard for one edit record. Returns flags + required confirmation."""
    prov = record.get("provenance", {})
    co_targets = bool(prov.get("co_targets_paralog"))
    in_promoter_window = record.get("gene") in ("HBG1", "HBG2")

    flags = []
    if in_promoter_window:
        flags.append("edits the duplicated HBG1/HBG2 proximal promoter (high paralog homology)")
    if co_targets:
        flags.append(f"guide co-targets the paralog ({prov.get('co_target_gene')})")

    hazard = bool(flags)
    assays = []
    if hazard:
        assays = [
            {"validates": "hbg_paralog_rearrangement",
             "assay": "long-read or amplicon sequencing spanning HBG1-HBG2 for the ~4.9 kb inter-HBG "
                      "deletion / rearrangement", "status": "required"},
            {"validates": "hbg_paralog_cotargeting",
             "assay": "paralog-resolving amplicon NGS quantifying editing at HBG1 and HBG2 separately",
             "status": "required"},
        ]
    return {
        "paralog_hazard": hazard,
        "paralog_offset_bp": anchor.HOMOLOGY_OFFSET,   # 4924 bp between the HBG1/HBG2 paralogs
        "flags": flags,
        "required_confirmation": assays,
        "note": "Locus-specific assay-planning flag for the HBG1/HBG2 duplication, not a clinical risk score.",
    }
