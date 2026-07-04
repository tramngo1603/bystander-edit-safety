"""Evidence tiers for HBG design candidates.

Every candidate carries one explicit evidence state. A missing measurement never silently becomes a
model value: a candidate is `predicted_only` unless it matches the curated HBG evidence registry
(`analog_supported`) or has real measurements attached (`empirically_measured`).
"""
from __future__ import annotations

PREDICTED_ONLY = "predicted_only"
ANALOG_SUPPORTED = "analog_supported"
EMPIRICALLY_MEASURED = "empirically_measured"
EVIDENCE_TIERS = (PREDICTED_ONLY, ANALOG_SUPPORTED, EMPIRICALLY_MEASURED)

# Registry relationships that count as analog support (see hbg_registry.relationship).
_ANALOG_RELATIONSHIPS = frozenset({"exact_edit", "same_motif", "nearby"})


def evidence_tier(registry_relationship: str | None, empirical_measurements) -> str:
    """Assign the evidence tier for one candidate.

    Empirical measurement wins; otherwise a registry analog; otherwise predicted-only. Missing
    empirical evidence stays predicted-only and is never upgraded to a model value.
    """
    if empirical_measurements:
        return EMPIRICALLY_MEASURED
    if registry_relationship in _ANALOG_RELATIONSHIPS:
        return ANALOG_SUPPORTED
    return PREDICTED_ONLY
