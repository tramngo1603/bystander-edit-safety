"""Write pipeline records to disk with the standard claim/disclaim framing."""
from __future__ import annotations

import dataclasses
import json
import os

from ..pipeline.model import EditRecord

DISCLAIMER = (
    "PREDICTED SHORTLIST FOR WET-LAB CONFIRMATION (e.g. CHANGE-seq-BE / amplicon NGS). "
    "Not a safety clearance. No edit is asserted to be safe or off-target; only that models "
    "predict/nominate it. Each prediction carries provenance and confidence."
)

# Full claim/disclaim framing carried on credibility-core artifacts.
CLAIM_DISCLAIM = {
    "what_this_is": (
        "A PREDICTED, LIKELIHOOD-WEIGHTED SAFETY SHORTLIST of base-editing guide/editor designs for "
        "wet-lab confirmation."
    ),
    "what_this_is_not": [
        "NOT a safety clearance: no design is asserted to be safe.",
        "NOT an efficacy ranking: ordering is safety-dominant; efficacy is a threshold, not the goal.",
        "NOT an individual-outcome prediction: the clonal-hematopoiesis layer is driver-gene-membership "
        "annotation/gate only; it produces no probability, time-to-event, or severity value.",
        "NOT measured: on-target rates are model-predicted and off-targets are reference-only "
        "nominations, none empirically confirmed.",
    ],
    "thresholds": (
        "Efficacy and weighting thresholds are PROVISIONAL calibration anchors (allelic-editing "
        "literature), to be re-estimated with the final guide/editor/process and a real potency assay; "
        "the flagged policy knobs are for clinician/expert sign-off."
    ),
    "confirmation_assays": "amplicon NGS (on-target), CHANGE-seq-BE / GUIDE-seq (off-target).",
}

OUTPUT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data",
    "outputs",
)


def output_path(filename: str) -> str:
    return os.path.join(OUTPUT_DIR, filename)


def enumeration_path() -> str:
    """Path to the Steps 2/4 enumeration artifact that Steps 3, 5, and 6 read back."""
    return output_path("enumeration_HBG1_HBG2.json")


def records_to_dicts(records: list[EditRecord]) -> list[dict]:
    return [dataclasses.asdict(r) for r in records]


def _write_json(filename: str, body: dict, meta: dict | None) -> str:
    """Write {disclaimer, claim_disclaim, meta, **body} as indented JSON under data/outputs."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    path = output_path(filename)
    payload = {"disclaimer": DISCLAIMER, "claim_disclaim": CLAIM_DISCLAIM, "meta": meta or {}, **body}
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)
    return path


def write_records(records: list[EditRecord], filename: str, meta: dict | None = None) -> str:
    """Write records as JSON under data/outputs with the disclaimer header. Returns the path."""
    return _write_json(filename, {"record_count": len(records), "records": records_to_dicts(records)}, meta)


def write_credibility(payload_body: dict, filename: str, meta: dict | None = None) -> str:
    """Write the Step 6 credibility-core artifact with the full claim/disclaim framing."""
    return _write_json(filename, payload_body, meta)


def write_ranking(ranked: list[dict], filename: str, meta: dict | None = None) -> str:
    """Write a ranked guide shortlist as JSON under data/outputs with the disclaimer header."""
    return _write_json(filename, {"guide_count": len(ranked), "ranking": ranked}, meta)


def write_design_records(record_set: dict, filename: str, meta: dict | None = None) -> str:
    """Write the HBG design-history artifact set with the claim/disclaim framing."""
    return _write_json(filename, record_set, meta)
