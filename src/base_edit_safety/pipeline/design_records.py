"""Emit the HBG1/HBG2 design-history + validation-planning artifact set (Steps 2/4 -> design records).

This is the primary product artifact: one signable design-rationale record per candidate, foregrounding
evidence vectors and required assays rather than a scalar rank.
"""
from __future__ import annotations

import json

from ..design_record import build_design_records
from ..io import outputs
from ..scoring.rubric import RubricWeights


def run_design_records(record_dicts: list[dict] | None = None,
                       weights: RubricWeights | None = None,
                       empirical: dict | None = None, write: bool = True) -> dict:
    """Build the HBG design-history artifact set. Loads the persisted Steps 2/4 output if none given."""
    weights = weights or RubricWeights()
    if record_dicts is None:
        with open(outputs.enumeration_path()) as handle:
            record_dicts = json.load(handle)["records"]

    record_set = build_design_records(record_dicts, weights, empirical)
    if write:
        outputs.write_design_records(
            record_set, filename="design_records_HBG1_HBG2.json",
            meta={"step": "design-history + validation-planning artifact generation",
                  "count": record_set["count"]},
        )
    return record_set


if __name__ == "__main__":
    rs = run_design_records()
    print(f"HBG design-history records: {rs['count']}\n")
    for r in rs["records"]:
        c = r["candidate"]
        edits = ",".join(f"-{e['promoter_offset']}" for e in c["target"]["intended_edits"])
        print(f"  {c['guide_id']:<22} {c['gene']} {edits:<10} "
              f"tier={r['evidence_tier']:<18} registry={r['registry_match']['relationship']:<12} "
              f"paralog_hazard={r['hbg_paralog_hazard']['paralog_hazard']!s:<5} "
              f"assays={len(r['required_validation'])}")
