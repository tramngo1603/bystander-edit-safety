"""Safety-ranking rubric (Step 5). Structure approved; WEIGHTS ARE PROVISIONAL (calibrated starts).

The rubric structure (CH gate, efficacy threshold, decomposition, active consequence slot) is approved.
The numbers below are evidence-anchored calibrated starting points, not final sign-off; several are
still unsettled and flagged as policy/clinical judgment. provisional stays True until a human reviews
the re-run. All knobs live in RubricWeights, so you can retune them without changing the pipeline.

This ranks by SAFETY, not efficacy. The design:

  1. Penalize bystanders by FUNCTION, not count. Mayuranathan 2023 (PMID 37400614) shows a
     higher-bystander editor (ABE8e) can be the better therapeutic, so "fewer bystanders" does not
     mean "safer." Each bystander is scored by what it does at the motif it hits (adverse /
     productive / neutral / unknown), using FUNCTIONAL_BYSTANDER below. EMERGING evidence; see
     per-entry sources.
  2. TIERED efficacy floor anchored to clinical/preclinical editing data (allelic-editing anchors;
     see EFFICACY_CAVEAT). Above the floor, a mild efficiency reward stops a barely-clearing guide
     from winning just for being empty.
  3. The CH gate is the dominating partition, framed as a NOMINATION, not a verdict: a predicted
     CH-driver off-target drops a guide to the worst tier "requires empirical confirmation" (per FDA
     guidance, a gene-name match alone is not enough). A functional-consequence field carries the
     AlphaGenome regulatory severity that Step 3 (annotate_consequence.py) fills in; it is empty until that
     step runs, so the term is zero by default. No clonal/lifetime-risk value is produced (G3).
  4. Fully decomposable: every guide's rank breaks into the terms that built it.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

WEIGHTS_STATUS = (
    "PROVISIONAL: evidence-anchored calibrated starting points, NOT final sign-off; "
    "several knobs remain policy/clinical judgment (flagged). Pending human review of the re-run."
)

# Allelic-editing context caveat carried into output so the floor is not misread as clinical proof.
EFFICACY_CAVEAT = (
    "Efficiency values are model predictions in cell-line/library context (BE-DICT). The tier "
    "anchors are ALLELIC editing in mixed peripheral-blood/marrow cells (not purified LT-HSCs), and "
    "allelic % does NOT convert to edited-HSC fraction. The efficacy floor is a relative-comparison "
    "gate for ranking, NOT a claim of clinical sufficiency."
)

# --- Functional consequence of a bystander edit at a promoter offset (negative HBG numbering). ---
# EMERGING evidence; functional class, not bystander count, determines harm. Offsets verified against
# the edited positions the pipeline actually produces for the canonical guides.
# Sources: Mayuranathan 2023 (PMID 37400614, ABE8e bystander vs productive editing); Martyn 2019
# (PMID 30617196, de novo GATA1 at -113); Wienert 2015 (PMID 25971621, TAL1 E-box at -175).
ADVERSE, PRODUCTIVE, NEUTRAL, UNKNOWN = "adverse", "productive", "neutral", "unknown"
FUNCTIONAL_BYSTANDER = {
    -112: (ADVERSE, "EMERGING: disrupts the de novo GATA1 motif created by the -113 edit"),
    -116: (PRODUCTIVE, "EMERGING: co-productive bystander near -113"),
    -114: (PRODUCTIVE, "EMERGING: intended distal-BCL11A-motif edit"),
    -115: (PRODUCTIVE, "EMERGING: part of the intended distal-BCL11A-motif disruption (with -114)"),
    -181: (ADVERSE, "EMERGING: disrupts the created TAL1 E-box near -175 (reasoned offset mapping; confirm)"),
}


def bystander_function(promoter_offset_magnitude: int) -> tuple[str, str]:
    """Functional class + note for a bystander at a promoter offset (positive magnitude in records)."""
    return FUNCTIONAL_BYSTANDER.get(-abs(promoter_offset_magnitude),
                                    (UNKNOWN, "no functional annotation; default mild penalty (flagged)"))


@dataclass(frozen=True)
class RubricWeights:
    """Tunable rubric configuration. Nothing downstream hard-codes these. provisional until sign-off."""

    # ===== EVIDENCE-ANCHORED efficacy tier boundaries (allelic-editing anchors; see EFFICACY_CAVEAT) =====
    clinical_floor: float = 0.30      # POLICY-CHOICE (not a proven minimum): lowest direct TDT clinical anchor (CS-101, n=5)
    preclinical_floor: float = 0.16   # evidence anchor: lower bound of preclinical-only activity
    strong_band_low: float = 0.60     # informational: well-validated clinical operating range 0.60-0.80
    strong_band_high: float = 0.80    # informational
    # ===== FUNCTIONAL bystander penalty magnitudes: UNSETTLED (no literature consensus; policy/clinical judgment) =====
    w_bystander: float = 1.0          # overall bystander-penalty multiplier
    penalty_adverse: float = 1.0      # per unit of adverse (motif-destroying) bystander edit prob
    penalty_unknown: float = 0.3      # mild penalty for unannotated bystanders (flagged)
    penalty_productive: float = 0.0   # productive bystanders are not penalized
    penalty_neutral: float = 0.0
    # ===== off-target penalty: MISMATCH-STRATIFIED by likelihood (not flat count). =====
    # The 0-mm paralog co-target is excluded upstream. Per-site weights fall steeply with mismatch
    # so a 4-mm site contributes far less than a 1-2-mm site (diagnostic: 60forw's 547 is ~91% 4-mm).
    # Magnitudes UNSETTLED (no literature consensus; policy/clinical judgment).
    w_offtarget: float = 1.0
    off_w_mm1: float = 8.0    # 1 mismatch (highest likelihood)
    off_w_mm2: float = 1.0    # 2 mismatches (reference unit)
    off_w_mm3: float = 0.1    # 3 mismatches
    off_w_mm4: float = 0.002  # 4 mismatches (diagnostic: low-likelihood, near-noise)
    # No hard off-target count/VAF gate: the diagnostic showed that would gate on noise.
    offtarget_vaf_threshold: float | None = None   # UNWIRED + flagged (policy/clinical judgment)
    # ===== efficiency reward above floor (false-cleanliness fix): magnitude UNSETTLED (safety-vs-efficiency balance) =====
    w_efficiency_reward: float = 0.5  # mild; applied to max(0, eff - preclinical_floor)
    # ===== consequence term (Step 3 ACTIVE: AlphaGenome regulatory severity): magnitude UNSETTLED =====
    # Penalizes the predicted regulatory severity of NON-productive bystanders (productive bystanders
    # are the intended effect and are not penalized). PROVISIONAL/flagged.
    w_consequence: float = 2.0
    # ===== CH gate conditioning: gate only on a LIKELY, FEATURE-RELEVANT CH-driver off-target =====
    # A CH-driver off-target vetoes (worst tier) only if it is BOTH low-mismatch AND in a coding/UTR
    # exon. A high-mismatch OR intronic/intergenic CH hit is DOWNGRADED to a visible noted flag.
    # Thresholds are UNSETTLED policy/clinical-judgment knobs.
    ch_gate_mm_threshold: int = 3          # max mismatches for a CH hit to count as gating (flagged)
    ch_gate_requires_exon: bool = True     # require exon/UTR feature for a CH hit to gate (flagged)
    ch_gate_trust_proxy_consequence: bool = False  # let the APPROXIMATE A>G proxy consequence gate (flagged)

    # ===== UNSETTLED knobs: present + flagged, NOT yet wired (no literature consensus; policy/clinical judgment) =====
    large_deletion_threshold: float | None = None   # large-deletion size cutoff (policy/clinical judgment)
    safety_vs_efficiency_weight: float | None = None # explicit composite balance (policy/clinical judgment)

    provisional: bool = True          # flips to False only after calibration is reviewed and signed off

    # Tier ordering (lower = better). CH gate is the dominating partition.
    tier_preferred: int = 0
    tier_preclinical: int = 1
    tier_below_activity: int = 2
    tier_ch_flagged: int = 3
    tier_labels = {
        0: "PREFERRED",
        1: "PRECLINICAL_ACTIVITY_ONLY",
        2: "BELOW_EVIDENCED_ACTIVITY",
        3: "CH_DRIVER_FLAGGED__REQUIRES_CONFIRMATION",
    }

    def as_dict(self) -> dict:
        # The tunable knobs, with the unsettled ones flagged; the tier/provisional fields are internal.
        internal = {"provisional", "tier_preferred", "tier_preclinical", "tier_below_activity",
                    "tier_ch_flagged"}
        d = {"status": WEIGHTS_STATUS if self.provisional else "calibrated"}
        d.update({k: v for k, v in asdict(self).items() if k not in internal})
        d["unsettled_policy_knobs"] = [
            "clinical_floor (policy-choice, not proven minimum)", "penalty_adverse", "penalty_unknown",
            "off_w_mm1..4 (off-target mismatch weighting magnitudes)",
            "w_efficiency_reward (safety-vs-efficiency balance)",
            "ch_gate_mm_threshold", "ch_gate_requires_exon",
            "ch_gate_trust_proxy_consequence (approximate A>G proxy off-target consequence)",
            "offtarget_vaf_threshold (unwired)", "large_deletion_threshold", "safety_vs_efficiency_weight",
        ]
        return d


@dataclass
class GuideScore:
    guide_id: str
    gene: str
    editor_class: str
    anchors: list
    on_target_efficiency: float
    efficacy_band: str
    meets_clinical_floor: bool
    bystander_penalty: float
    bystander_breakdown: dict              # rate summed per functional class
    offtarget_count: int                   # raw, includes paralog (reference only)
    offtarget_mm_histogram: dict           # mismatches (>=1) -> site count, paralog excluded
    offtarget_cfd_burden: float            # sum CFD over non-paralog sites (alternative view)
    mit_specificity: str
    ch_gate: bool                          # CH-driver off-target vetoes to worst tier (likely + exonic)
    ch_flag_noted: bool                    # CH-driver off-target present but NOT gating (high-mm or intronic)
    ch_flag_reason: object
    ch_genes: list
    ch_functional_consequence: object      # populated by Step 3 (annotate_consequence.py); weigh LoF > synonymous later
    consequence_severity: float
    tier: int
    tier_label: str
    safety_score: float
    provisional_score_label: str = (
        "provisional within-tier ordering signal; not an absolute safety measure or clearance")
    terms: dict = field(default_factory=dict)
    rank: int = 0

    def as_dict(self) -> dict:
        d = asdict(self)
        # Foreground the primary read (tier + decomposed terms + provisional label); the exact numeric
        # safety_score is a secondary, provisional ordering signal, so it is emitted after those.
        primary = ["guide_id", "gene", "editor_class", "anchors", "tier", "tier_label",
                   "efficacy_band", "provisional_score_label", "terms"]
        ordered = {k: d[k] for k in primary if k in d}
        ordered.update({k: v for k, v in d.items() if k not in ordered})
        return ordered


def _to_int(value, default=0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _efficacy_band(eff: float, w: RubricWeights) -> str:
    if eff >= w.strong_band_low:
        return "STRONG_CLINICAL_RANGE"
    if eff >= w.clinical_floor:
        return "MEETS_CLINICAL_PRECEDENT"
    if eff >= w.preclinical_floor:
        return "PRECLINICAL_ACTIVITY_ONLY"
    return "BELOW_EVIDENCED_ACTIVITY"


def score_guides(records: list[dict], weights: RubricWeights | None = None) -> list[GuideScore]:
    """Rank therapeutic designs (guides that install an in-scope anchor) by safety."""
    w = weights or RubricWeights()
    pen_by_class = {
        ADVERSE: w.penalty_adverse, PRODUCTIVE: w.penalty_productive,
        NEUTRAL: w.penalty_neutral, UNKNOWN: w.penalty_unknown,
    }

    by_guide: dict[str, list[dict]] = {}
    for r in records:
        by_guide.setdefault(r["provenance"]["guide_id"], []).append(r)

    scores: list[GuideScore] = []
    for guide_id, recs in by_guide.items():
        on_target = [r for r in recs if r["on_target_or_bystander"] == "on_target"]
        if not on_target:
            continue
        bystanders = [r for r in recs if r["on_target_or_bystander"] == "bystander"]
        prov = recs[0]["provenance"]
        offtarget = prov.get("offtarget", {})
        ch_hits = prov.get("ch_driver_offtargets", [])

        on_target_efficiency = max((r["predicted_rate"] or 0.0) for r in on_target)

        # Functional bystander penalty: rate weighted by functional class, not bystander count.
        breakdown = {ADVERSE: 0.0, PRODUCTIVE: 0.0, NEUTRAL: 0.0, UNKNOWN: 0.0}
        bystander_penalty = 0.0
        for r in bystanders:
            rate = r["predicted_rate"] or 0.0
            if rate <= 0:
                continue
            cls, _ = bystander_function(r["promoter_offset"])
            breakdown[cls] += rate
            bystander_penalty += pen_by_class[cls] * rate
        bystander_penalty *= w.w_bystander

        offtarget_count = _to_int(offtarget.get("offtarget_count"))
        # Mismatch-stratified off-target penalty (paralog already excluded from mm_histogram).
        mm_weight = {1: w.off_w_mm1, 2: w.off_w_mm2, 3: w.off_w_mm3, 4: w.off_w_mm4}
        histogram = offtarget.get("mm_histogram", {}) or {}
        offt_tier_contrib = {}
        offt_weighted = 0.0
        for k, n in histogram.items():
            try:
                m, count = int(k), int(n)
            except (TypeError, ValueError):
                raise ValueError(f"malformed off-target histogram for {guide_id}: {k!r} -> {n!r} "
                                 f"(expected int mismatch-count -> int site-count)")
            if m < 0 or count < 0:
                raise ValueError(f"malformed off-target histogram for {guide_id}: negative value {m}:{count}")
            wt = mm_weight.get(m, w.off_w_mm4)  # >4 shouldn't occur (search cap=4); treat as 4-mm
            offt_tier_contrib[str(m)] = round(wt * count, 6)
            offt_weighted += wt * count
        cfd_burden = float(offtarget.get("cfd_weighted_burden") or 0.0)

        # CH gate conditioned on likelihood AND functional consequence (AlphaGenome when available,
        # else falling back to the coordinate exon/intron feature).
        gating_hits, noted_hits = [], []
        for h in ch_hits:
            mm = _to_int(h.get("mismatch_count"), 99)
            cfc = h.get("ch_functional_consequence") or {}
            classification = (cfc.get("classification") or "").lower()
            # An APPROXIMATE A>G proxy consequence does not gate on its own unless explicitly trusted;
            # otherwise fall back to the honest coordinate feature (exon/UTR).
            trust_call = bool(classification) and not (cfc.get("is_proxy") and not w.ch_gate_trust_proxy_consequence)
            if trust_call:
                functionally_adverse = classification == "functionally-consequential"
            else:  # no call, or an untrusted proxy call: fall back to the coordinate feature
                functionally_adverse = ((h.get("feature") or "").lower() in ("exon", "utr")
                                        or not w.ch_gate_requires_exon)
            gates = (mm <= w.ch_gate_mm_threshold) and functionally_adverse
            (gating_hits if gates else noted_hits).append(h)
        ch_gate = len(gating_hits) > 0
        ch_flag_noted = (not ch_gate) and len(noted_hits) > 0
        ch_functional_consequence = ch_hits[0].get("ch_functional_consequence") if ch_hits else None
        ch_flag_reason = None
        if ch_flag_noted:
            h = noted_hits[0]
            cls = ((h.get("ch_functional_consequence") or {}).get("classification") or "n/a")
            ch_flag_reason = (f"CH-driver off-target in {h['ch_driver_gene']} noted, NOT gating "
                              f"({h.get('mismatch_count')}-mm, {h.get('feature')}, consequence={cls}); "
                              f"requires empirical confirmation")
        ch_genes = sorted({h["ch_driver_gene"] for h in ch_hits})

        # Consequence penalty: predicted regulatory severity of NON-productive bystanders.
        consequence_severity = 0.0
        for r in bystanders:
            rate = r["predicted_rate"] or 0.0
            cons = r.get("consequence")
            sev = cons.get("severity") if isinstance(cons, dict) else None
            if sev is None or rate <= 0:
                continue
            if bystander_function(r["promoter_offset"])[0] == "productive":
                continue  # intended functional effect, not a safety penalty
            consequence_severity += sev * rate

        band = _efficacy_band(on_target_efficiency, w)
        meets_clinical_floor = on_target_efficiency >= w.clinical_floor

        if ch_gate:
            tier = w.tier_ch_flagged
        elif on_target_efficiency >= w.clinical_floor:
            tier = w.tier_preferred
        elif on_target_efficiency >= w.preclinical_floor:
            tier = w.tier_preclinical
        else:
            tier = w.tier_below_activity

        offt_pen = w.w_offtarget * offt_weighted
        cons_pen = w.w_consequence * consequence_severity
        eff_reward = w.w_efficiency_reward * max(0.0, on_target_efficiency - w.preclinical_floor)
        safety_score = -(bystander_penalty + offt_pen + cons_pen) + eff_reward

        scores.append(GuideScore(
            guide_id=guide_id,
            gene=recs[0]["gene"],
            editor_class=recs[0]["editor_class"],
            anchors=sorted({-abs(r["promoter_offset"]) for r in on_target}),
            on_target_efficiency=round(on_target_efficiency, 6),
            efficacy_band=band,
            meets_clinical_floor=meets_clinical_floor,
            bystander_penalty=round(bystander_penalty, 6),
            bystander_breakdown={k: round(v, 6) for k, v in breakdown.items() if v > 0},
            offtarget_count=offtarget_count,
            offtarget_mm_histogram=histogram,
            offtarget_cfd_burden=round(cfd_burden, 4),
            mit_specificity=offtarget.get("mit_specificity"),
            ch_gate=ch_gate,
            ch_flag_noted=ch_flag_noted,
            ch_flag_reason=ch_flag_reason,
            ch_genes=ch_genes,
            ch_functional_consequence=ch_functional_consequence,  # AlphaGenome regulatory call (Step 3)
            consequence_severity=consequence_severity,
            tier=tier,
            tier_label=RubricWeights.tier_labels[tier],
            # PROVISIONAL within-tier ordering signal, rounded to avoid false precision; the tier and
            # the decomposed `terms` are the primary read, not this exact number (see rank_designs meta).
            safety_score=round(safety_score, 3),
            terms={
                "bystander_penalty": round(bystander_penalty, 6),
                "bystander_breakdown": {k: round(v, 6) for k, v in breakdown.items() if v > 0},
                "offtarget_penalty": round(offt_pen, 6),
                "offtarget_mm_histogram": histogram,          # paralog-excluded raw counts (visible)
                "offtarget_penalty_by_mm": offt_tier_contrib,  # how each mismatch tier contributes
                "offtarget_cfd_burden": round(cfd_burden, 4),  # alternative CFD-weighted view
                "offtarget_paralog_excluded": True,
                "consequence_penalty": round(cons_pen, 6),
                "consequence_severity_nonproductive_bystanders": round(consequence_severity, 6),
                "efficiency_reward": round(eff_reward, 6),
                "efficacy_band": band,
                "ch_gate": ch_gate,
                "ch_flag_noted": ch_flag_noted,
                "ch_flag_reason": ch_flag_reason,
            },
        ))

    # Sort: tier (asc, CH gate dominates) -> safety_score (desc) -> efficiency (desc, tie-break).
    scores.sort(key=lambda s: (s.tier, -s.safety_score, -s.on_target_efficiency))
    for i, s in enumerate(scores, start=1):
        s.rank = i
    return scores
