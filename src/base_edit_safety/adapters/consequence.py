"""Functional-consequence annotation (Step 3).

Two paths live here. `annotate_all` is the cheap default the spine and enumeration steps use: it
marks every record's consequence as 'deferred' so records stay well-formed without calling any
external model. `AlphaGenomeConsequenceAnnotator` is the real thing, run later by annotate_consequence.py
to fill the consequence in with actual regulatory-perturbation scores.
"""
from __future__ import annotations

import json
import os

from ..pipeline.model import EditRecord

DEFERRED_VALUE = "deferred"
_DEFERRED_SOURCE = "deferred-consequence-passthrough"


def annotate_all(records: list[EditRecord]) -> list[EditRecord]:
    """Mark every record's consequence as deferred (the Step 3 seam, mutates in place)."""
    for record in records:
        record.consequence = DEFERRED_VALUE
        record.provenance = {
            **record.provenance,
            "consequence_source": _DEFERRED_SOURCE,
            "consequence_note": ("regulatory/coding consequence not evaluated in this build; "
                                 "slot reserved for external annotator"),
        }
    return records


# --- Real regulatory-consequence annotator (AlphaGenome) ------------------------------------------
# Upstream license: NON-COMMERCIAL (AlphaGenome Terms of Service). The user installs AlphaGenome as an
# external dependency, and this adapter calls it through its client; this repository never copies its
# source or weights. Its terms say AlphaGenome outputs MUST NOT be used to train any model. The API key
# is read only from ALPHAGENOME_API_KEY and is never logged, printed, or written to disk.
#
# Severity for a variant is the mean |quantile_score| across the model's recommended output tracks. The
# per-track max sits near 1.0 and cannot tell tracks apart; the mean can (functional HBG promoter edits
# score ~0.65-0.75 vs ~0.39 for a distant control). Severity is a magnitude of predicted regulatory
# perturbation, not a direction, and NOT a clonal/lifetime-risk value (G3).

DEFAULT_FUNCTIONAL_THRESHOLD = 0.5
_SEQUENCE_LENGTH = "SEQUENCE_LENGTH_16KB"
_CACHE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data", "inputs", "alphagenome_cache.json",
)


class AlphaGenomeConsequenceAnnotator:
    """Scores one variant's regulatory-consequence severity via the external AlphaGenome client."""

    name = "AlphaGenome"

    def __init__(self, functional_threshold: float = DEFAULT_FUNCTIONAL_THRESHOLD, use_cache: bool = True):
        self.functional_threshold = functional_threshold  # PROVISIONAL/flagged
        self.use_cache = use_cache
        self._model = None
        self._cache = self._load_cache() if use_cache else {}

    # --- model + cache plumbing -------------------------------------------------
    def _ensure_model(self):
        if self._model is None:
            from alphagenome.models import dna_client  # lazy: external dep
            key = os.environ.get("ALPHAGENOME_API_KEY")
            if not key:
                raise RuntimeError("ALPHAGENOME_API_KEY not set; cannot reach the consequence model.")
            self._model = dna_client.create(key)  # key stays in memory only
        return self._model

    def _load_cache(self) -> dict:
        try:
            with open(_CACHE_PATH) as handle:
                return json.load(handle)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def _save_cache(self):
        if not self.use_cache:
            return
        os.makedirs(os.path.dirname(_CACHE_PATH), exist_ok=True)
        with open(_CACHE_PATH, "w") as handle:
            json.dump(self._cache, handle)  # severity/classification only; never the key

    # --- scoring ----------------------------------------------------------------
    def score(self, chrom: str, pos_1based: int, ref: str, alt: str) -> dict:
        """Return regulatory-consequence severity + classification for one variant (cached)."""
        vkey = f"{chrom}:{pos_1based}:{ref}>{alt}"
        if self.use_cache and vkey in self._cache:
            return self._cache[vkey]

        import numpy as np
        from alphagenome.models import dna_client, variant_scorers
        from alphagenome.data import genome

        model = self._ensure_model()
        variant = genome.Variant(chromosome=chrom, position=pos_1based,
                                 reference_bases=ref, alternate_bases=alt)
        interval = variant.reference_interval.resize(getattr(dna_client, _SEQUENCE_LENGTH))
        scores = model.score_variant(interval=interval, variant=variant, variant_scorers=None)
        df = variant_scorers.tidy_scores(scores)
        q = np.abs(df["quantile_score"].astype(float).to_numpy())
        q = q[~np.isnan(q)]
        severity = float(q.mean()) if q.size else None
        top = df.iloc[int(np.nanargmax(np.abs(df["quantile_score"].astype(float))))]
        result = {
            "source": self.name,
            "license": "non-commercial; outputs not used to train models",
            "model_predicted": True,
            "measured": False,
            "severity": round(severity, 4) if severity is not None else None,
            "severity_metric": "mean|quantile_score| over recommended tracks",
            "classification": self._classify(severity),
            "top_output_type": str(top.get("output_type")),
            "top_biosample": str(top.get("biosample_name")),
            "n_tracks": int(len(df)),
        }
        if self.use_cache:
            self._cache[vkey] = result
            self._save_cache()
        return result

    def _classify(self, severity) -> str:
        if severity is None:
            return "unknown"
        return "functionally-consequential" if severity >= self.functional_threshold else "low-consequence"

    def annotate(self, record: EditRecord) -> EditRecord:
        """Populate a record's consequence from the genomic-plus allele (regulatory, noncoding)."""
        result = self.score(record.chrom, record.pos_1based,
                            record.genomic_plus_ref, record.genomic_plus_alt)
        record.consequence = result
        record.provenance = {
            **record.provenance,
            "consequence_source": self.name,
            "consequence_note": "regulatory consequence MODEL-PREDICTED (AlphaGenome); requires wet-lab confirmation",
        }
        return record
