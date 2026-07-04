"""Adapter for the external base-editing bystander outcome predictor.

Upstream tool license: MIT. The user installs the tool in its own isolated environment. This adapter
runs it in a separate process, because its environment differs from the pipeline's, and never copies
its source or weights into this repository. Set its location with the environment variables documented
in config.py.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile

from ..config import BedictConfig, load_bedict_config
from ..pipeline.model import Guide
from ..pipeline.model import EditSpectrum, OutcomeHaplotype

_RUNNER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_bedict_runner.py")

# Default editor identifiers per class when a guide does not name one explicitly.
DEFAULT_EDITOR = {"ABE": "ABE8e", "CBE": "BE4max"}


def _validated_haplotypes(guide: Guide, raw: list[dict]) -> list[OutcomeHaplotype]:
    """Parse + validate predicted haplotypes: right length, probabilities in [0, 1], mass <= 1.

    The external predictor is trusted for the numbers but not for well-formedness; a malformed
    response should fail loudly here rather than produce silently wrong edit rates downstream.
    """
    n = len(guide.protospacer)
    haplotypes: list[OutcomeHaplotype] = []
    mass = 0.0
    for o in raw:
        seq = o.get("output_seq")
        try:
            score = float(o["pred_score"])
        except (KeyError, TypeError, ValueError):
            raise ValueError(f"{guide.guide_id}: non-numeric haplotype probability {o.get('pred_score')!r}")
        if not isinstance(seq, str) or len(seq) != n:
            raise ValueError(f"{guide.guide_id}: haplotype length {len(seq) if isinstance(seq, str) else '?'} "
                             f"!= protospacer length {n}")
        if not 0.0 <= score <= 1.0:
            raise ValueError(f"{guide.guide_id}: haplotype probability {score} outside [0, 1]")
        haplotypes.append(OutcomeHaplotype(output_seq=seq, pred_score=score))
        mass += score
    if mass > 1.0 + 1e-6:
        raise ValueError(f"{guide.guide_id}: haplotype probabilities sum to {mass:.4f} (> 1); malformed distribution")
    return haplotypes


class BedictBystanderPredictor:
    """EditSpectrumPredictor backed by the external bystander model ensemble."""

    name = "BE-DICT-bystander"

    def __init__(self, config: BedictConfig | None = None, runs: int = 5):
        self.config = config or load_bedict_config()
        self.runs = runs

    def predict(self, guide: Guide) -> EditSpectrum:
        """Predict the edit spectrum for a single guide."""
        return self.predict_many([guide])[guide.guide_id]

    def predict_many(self, guides: list[Guide]) -> dict[str, EditSpectrum]:
        """Predict spectra for many guides, with one external call per editor.

        Returns a mapping of guide_id -> EditSpectrum. This groups guides by editor because each
        editor uses its own model set; the external runner then scores all of an editor's sequences
        in one pass.
        """
        by_editor: dict[str, list[Guide]] = {}
        for guide in guides:
            editor = guide.editor_name or DEFAULT_EDITOR[guide.editor_class]
            by_editor.setdefault(editor, []).append(guide)

        spectra: dict[str, EditSpectrum] = {}
        for editor, group in by_editor.items():
            # Within an editor, conversion bases are identical across guides of the same class.
            request = {
                "editor": editor,
                "edit_from": group[0].edit_from,
                "edit_to": group[0].edit_to,
                "runs": self.runs,
                "sequences": [{"seq_id": g.guide_id, "inp_seq": g.protospacer} for g in group],
            }
            response = self._run(request)
            for guide in group:
                spectra[guide.guide_id] = self.spectrum_from_response(guide, editor, response)
        return spectra

    def spectrum_from_response(self, guide: Guide, editor: str, response: dict) -> EditSpectrum:
        """Build a validated EditSpectrum from one predictor response (no subprocess; fixture-testable)."""
        raw = response.get("sequences", {}).get(guide.guide_id, [])
        return EditSpectrum(
            guide=guide,
            model_name=f"{self.name}:{editor}",
            haplotypes=_validated_haplotypes(guide, raw),
            provenance={
                "predictor": self.name,
                "editor": editor,
                "ensemble_runs": response.get("ensemble_runs", self.runs),
                "license": "MIT",
                "execution": "out-of-process; external isolated environment",
            },
        )

    def tool_metadata(self) -> dict:
        """Reproducibility record for this predictor run (paths, run count, git commit if available)."""
        from ..pipeline import environment
        return environment.tool_environment(
            self.name, repo_root=self.config.repo_root,
            python_executable=self.config.python_executable, ensemble_runs=self.runs,
        )

    def _run(self, request: dict) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            req_path = os.path.join(tmp, "request.json")
            resp_path = os.path.join(tmp, "response.json")
            with open(req_path, "w") as handle:
                json.dump(request, handle)
            completed = subprocess.run(
                [self.config.python_executable, _RUNNER, self.config.repo_root, req_path, resp_path],
                capture_output=True,
                text=True,
            )
            if completed.returncode != 0 or not os.path.isfile(resp_path):
                raise RuntimeError(
                    "External bystander predictor failed.\n"
                    f"stderr (tail):\n{completed.stderr[-2000:]}"
                )
            with open(resp_path) as handle:
                return json.load(handle)
