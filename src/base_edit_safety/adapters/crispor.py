"""Adapter for the external guide-design / off-target tool.

Upstream tool license: see the tool's LICENSE (CRISPOR core is under a tool-specific license;
bundled aligners are GPL). The user installs the tool with its own indexed genome. This adapter
runs it in a separate process and never copies its source, binaries, or genome into this
repository. Set its location with the environment variables documented in config.py.

The tool returns guides as editor-agnostic CandidateGuides. For each guide, the adapter finds its
protospacer+PAM in the reference window to set its position and strand, then attaches the tool's
off-target summary (specificity scores, off-target count, locus).
"""
from __future__ import annotations

import csv
import os
import subprocess
import tempfile

from ..config import CrisporConfig, load_crispor_config
from ..io.reference import fetch_window
from ..pipeline.model import reverse_complement
from ..pipeline.model import CandidateGuide

PROTO_LEN = 20


class AmbiguousTargetMapping(ValueError):
    """Raised when a CRISPOR target maps to more than one site in the reference window."""


def _find_all(haystack: str, needle: str) -> list[int]:
    """All start indices where needle occurs in haystack (non-overlapping)."""
    hits, idx = [], haystack.find(needle)
    while idx != -1:
        hits.append(idx)
        idx = haystack.find(needle, idx + 1)
    return hits


def _offtarget_likelihood_summary(sites: list[dict]) -> dict:
    """Likelihood facts about a guide's off-targets, for the rubric's mismatch-weighted penalty.

    The 0-mismatch site(s) are the intended HBG1/HBG2 paralog co-target, not off-targets. This counts
    them separately (paralog_count) and leaves them out of the off-target burden inputs. Returns the
    raw per-mismatch histogram and the CFD-weighted burden (both are facts; the rubric applies the
    policy weights). Mismatch and CFD come straight from the tool's per-site output.
    """
    histogram: dict[str, int] = {}
    cfd_burden = 0.0
    paralog_count = 0
    for s in sites:
        try:
            mm = int(s.get("mismatchCount"))
        except (TypeError, ValueError):
            continue
        if mm == 0:
            paralog_count += 1
            continue  # intended paralog co-target, not an off-target
        histogram[str(mm)] = histogram.get(str(mm), 0) + 1
        try:
            cfd_burden += float(s.get("cfdOfftargetScore"))
        except (TypeError, ValueError):
            pass
    return {
        "mm_histogram": histogram,                 # mismatches (>=1) -> site count, paralog excluded
        "cfd_weighted_burden": round(cfd_burden, 4),
        "paralog_count": paralog_count,
        "offtarget_count_excl_paralog": sum(histogram.values()),
    }


class CrisporGuideDesigner:
    """GuideDesigner backed by the external off-target-aware guide-design tool."""

    name = "CRISPOR"

    def __init__(self, config: CrisporConfig | None = None, pam: str = "NGG",
                 tss_for_gene=None):
        self.config = config or load_crispor_config()
        self.pam = pam
        # Map a gene symbol to its TSS; defaults to the Step 0 anchor values.
        if tss_for_gene is None:
            from ..io import step0_anchor as anchor
            tss_for_gene = {"HBG1": anchor.HBG1_TSS, "HBG2": anchor.HBG2_TSS}
        self.tss_for_gene = tss_for_gene

    def design(self, chrom: str, start_1based: int, end_1based: int, gene: str) -> list[CandidateGuide]:
        window = fetch_window(chrom, start_1based, end_1based)
        rows, offtarget_rows = self._run_crispor(
            window.plus_sequence, f"{gene}_{chrom}_{start_1based}_{end_1based}"
        )
        sites_by_guide: dict[str, list[dict]] = {}
        for site in offtarget_rows:
            sites_by_guide.setdefault(site.get("guideId", ""), []).append(site)
        guides: list[CandidateGuide] = []
        self.last_ambiguous: list[str] = []   # flagged, not silently mapped to the first hit
        for row in rows:
            try:
                guide = self._row_to_candidate(row, window, gene, sites_by_guide)
            except AmbiguousTargetMapping as exc:
                self.last_ambiguous.append(str(exc))
                continue
            if guide is not None:
                guides.append(guide)
        return guides

    def tool_metadata(self) -> dict:
        """Reproducibility record for this tool run (paths, genome index, git commit if available)."""
        from ..pipeline import environment
        return environment.tool_environment(
            self.name, repo_root=self.config.repo_root,
            python_executable=self.config.python_executable, genome=self.config.genome, pam=self.pam,
        )

    def _run_crispor(self, window_seq: str, seq_id: str) -> list[dict]:
        with tempfile.TemporaryDirectory() as tmp:
            in_fa = os.path.join(tmp, "in.fa")
            out_tsv = os.path.join(tmp, "out.tsv")
            off_tsv = os.path.join(tmp, "off.tsv")
            with open(in_fa, "w") as handle:
                handle.write(f">{seq_id}\n{window_seq}\n")
            cmd = [
                self.config.python_executable, "crispor.py", self.config.genome,
                in_fa, out_tsv, "-p", self.pam, "--noEffScores",
                "--genomeDir", self.config.genome_dir, "-o", off_tsv,
            ]
            completed = subprocess.run(
                cmd, cwd=self.config.repo_root, capture_output=True, text=True
            )
            if completed.returncode != 0 or not os.path.isfile(out_tsv):
                raise RuntimeError(
                    "External guide-design tool failed.\n"
                    f"stderr (tail):\n{completed.stderr[-2000:]}"
                )
            with open(out_tsv) as handle:
                guide_rows = [dict(r) for r in csv.DictReader(handle, delimiter="\t")]
            offtarget_rows = []
            if os.path.isfile(off_tsv):
                with open(off_tsv) as handle:
                    offtarget_rows = [dict(r) for r in csv.DictReader(handle, delimiter="\t")]
            return guide_rows, offtarget_rows

    def _map_target(self, window, target_seq: str):
        """Locate target_seq in the window. Returns (strand, 5'-genomic-coord), or (None, None) if it
        is absent. Raises AmbiguousTargetMapping if it occurs at more than one position or strand, so
        an ambiguous mapping is never silently resolved to the first hit."""
        plus = window.plus_sequence
        plus_hits = _find_all(plus, target_seq)
        minus_hits = _find_all(reverse_complement(plus), target_seq)
        total = len(plus_hits) + len(minus_hits)
        if total == 0:
            return None, None
        if total > 1:
            raise AmbiguousTargetMapping(
                f"target {target_seq} maps to {total} sites in the reference window "
                f"(plus={plus_hits}, minus={minus_hits}); refusing to pick one silently."
            )
        if plus_hits:
            return "+", window.start_1based + plus_hits[0]   # 5' end at low coord on plus strand
        return "-", window.end_1based - minus_hits[0]         # 5' end at high coord on minus strand

    def _row_to_candidate(self, row: dict, window, gene: str,
                          sites_by_guide: dict | None = None) -> CandidateGuide | None:
        target_seq = row.get("targetSeq")
        if not target_seq or len(target_seq) < PROTO_LEN + 1:
            return None
        protospacer = target_seq[:PROTO_LEN]
        pam = target_seq[PROTO_LEN:]

        strand, proto_5p_genomic = self._map_target(window, target_seq)
        if strand is None:
            return None  # not present in this window (e.g. the guide spans the edge)

        tool_guide_id = row.get("guideId", "")
        sites = (sites_by_guide or {}).get(tool_guide_id, [])
        offtarget = {
            "source": self.name,
            "mit_specificity": row.get("mitSpecScore"),
            "cfd_specificity": row.get("cfdSpecScore"),
            "offtarget_count": row.get("offtargetCount"),  # raw count, includes the paralog site
            "genome_locus": row.get("targetGenomeGeneLocus"),
            "tool_guide_id": tool_guide_id,
            **_offtarget_likelihood_summary(sites),
        }
        return CandidateGuide(
            guide_id=f"{gene}_{tool_guide_id or 'NA'}",
            protospacer=protospacer,
            pam=pam,
            strand=strand,
            genome_build=window.genome_build,
            chrom=window.chrom,
            gene=gene,
            gene_strand="-",
            tss_1based=self.tss_for_gene[gene],
            proto_5p_genomic=proto_5p_genomic,
            offtarget=offtarget,
            offtarget_sites=sites,
        )
