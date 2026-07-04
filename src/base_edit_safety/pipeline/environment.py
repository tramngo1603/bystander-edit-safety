"""Best-effort environment + reproducibility metadata carried into pipeline outputs.

Records what a run depended on (tool paths and git commits when available, genome build/index,
editors, ensemble run count, reference-window checksums) so a written artifact can be traced back to
the exact inputs that produced it. Everything here is best-effort and never raises.
"""
from __future__ import annotations

import datetime
import os
import subprocess

from ..io import step0_anchor as anchor

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def git_commit(path: str | None) -> str | None:
    """Short git commit of a checkout, or None. Best-effort; never raises."""
    if not path:
        return None
    try:
        done = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=5)
    except Exception:
        return None
    return (done.stdout.strip() or None) if done.returncode == 0 else None


def tool_environment(name: str, repo_root: str | None = None,
                     python_executable: str | None = None, genome: str | None = None,
                     **extra) -> dict:
    """Reproducibility record for one external tool: exact paths, genome index, and a best-effort git
    commit. Version lookup never breaks a run; when it fails, that is recorded visibly (not dropped)."""
    env: dict = {"name": name}
    if repo_root is not None:
        env["repo_root"] = repo_root
        commit = git_commit(repo_root)
        env["commit"] = commit
        env["commit_lookup"] = "ok" if commit else "unavailable"
    if python_executable is not None:
        env["python_executable"] = python_executable
    if genome is not None:
        env["genome_index"] = genome
    env.update(extra)
    return env


def reproducibility(tools: list[dict], editors: dict, ensemble_runs: int | None,
                    reference_windows: dict) -> dict:
    """Assemble the reproducibility block written into an output's meta."""
    return {
        "genome_build": anchor.GENOME_BUILD,
        "tools": tools,
        "editors": editors,
        "ensemble_runs": ensemble_runs,
        "reference_windows": reference_windows,
    }


def repo_state() -> dict:
    """This repository's git commit and dirty state (best-effort; never raises)."""
    commit = git_commit(_REPO_ROOT)
    dirty = None
    try:
        done = subprocess.run(["git", "-C", _REPO_ROOT, "status", "--porcelain"],
                              capture_output=True, text=True, timeout=5)
        dirty = bool(done.stdout.strip()) if done.returncode == 0 else None
    except Exception:
        dirty = None
    return {"commit": commit, "dirty": dirty, "commit_lookup": "ok" if commit else "unavailable"}


def run_stamp() -> str:
    """UTC timestamp for a run artifact (ISO 8601)."""
    return datetime.datetime.now(datetime.timezone.utc).isoformat()
