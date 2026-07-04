"""Read reference sequence for GRCh38 promoter windows.

Fetches the plus-strand reference sequence for a genomic window and caches it. The pipeline uses
it to read flanking context when building guides and to set reference alleles. The cache lives
under data/inputs/ (gitignored), so repeated runs do not re-download.
"""
from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from dataclasses import dataclass

_UCSC_SEQUENCE_API = "https://api.genome.ucsc.edu/getData/sequence"
_CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data",
    "inputs",
    "reference_cache",
)


@dataclass(frozen=True)
class ReferenceWindow:
    """Plus-strand reference sequence for chrom:[start_1based, end_1based] inclusive."""

    genome_build: str
    chrom: str
    start_1based: int
    end_1based: int
    plus_sequence: str


def _ucsc_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def fetch_window(
    chrom: str,
    start_1based: int,
    end_1based: int,
    genome_build: str = "hg38",
    use_cache: bool = True,
    timeout: int = 30,
) -> ReferenceWindow:
    """Fetch plus-strand reference for an inclusive 1-based window, with on-disk caching."""
    key = f"{genome_build}_{_ucsc_chrom(chrom)}_{start_1based}_{end_1based}.json"
    cache_path = os.path.join(_CACHE_DIR, key)
    if use_cache and os.path.isfile(cache_path):
        with open(cache_path) as handle:
            payload = json.load(handle)
        return ReferenceWindow(**payload)

    start0 = start_1based - 1  # UCSC API uses 0-based, half-open coordinates
    url = (
        f"{_UCSC_SEQUENCE_API}?genome={genome_build};chrom={_ucsc_chrom(chrom)}"
        f";start={start0};end={end_1based}"
    )
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = json.load(response)
    seq = data["dna"].upper()
    window = ReferenceWindow(
        genome_build=genome_build,
        chrom=_ucsc_chrom(chrom),
        start_1based=start_1based,
        end_1based=end_1based,
        plus_sequence=seq,
    )
    if use_cache:
        os.makedirs(_CACHE_DIR, exist_ok=True)
        with open(cache_path, "w") as handle:
            json.dump(window.__dict__, handle)
    return window


def sequence_checksum(seq: str) -> str:
    """Short sha256 of a sequence, for reproducibility provenance."""
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def window_provenance(chrom: str, start_1based: int, end_1based: int,
                      genome_build: str = "hg38") -> dict:
    """Reproducibility facts about a cached reference window (no network).

    Returns the cache key, whether it is cached, and (if cached) the sequence checksum and length.
    """
    key = f"{genome_build}_{_ucsc_chrom(chrom)}_{start_1based}_{end_1based}.json"
    path = os.path.join(_CACHE_DIR, key)
    if not os.path.isfile(path):
        return {"cache_key": key, "cached": False}
    with open(path) as handle:
        seq = json.load(handle).get("plus_sequence", "")
    return {"cache_key": key, "cached": True,
            "checksum_sha256_16": sequence_checksum(seq), "length": len(seq)}
