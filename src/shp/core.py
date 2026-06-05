"""Core SHP implementation.

SHP compares two binary views of the same local symbol window:

- chroma: which n-grams are present;
- rhythm: which adjacent n-gram transitions are present.

The Jaccard distance between those views is the local cross-harm trace. Changes
in cross-harm above a calibrated fair-IID threshold define structural events.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
import math
from pathlib import Path
from typing import Iterable, Iterator


DEFAULT_THETA0 = 0.0999


@dataclass(frozen=True)
class SHPResult:
    """Summary of one SHP sequence scan."""

    length: int
    n_windows: int
    n_transitions: int
    mean_h: float
    mean_d: float
    fixed_wit: float
    tail_energy: float
    skew_d: float
    kurt_d: float


def _hash_bucket(text: str, dim: int) -> int:
    digest = hashlib.md5(text.encode("utf-8")).hexdigest()
    return int(digest, 16) % dim


def _windows(seq: str, window: int, stride: int) -> Iterator[str]:
    if len(seq) < window:
        return
    for start in range(0, len(seq) - window + 1, stride):
        yield seq[start : start + window]


def _jaccard_distance(a: set[int], b: set[int]) -> float:
    union = a | b
    if not union:
        return 0.0
    return 1.0 - len(a & b) / len(union)


def _moments(values: list[float]) -> tuple[float, float]:
    if len(values) < 2:
        return 0.0, 0.0
    mean = sum(values) / len(values)
    var = sum((x - mean) ** 2 for x in values) / len(values)
    if var <= 0:
        return 0.0, -3.0
    sd = math.sqrt(var)
    skew = sum(((x - mean) / sd) ** 3 for x in values) / len(values)
    kurt = sum(((x - mean) / sd) ** 4 for x in values) / len(values) - 3.0
    return skew, kurt


def normalize_dna(seq: str) -> str:
    """Keep only A/C/G/T symbols and uppercase them."""

    return "".join(ch for ch in seq.upper() if ch in "ACGT")


def compute_shp(
    seq: str,
    *,
    ngram: int = 3,
    dim: int = 64,
    window: int = 128,
    stride: int | None = None,
    theta0: float = DEFAULT_THETA0,
    alphabet: str = "ACGT",
) -> SHPResult:
    """Compute SHP summary metrics for one sequence.

    Parameters are deliberately explicit so the calibration setting is visible.
    The default DNA setting is k=4, n=3, D=64, W=128, theta0=0.0999.
    """

    if stride is None:
        stride = max(1, window // 5)
    allowed = set(alphabet.upper())
    clean = "".join(ch for ch in seq.upper() if ch in allowed)
    hs: list[float] = []

    for win in _windows(clean, window, stride):
        kmers = [win[i : i + ngram] for i in range(0, len(win) - ngram + 1)]
        chroma = {_hash_bucket("C:" + k, dim) for k in kmers}
        rhythm = {_hash_bucket("R:" + a + ">" + b, dim) for a, b in zip(kmers[:-1], kmers[1:])}
        hs.append(_jaccard_distance(chroma, rhythm))

    ds = [abs(hs[i] - hs[i - 1]) for i in range(1, len(hs))]
    fixed_wit = sum(1 for d in ds if d > theta0) / len(ds) if ds else 0.0
    tail_energy = sum(max(0.0, d - theta0) for d in ds) / len(ds) if ds else 0.0
    skew, kurt = _moments(ds)

    return SHPResult(
        length=len(clean),
        n_windows=len(hs),
        n_transitions=len(ds),
        mean_h=sum(hs) / len(hs) if hs else 0.0,
        mean_d=sum(ds) / len(ds) if ds else 0.0,
        fixed_wit=fixed_wit,
        tail_energy=tail_energy,
        skew_d=skew,
        kurt_d=kurt,
    )


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (record_id, sequence) pairs from plain or gzipped FASTA."""

    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8", errors="ignore") as f:
        name: str | None = None
        parts: list[str] = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            yield name, "".join(parts)


def scan_fasta(path: str | Path, **kwargs: object) -> Iterator[tuple[str, SHPResult]]:
    """Compute SHP metrics for every FASTA record."""

    for name, seq in read_fasta(path):
        yield name, compute_shp(seq, **kwargs)
