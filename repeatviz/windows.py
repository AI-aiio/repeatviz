"""Per-window tracks: multi-copy k-mer coverage and tandem-repeat coverage by class."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from .tandem import CLASSES


def read_bed_runs(path: Path) -> dict[str, list[tuple[int, int]]]:
    """meryl-lookup -bed-runs output (0-based half-open; may carry a label column)."""
    runs: dict[str, list[tuple[int, int]]] = {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3:
                continue
            runs.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    for v in runs.values():
        v.sort()
    return runs


def merge_intervals(iv: list[tuple[int, int]], gap: int = 0) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def binned_coverage(iv: list[tuple[int, int]], length: int, win: int) -> np.ndarray:
    """bp covered per window by (merged) intervals."""
    n = (length + win - 1) // win
    cov = np.zeros(n, dtype=np.int64)
    for s, e in merge_intervals(iv):
        s, e = max(0, s), min(length, e)
        while s < e:
            w = s // win
            we = min(e, (w + 1) * win)
            cov[w] += we - s
            s = we
    return cov


def window_sizes(length: int, win: int) -> np.ndarray:
    n = (length + win - 1) // win
    sz = np.full(n, win, dtype=np.int64)
    if n:
        sz[-1] = length - (n - 1) * win
    return sz


def high_copy_regions(runs: dict[str, list[tuple[int, int]]], lengths: dict[str, int],
                      merge_gap: int = 1000, min_len: int = 10000, min_frac: float = 0.8,
                      pad: int = 1000) -> dict[str, list[tuple[int, int]]]:
    """Long, densely high-copy stretches (satellite arrays, collapsed SDs).

    Runs closer than `merge_gap` are joined; a region is kept if it is at least
    `min_len` long and >= `min_frac` covered by high-copy k-mers (this excludes
    dispersed SINE/LINE copies). Kept regions are padded by `pad`.
    """
    out: dict[str, list[tuple[int, int]]] = {}
    for c, iv in runs.items():
        merged: list[list[int]] = []
        for s, e in iv:
            if merged and s <= merged[-1][1] + merge_gap:
                merged[-1][1] = max(merged[-1][1], e)
                merged[-1][2] += e - s
            else:
                merged.append([s, e, e - s])
        L = lengths.get(c, 0)
        keep = [(max(0, s - pad), min(L, e + pad)) for s, e, cov in merged
                if e - s >= min_len and cov / (e - s) >= min_frac]
        if keep:
            out[c] = merge_intervals(keep)
    return out


def build_windows(lengths: dict[str, int], win: int,
                  ge2: dict[str, list[tuple[int, int]]],
                  ge10: dict[str, list[tuple[int, int]]],
                  arrays: pd.DataFrame) -> pd.DataFrame:
    """One row per window: fraction low-copy (k-mer count >=2), high-copy (>=10),
    and tandem coverage per class."""
    by_cls: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for r in arrays.itertuples(index=False):
        by_cls.setdefault((r.contig, r[6]), []).append((r.start, r.end))
    rows = []
    for c, L in lengths.items():
        sz = window_sizes(L, win)
        d = {
            "contig": c,
            "start": np.arange(len(sz)) * win,
        }
        d["end"] = d["start"] + sz
        d["frac_ge2"] = binned_coverage(ge2.get(c, []), L, win) / sz
        d["frac_ge10"] = binned_coverage(ge10.get(c, []), L, win) / sz
        tot = np.zeros(len(sz), dtype=np.int64)
        for k in CLASSES:
            cov = binned_coverage(by_cls.get((c, k), []), L, win)
            d[f"tandem_{k}"] = cov / sz
            tot = tot + cov
        d["tandem_all"] = np.minimum(tot / sz, 1.0)
        rows.append(pd.DataFrame(d))
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return df.round(4)
