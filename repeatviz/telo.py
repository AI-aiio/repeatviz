"""Telomere repeats by motif density (v0.3).

ULTRA reports telomere arrays with a period of 6, 12 or 18 bp and often a variant consensus
(TTGGGG-type, AGGGTTCGGGTT, ...), so classifying by the consensus missed most real chromosome
ends that tidk calls. Here the assembly is scanned directly for the exact motif
(TTAGGG on one strand, CCCTAA on the other), like `tidk search`: hits closer than `max_gap` bp
are chained into runs, and a run is a telomere array when it is >= `min_len` bp and at least
`min_density` of it is exact motif. Interstitial runs are kept (repeatviz shows all repeats).
"""
from __future__ import annotations

import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable

MOTIFS = {"+": "TTAGGG", "-": "CCCTAA"}
_RE = {s: re.compile(m) for s, m in MOTIFS.items()}


def scan_seq(name: str, seq: str, max_gap: int = 18, min_len: int = 100,
             min_density: float = 0.5) -> list[tuple[str, int, int, int, str]]:
    """Runs of the telomere motif in one sequence: (contig, start, end, n_motifs, strand)."""
    seq = seq.upper()
    out = []
    for strand, rx in _RE.items():
        s = e = n = 0
        for m in rx.finditer(seq):
            a = m.start()
            if n and a - e <= max_gap:
                e, n = a + 6, n + 1
                continue
            if n and e - s >= min_len and 6 * n / (e - s) >= min_density:
                out.append((name, s, e, n, strand))
            s, e, n = a, a + 6, 1
        if n and e - s >= min_len and 6 * n / (e - s) >= min_density:
            out.append((name, s, e, n, strand))
    out.sort(key=lambda r: r[1])
    return out


def _scan_one(args):
    return scan_seq(*args)


def scan_fasta(records: Iterable[tuple[str, str]], threads: int = 1, **kw) -> list[tuple]:
    jobs = ((n, s, kw.get("max_gap", 18), kw.get("min_len", 100), kw.get("min_density", 0.5))
            for n, s in records)
    runs: list[tuple] = []
    if threads <= 1:
        for j in jobs:
            runs += _scan_one(j)
    else:
        with ProcessPoolExecutor(threads) as ex:
            for r in ex.map(_scan_one, jobs):
                runs += r
    return runs


def write_runs(runs: list[tuple], path: Path) -> None:
    with open(path, "w") as o:
        for c, s, e, n, st in runs:
            o.write(f"{c}\t{s}\t{e}\t{n}\t{st}\n")


def read_runs(path: Path) -> dict[str, list[tuple[int, int, int, str]]]:
    by: dict[str, list] = {}
    for line in open(path):
        c, s, e, n, st = line.rstrip("\n").split("\t")
        by.setdefault(c, []).append((int(s), int(e), int(n), st))
    for v in by.values():
        v.sort()
    return by


def apply(arrays: list, runs: dict[str, list[tuple[int, int, int, str]]], cover: float = 0.5,
          make_call=None) -> tuple[list, int, int]:
    """Label arrays telomere when >= `cover` of their length is telomere-motif runs; add runs that
    no array covers as new arrays (via `make_call(contig, start, end, n, strand)`).
    Returns (arrays, n_relabelled, n_added)."""
    import bisect
    relab = 0
    covered: dict[str, list[tuple[int, int]]] = {}
    for a in arrays:
        iv = runs.get(a.contig)
        if not iv:
            continue
        i = max(0, bisect.bisect_right(iv, (a.start, 10**12)) - 1)
        bp = 0
        while i < len(iv) and iv[i][0] < a.end:
            bp += max(0, min(iv[i][1], a.end) - max(iv[i][0], a.start))
            i += 1
        if bp >= cover * a.length:
            if a.cls != "telomere":
                relab += 1
            a.cls_override = "telomere"
            covered.setdefault(a.contig, []).append((a.start, a.end))
    added = 0
    if make_call is not None:
        for c, iv in runs.items():
            cv = sorted(covered.get(c, []))
            for s, e, n, st in iv:
                j = max(0, bisect.bisect_right(cv, (s, 10**12)) - 1)
                inside = 0
                while j < len(cv) and cv[j][0] < e:
                    inside += max(0, min(cv[j][1], e) - max(cv[j][0], s))
                    j += 1
                if inside < 0.5 * (e - s):
                    arrays.append(make_call(c, s, e, n, st))
                    added += 1
    return arrays, relab, added
