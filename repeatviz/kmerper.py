"""Light long pass: tandem-repeat period from k-mer return distances (no ULTRA).

For every position of a high-copy region the distance to the next (or, for the
last copy, the previous) occurrence of the same k-mer is computed. Inside a
tandem array with unit length P most of these distances are P (or 2P, 3P when a
copy is mutated at that k-mer). The region is cut into tiles; a tile whose
distance mode P explains >= `min_frac` of its positions (as P, 2P, 3P, 4P) is
tandem with period P. Adjacent tandem tiles of similar period are merged into
arrays, whose ends are then trimmed to the first/last supported position.

Cost: one sort per region (O(n log n)), independent of the unit length, so there
is no maximum period (a 40 kb unit is found as easily as a 171 bp one).
What is lost compared with ULTRA: a real consensus (the unit reported is one raw
copy from the middle of the array) and bp-exact boundaries in mixed arrays.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

_LUT = np.full(256, 255, dtype=np.uint8)
for _b, _v in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    _LUT[_b] = _v


@dataclass
class KArray:
    contig: str
    start: int          # 0-based
    end: int            # exclusive
    period: int
    tandem_frac: float  # share of valid positions in the array explained by period (x1..x4)
    period2: int = 0    # other, non-harmonic return-distance peaks (HOR / mixed arrays)
    period3: int = 0
    unit: str = ""

    @property
    def length(self) -> int:
        return self.end - self.start


def kmer_codes(seq: str, k: int) -> tuple[np.ndarray, np.ndarray]:
    """2-bit packed forward k-mers (uint64) and a mask of k-mers without N."""
    a = _LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]
    m = len(a) - k + 1
    if m <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, bool)
    bad = a == 255
    a64 = np.where(bad, 0, a).astype(np.uint64)
    codes = np.zeros(m, np.uint64)
    two = np.uint64(2)
    for j in range(k):
        np.left_shift(codes, two, out=codes)
        np.bitwise_or(codes, a64[j:j + m], out=codes)
    cs = np.concatenate(([0], np.cumsum(bad, dtype=np.int64)))
    valid = (cs[k:k + m] - cs[:m]) == 0
    return codes, valid


def return_distance(codes: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Distance to the next occurrence of the same k-mer; for the last occurrence,
    the distance to the previous one; 0 = k-mer occurs once (or invalid)."""
    n = len(codes)
    d = np.zeros(n, np.int64)
    idx = np.flatnonzero(valid)
    if len(idx) < 2:
        return d
    c = codes[idx]
    order = np.argsort(c, kind="stable")       # stable: equal codes stay in position order
    sp = idx[order]
    same = c[order][1:] == c[order][:-1]
    gap = sp[1:] - sp[:-1]
    d[sp[1:][same]] = gap[same]                 # previous distance first ...
    d[sp[:-1][same]] = gap[same]                # ... overwritten by next distance where it exists
    return d


def _logbins(x: np.ndarray, tol: float) -> np.ndarray:
    return np.round(np.log(x) / math.log1p(tol)).astype(np.int64)


def _mode(dd: np.ndarray, tol: float) -> int:
    b = _logbins(dd, tol)
    b0 = b.min()
    cnt = np.bincount(b - b0)
    # smooth over 3 bins so an indel-blurred peak beats a sharp minor one
    sm = cnt.copy()
    sm[1:] += cnt[:-1]
    sm[:-1] += cnt[1:]
    mb = int(np.argmax(sm)) + b0
    sel = np.abs(b - mb) <= 1
    return int(np.median(dd[sel]))


def _harmonic_mask(dd: np.ndarray, p: int, tol: float, mmax: int = 4) -> np.ndarray:
    m = np.clip(np.rint(dd / p), 1, mmax)
    return (np.abs(dd - m * p) <= np.maximum(2, tol * m * p)) & (dd <= (mmax + 0.5) * p)


def _other_peaks(dd: np.ndarray, p: int, tol: float, min_share: float, n: int = 2) -> list[int]:
    """Up to `n` further peaks of the distance histogram that are not harmonics of
    `p` (nor of each other): monomer of a HOR, or a second unit in a mixed array."""
    rest = dd[~_harmonic_mask(dd, p, tol)]
    out: list[int] = []
    need = max(20, int(min_share * len(dd)))
    for _ in range(n):
        if len(rest) < need:
            break
        q = _mode(rest, tol)
        sel = _harmonic_mask(rest, q, tol)
        if sel.sum() < need:
            break
        out.append(q)
        rest = rest[~sel]
    return out


def analyse_region(contig: str, offset: int, seq: str, k: int = 16, tile: int = 5000,
                   tol: float = 0.03, min_frac: float = 0.3, min_period: int = 101,
                   min_copies: float = 2.0, min_len: int = 1000) -> list[KArray]:
    """Tandem arrays with period >= `min_period` in one region (coordinates on the contig)."""
    codes, valid = kmer_codes(seq, k)
    n = len(codes)
    if n < 2 * k:
        return []
    d = return_distance(codes, valid)
    # per tile: dominant distance and its support
    tiles = []   # (start, end, period) of tandem tiles
    for s in range(0, n, tile):
        e = min(n, s + tile)
        nv = int(valid[s:e].sum())
        if nv < 0.5 * (e - s):
            continue
        dd = d[s:e]
        dd = dd[dd > 0]
        if len(dd) < min_frac * nv or len(dd) < 20:
            continue
        p = _mode(dd, tol)
        if p < min_period:
            continue
        if _harmonic_mask(dd, p, tol).sum() >= min_frac * nv:
            tiles.append([s, e, p])
    # merge tiles of similar period (one non-tandem tile in between is tolerated)
    arrays: list[list[int]] = []
    for s, e, p in tiles:
        if arrays:
            a = arrays[-1]
            if s <= a[1] + tile and abs(p - a[2]) <= max(2, tol * max(p, a[2])):
                a[1] = e
                a[3].append(p)
                continue
        arrays.append([s, e, p, [p]])
    out: list[KArray] = []
    for s, e, _, ps in arrays:
        p = int(np.median(ps))
        # a long unit spans several tiles: extend by one unit on each side and trim to support
        s2, e2 = max(0, s - p), min(n, e + p)
        dd = d[s2:e2]
        sup = (dd > 0) & _harmonic_mask(np.maximum(dd, 1), p, tol)
        pos = np.flatnonzero(sup)
        if len(pos) == 0:
            continue
        a0, a1 = s2 + int(pos[0]), s2 + int(pos[-1]) + k
        if out and a0 < out[-1].end - offset:          # overlaps the previous array: clip
            a0 = out[-1].end - offset
        length = a1 - a0
        if length < min_len or length / p < min_copies:
            continue
        sup = sup[a0 - s2:a1 - k + 1 - s2]
        nv = max(1, int(valid[a0:a1 - k + 1].sum()))
        frac = float(sup.sum()) / nv
        if frac < min_frac:
            continue
        dd = d[a0:a1 - k + 1]
        dd = dd[(dd > 0) & (dd <= length)]           # distances to other arrays are not peaks
        others = _other_peaks(dd, p, tol, 0.05)
        mid = a0 + (length - p) // 2 if length > p else a0
        out.append(KArray(contig, offset + a0, offset + a1, p, round(min(frac, 1.0), 3),
                          others[0] if others else 0, others[1] if len(others) > 1 else 0,
                          seq[mid:mid + p].upper()))
    return out


def _worker(job):
    contig, offset, seq, kw = job
    return analyse_region(contig, offset, seq, **kw)


def run_regions(fasta_iter, regions: dict[str, list[tuple[int, int]]], threads: int = 1,
                **kw) -> list[KArray]:
    """Analyse every region. `fasta_iter` yields (name, seq). Big regions are sent first."""
    jobs = []
    for name, seq in fasta_iter:
        for s, e in regions.get(name, ()):
            jobs.append((name, s, seq[s:e], kw))
    jobs.sort(key=lambda j: -len(j[2]))
    if threads <= 1 or len(jobs) < 2:
        res = [_worker(j) for j in jobs]
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(threads) as pool:
            res = list(pool.imap_unordered(_worker, jobs, chunksize=1))
    out = [a for r in res for a in r]
    out.sort(key=lambda a: (a.contig, a.start))
    return out


KMER_COLS = ["contig", "start", "end", "period", "tandem_frac", "period2", "period3", "unit"]


def write_karrays(arrs: list[KArray], path) -> None:
    with open(path, "w") as o:
        o.write("\t".join(KMER_COLS) + "\n")
        for a in arrs:
            o.write(f"{a.contig}\t{a.start}\t{a.end}\t{a.period}\t{a.tandem_frac}\t"
                    f"{a.period2}\t{a.period3}\t{a.unit}\n")


def read_karrays(path) -> list[KArray]:
    out = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append(KArray(f[0], int(f[1]), int(f[2]), int(f[3]), float(f[4]),
                              int(f[5]), int(f[6]), f[7] if len(f) > 7 else ""))
    return out


# ---------------------------------------------------------------- short-pass prefilter
def _near_recurrence(codes: np.ndarray, valid: np.ndarray, maxd: int) -> np.ndarray:
    """True where the same k-mer occurs again within `maxd` bp (before or after)."""
    f = np.zeros(len(codes), bool)
    idx = np.flatnonzero(valid)
    if len(idx) < 2:
        return f
    order = np.argsort(codes[idx], kind="stable")
    sp = idx[order]
    near = (codes[idx][order][1:] == codes[idx][order][:-1]) & (sp[1:] - sp[:-1] <= maxd)
    f[sp[1:][near]] = True
    f[sp[:-1][near]] = True
    return f


def short_candidates(seq: str, k: int = 10, maxd: int = 100, chunk: int = 4_000_000,
                     gap: int = 50, min_count: int = 5, min_span: int = 30, pad: int = 150,
                     merge: int = 300) -> list[tuple[int, int]]:
    """Windows of `seq` where a tandem repeat with period <= `maxd` can be: runs of
    positions whose k-mer recurs within `maxd` bp (padded and merged). Only these are
    given to the ULTRA short pass. This typically keeps a few % of the sequence and nearly every
    ULTRA array >= 500 bp (what is lost are weak, low-complexity calls of a few hundred bp)."""
    n = len(seq)
    flag = np.zeros(n, bool)
    for c0 in range(0, n, chunk):
        sub = seq[c0:c0 + chunk + maxd + k]
        codes, valid = kmer_codes(sub, k)
        if len(codes) == 0:
            continue
        f = _near_recurrence(codes, valid, maxd)
        m = min(len(f), n - c0)
        flag[c0:c0 + m] |= f[:m]
    idx = np.flatnonzero(flag)
    if len(idx) == 0:
        return []
    br = np.flatnonzero(np.diff(idx) > gap)
    st = np.r_[idx[0], idx[br + 1]]
    en = np.r_[idx[br], idx[-1]] + k
    cnt = np.diff(np.r_[0, br + 1, len(idx)])
    keep = (en - st >= min_span) & (cnt >= min_count)
    out: list[list[int]] = []
    for a, b in zip(st[keep] - pad, en[keep] + pad):
        a, b = max(0, int(a)), min(n, int(b))
        if out and a <= out[-1][1] + merge:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def _short_worker(job):
    name, offset, seq, kw = job
    return name, [(offset + s, offset + e) for s, e in short_candidates(seq, **kw)]


def run_short_candidates(fasta_iter, threads: int = 1, piece: int = 20_000_000,
                         **kw) -> dict[str, list[tuple[int, int]]]:
    """short_candidates for a whole assembly, contigs cut into `piece` bp jobs
    (overlapping by `maxd + k` so no recurrence across a cut is lost)."""
    over = kw.get("maxd", 100) + kw.get("k", 10)

    def jobs():
        for name, seq in fasta_iter:
            for s in range(0, len(seq), piece):
                yield name, s, seq[s:s + piece + over], kw
    res: dict[str, list[tuple[int, int]]] = {}
    if threads <= 1:
        it = map(_short_worker, jobs())
        for name, iv in it:
            res.setdefault(name, []).extend(iv)
    else:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(threads) as pool:
            for name, iv in pool.imap(_short_worker, jobs(), chunksize=1):
                res.setdefault(name, []).extend(iv)
    merge = kw.get("merge", 300)
    out = {}
    for name, iv in res.items():
        m: list[list[int]] = []
        for a, b in sorted(iv):
            if m and a <= m[-1][1] + merge:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        out[name] = [(a, b) for a, b in m]
    return out


# ---------------------------------------------------------------- short-pass prefilter: tantan
def merge_runs(iv: list[tuple[int, int]], length: int, gap: int = 20, min_span: int = 40,
               pad: int = 100, merge: int = 300) -> list[tuple[int, int]]:
    """tantan BED runs -> padded, merged ULTRA windows."""
    m: list[list[int]] = []
    for s, e in sorted(iv):
        if m and s <= m[-1][1] + gap:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    out: list[list[int]] = []
    for s, e in m:
        if e - s < min_span:
            continue
        s, e = max(0, s - pad), min(length, e + pad)
        if out and s <= out[-1][1] + merge:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def _tantan_worker(job):
    import shlex
    import subprocess
    name, seq, cmd, maxd = job
    r = subprocess.run(shlex.split(cmd) + ["-f", "3", "-w", str(maxd), "-"],
                       input=f">{name}\n{seq}\n".encode(), capture_output=True, check=True)
    iv = []
    for line in r.stdout.decode().splitlines():
        f = line.split("\t")
        iv.append((int(f[1]), int(f[2])))
    return name, merge_runs(iv, len(seq))


def run_tantan_candidates(fasta_iter, tantan_cmd: str = "tantan", threads: int = 1,
                          maxd: int = 100) -> dict[str, list[tuple[int, int]]]:
    """ULTRA short-pass windows from tantan (-f 3 BED, -w maxd), one tantan per contig.
    tantan is a fast HMM tandem-repeat masker; unlike exact k-mer recurrence it also
    flags diverged arrays (typically < 10% of a mammalian genome is kept, with nearly all ULTRA array bp >= 500 bp)."""
    jobs = ((n, s, tantan_cmd, maxd) for n, s in fasta_iter)
    if threads <= 1:
        return dict(map(_tantan_worker, jobs))
    import multiprocessing as mp
    with mp.get_context("fork").Pool(threads) as pool:
        return dict(pool.imap_unordered(_tantan_worker, jobs, chunksize=1))
