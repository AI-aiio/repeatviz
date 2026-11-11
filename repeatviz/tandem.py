"""Parse ULTRA tandem-repeat calls, merge window-split calls into arrays, classify.

ULTRA (--tsv) coordinates are 0-based, end-exclusive (verified: a 2,000 bp array
inserted at 1000 is reported as 1000..3005 incl. a few bases of extension).
ULTRA analyses the sequence in windows (~25-50 kb), so one long array is
reported as several adjacent calls; `merge_calls` joins them back.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# unit-length classes, in drawing order
CLASSES = ("telomere", "micro", "mini", "satellite", "macro")
CLASS_LABEL = {
    "telomere": "telomere (TTAGGG)",
    "micro": "unit <10 bp",
    "mini": "unit 10-99 bp",
    "satellite": "unit 100-999 bp",
    "macro": "unit >=1000 bp",
}

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
_TELO = {"TTAGGG", "CCCTAA"}
_TELO_ROT = {m[i:] + m[:i] for m in _TELO for i in range(6)}


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def is_telomere_unit(unit: str) -> bool:
    """True if `unit` is (a rotation of) TTAGGG / CCCTAA, or 2 copies of it."""
    u = unit.upper()
    if len(u) == 12 and u[:6] == u[6:]:
        u = u[:6]
    return len(u) == 6 and u in _TELO_ROT


def classify(period: int, unit: str) -> str:
    if is_telomere_unit(unit):
        return "telomere"
    if period < 10:
        return "micro"
    if period < 100:
        return "mini"
    if period < 1000:
        return "satellite"
    return "macro"


@dataclass
class Call:
    contig: str
    start: int  # 0-based
    end: int    # exclusive
    period: int
    score: float
    consensus: str
    copies: float = 0.0
    members: int = 1
    best_len: int = field(default=0)
    method: str = "ultra"          # ultra | kmer | kmer+ultra
    tandem_frac: float = float("nan")
    period2: int = 0
    period3: int = 0
    ultra_period: int = 0
    cls_override: str = ""         # v0.3: "telomere" when set by motif density (telo.apply)

    def __post_init__(self):
        if not self.best_len:
            self.best_len = self.end - self.start

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def cls(self) -> str:
        return self.cls_override or classify(self.period, self.consensus)


def parse_ultra_tsv(path: Path, name_map=None) -> list[Call]:
    """Read ULTRA --tsv output (settings/summary lines are skipped).

    name_map: optional function seqid -> (contig, offset) for calls made on
    extracted sub-sequences.
    """
    calls: list[Call] = []
    cols = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("SeqID\t"):
                cols = line.rstrip("\n").split("\t")
                continue
            if cols is None or not line.strip() or line.startswith(("-", "=", "*")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            d = dict(zip(cols, f))
            try:
                s, e, p = int(d["Start"]), int(d["End"]), int(d["Period"])
            except (KeyError, ValueError):
                continue
            contig, off = (name_map(d["SeqID"]) if name_map else (d["SeqID"], 0))
            cp = float(d.get("#copies", 0) or 0) or (e - s) / max(p, 1)
            calls.append(Call(contig, s + off, e + off, p, float(d["Score"]), d["Consensus"], cp))
    return calls


def _similar(p1: int, p2: int, tol: float) -> bool:
    return abs(p1 - p2) <= max(2, tol * max(p1, p2))


def merge_calls(calls: list[Call], tol: float = 0.02, min_gap: int = 50) -> list[Call]:
    """Join adjacent/overlapping calls of similar period on the same contig.

    Two calls merge if the next starts within max(2*period, min_gap) bp of the
    current end and the periods differ by <= max(2 bp, tol*period). The merged
    array keeps the consensus of its longest member.
    """
    out: list[Call] = []
    for c in sorted(calls, key=lambda x: (x.contig, x.start, x.end)):
        if out:
            m = out[-1]
            if (m.contig == c.contig and _similar(m.period, c.period, tol)
                    and c.start <= m.end + max(2 * m.period, min_gap)):
                if c.length > m.best_len:
                    m.consensus, m.period, m.best_len = c.consensus, c.period, c.length
                m.copies += c.copies * max(0, c.end - max(c.start, m.end)) / max(c.length, 1)
                m.end = max(m.end, c.end)
                m.score += c.score
                m.members += 1
                continue
        out.append(Call(c.contig, c.start, c.end, c.period, c.score, c.consensus,
                        c.copies, 1, c.length))
    return out


def from_karray(a) -> Call:
    """kmerper.KArray -> Call (score = bp explained by the period)."""
    return Call(a.contig, a.start, a.end, a.period, a.tandem_frac * a.length, a.unit,
                a.length / a.period, 1, a.length, method="kmer", tandem_frac=a.tandem_frac,
                period2=a.period2, period3=a.period3)


def drop_inside(calls: list[Call], regions: dict[str, list[tuple[int, int]]]) -> list[Call]:
    """Remove calls whose midpoint lies inside a scanned region (the long pass
    re-annotates those regions with a larger maximum period)."""
    if not regions:
        return calls
    keep = []
    for c in calls:
        mid = (c.start + c.end) // 2
        if not any(s <= mid < e for s, e in regions.get(c.contig, ())):
            keep.append(c)
    return keep


def drop_covered(calls: list[Call], spans: dict[str, list[tuple[int, int]]],
                 frac: float = 0.8) -> list[Call]:
    """Remove calls that are >= `frac` covered by `spans` (merged, sorted intervals per contig).
    Used in the k-mer long mode: a short-period call that runs on beyond a k-mer array (e.g. the
    short monomers of a higher-order repeat, called over a long stretch by the short pass) is kept."""
    import bisect
    keep = []
    for c in calls:
        iv = spans.get(c.contig, ())
        i = max(0, bisect.bisect_right(iv, (c.start, float("inf"))) - 1)
        cov = 0
        while i < len(iv) and iv[i][0] < c.end:
            cov += max(0, min(iv[i][1], c.end) - max(iv[i][0], c.start))
            i += 1
        if cov < frac * c.length:
            keep.append(c)
    return keep


# v0.2: columns after `consensus` are new (older readers that index by position still work)
ARRAY_COLS = ["contig", "start", "end", "length", "period", "copies", "class",
              "score", "n_calls", "consensus", "method", "tandem_frac", "period2", "period3",
              "ultra_period"]


def write_arrays(calls: list[Call], path: Path, min_len: int = 0) -> int:
    n = 0
    with open(path, "w") as o:
        o.write("\t".join(ARRAY_COLS) + "\n")
        for c in sorted(calls, key=lambda x: (x.contig, x.start)):
            if c.length < min_len:
                continue
            o.write(f"{c.contig}\t{c.start}\t{c.end}\t{c.length}\t{c.period}\t{c.copies:.1f}\t"
                    f"{c.cls}\t{c.score:.1f}\t{c.members}\t{c.consensus}\t{c.method}\t"
                    f"{'' if c.tandem_frac != c.tandem_frac else c.tandem_frac}\t"
                    f"{c.period2 or ''}\t{c.period3 or ''}\t{c.ultra_period or ''}\n")
            n += 1
    return n
