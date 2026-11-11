"""Chromosome layouts from a synscaf/YaHS-style AGP or a curation edits.tsv.

Internally every placed piece is a Part with 0-based half-open contig
coordinates and its position on the chromosome (also 0-based half-open).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Part:
    contig: str
    cstart: int  # 0-based
    cend: int    # exclusive
    orient: str  # + / -
    pos: int     # chromosome start (0-based)

    @property
    def length(self) -> int:
        return self.cend - self.cstart

    def to_chrom(self, x: int) -> int | None:
        """Contig coordinate -> chromosome coordinate (None if outside this part)."""
        if not (self.cstart <= x <= self.cend):
            return None
        return self.pos + (x - self.cstart if self.orient == "+" else self.cend - x)


@dataclass
class Chrom:
    name: str
    parts: list[Part] = field(default_factory=list)
    length: int = 0
    status: str = ""


@dataclass
class Layout:
    chroms: list[Chrom]
    unplaced: list[tuple[str, int]]  # (contig, length)
    removed: list[tuple[str, int]]


def natural_key(name: str):
    n = re.sub(r"^(scaf_|scaffold_|super_)", "", name, flags=re.I)
    m = re.match(r"^(chr)?(\d+)([A-Za-z]?)$", n, flags=re.I)
    if m:
        return (0, int(m.group(2)), m.group(3))
    m = re.match(r"^(chr)?([A-Za-z]+)$", n, flags=re.I)
    if m:
        order = {"X": 0, "Y": 1, "Z": 2, "W": 3, "M": 9, "MT": 9}
        return (1, order.get(m.group(2).upper(), 5), m.group(2))
    return (2, 0, n)


def display_name(name: str) -> str:
    return re.sub(r"^(scaf_)", "", name)


def from_edits(path: Path, lengths: dict[str, int]) -> Layout:
    """edits.tsv: object part contig start end orient gap_after status note (1-based incl.)."""
    chroms: dict[str, Chrom] = {}
    unplaced, removed = [], []
    rows = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "object" or len(f) < 7:
            continue
        rows.append(f)
    rows.sort(key=lambda f: (f[0], int(f[1]) if f[1].isdigit() else 0))
    for f in rows:
        obj, contig, s, e, o, gap = f[0], f[2], int(f[3]), int(f[4]), f[5], f[6]
        status = f[7] if len(f) > 7 else ""
        if obj == "UNPLACED":
            unplaced.append((contig, e - s + 1))
            continue
        if obj == "REMOVED":
            removed.append((contig, e - s + 1))
            continue
        ch = chroms.setdefault(obj, Chrom(obj))
        ch.parts.append(Part(contig, s - 1, e, o, ch.length))
        ch.length += e - s + 1 + (int(gap) if gap.strip().isdigit() else 0)
        if status == "draft" or not ch.status:
            ch.status = status
    for ch in chroms.values():  # last part has no trailing gap
        last = ch.parts[-1]
        ch.length = last.pos + last.length
    placed = {p.contig for c in chroms.values() for p in c.parts} | {u for u, _ in unplaced + removed}
    unplaced += [(c, L) for c, L in lengths.items() if c not in placed]
    return Layout(sorted(chroms.values(), key=lambda c: natural_key(c.name)), unplaced, removed)


def from_agp(path: Path, lengths: dict[str, int], min_parts_chrom: int = 1,
             chrom_regex: str = r"^(scaf_)?(chr)?(\d+|[XYZWxyzw])$") -> Layout:
    """AGP v2. Objects whose name matches `chrom_regex` are chromosomes; other
    objects made of a single component named like the object are 'unplaced'."""
    objs: dict[str, Chrom] = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        obj, ob, oe, typ = f[0], int(f[1]), int(f[2]), f[4]
        ch = objs.setdefault(obj, Chrom(obj))
        ch.length = max(ch.length, oe)
        if typ in ("N", "U"):
            continue
        ch.parts.append(Part(f[5], int(f[6]) - 1, int(f[7]), f[8] if f[8] in "+-" else "+", ob - 1))
    rx = re.compile(chrom_regex, re.I)
    chroms, unplaced = [], []
    for name, ch in objs.items():
        if rx.match(name) and len(ch.parts) >= min_parts_chrom:
            chroms.append(ch)
        else:
            unplaced += [(p.contig, p.length) for p in ch.parts]
    placed = {p.contig for c in objs.values() for p in c.parts}
    unplaced += [(c, L) for c, L in lengths.items() if c not in placed]
    return Layout(sorted(chroms, key=lambda c: natural_key(c.name)), unplaced, [])


def contig_index(layout: Layout) -> dict[str, list[tuple[Chrom, Part]]]:
    idx: dict[str, list[tuple[Chrom, Part]]] = {}
    for ch in layout.chroms:
        for p in ch.parts:
            idx.setdefault(p.contig, []).append((ch, p))
    return idx
