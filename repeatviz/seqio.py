"""Minimal FASTA helpers (streaming; no samtools needed)."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator


def _open(path: Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (name, sequence) one record at a time (name = first word)."""
    name, buf = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def read_lengths(path: Path) -> dict[str, int]:
    """Contig lengths, from `<fasta>.fai` if present, else by streaming the FASTA."""
    fai = Path(str(path) + ".fai")
    if fai.exists() and Path(path).suffix != ".gz":
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai) if l.strip()}
    return {n: len(s) for n, s in iter_fasta(path)}


def write_lengths(lengths: dict[str, int], path: Path) -> None:
    with open(path, "w") as o:
        for n, L in lengths.items():
            o.write(f"{n}\t{L}\n")


def load_lengths(path: Path) -> dict[str, int]:
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(path) if l.strip()}


SEP = "__rvoff__"


def extract_regions(fasta: Path, regions: dict[str, list[tuple[int, int]]], out: Path) -> int:
    """Write sub-sequences (0-based half-open) named `<contig>__rvoff__<start>`."""
    n = 0
    with open(out, "w") as o:
        for name, seq in iter_fasta(fasta):
            for s, e in regions.get(name, ()):
                o.write(f">{name}{SEP}{s}\n")
                sub = seq[s:e]
                for i in range(0, len(sub), 80):
                    o.write(sub[i:i + 80] + "\n")
                n += 1
    return n


def extract_regions_multi(fasta: Path, jobs: dict[Path, dict[str, list[tuple[int, int]]]]) -> None:
    """Like extract_regions for several outputs, reading the FASTA once."""
    fhs = {out: open(out, "w") for out in jobs}
    try:
        for name, seq in iter_fasta(fasta):
            for out, regions in jobs.items():
                for s, e in regions.get(name, ()):
                    o = fhs[out]
                    o.write(f">{name}{SEP}{s}\n")
                    sub = seq[s:e]
                    for i in range(0, len(sub), 80):
                        o.write(sub[i:i + 80] + "\n")
    finally:
        for o in fhs.values():
            o.close()


def region_name_map(seqid: str) -> tuple[str, int]:
    if SEP in seqid:
        c, off = seqid.rsplit(SEP, 1)
        return c, int(off)
    return seqid, 0


def decompress_if_needed(fasta: Path, workdir: Path) -> Path:
    fasta = Path(fasta)
    if fasta.suffix != ".gz":
        return fasta
    out = workdir / fasta.with_suffix("").name
    if not out.exists():
        with gzip.open(fasta, "rt") as i, open(out, "w") as o:
            for line in i:
                o.write(line)
    return out
