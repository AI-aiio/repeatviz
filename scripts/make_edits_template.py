#!/usr/bin/env python3
"""Write a starting edits.tsv for repeatviz.

  python3 scripts/make_edits_template.py assembly.fa(.fai)            > edits.tsv   # every contig UNPLACED
  python3 scripts/make_edits_template.py assembly.fa --agp scaf.agp    > edits.tsv   # start from an AGP

Contig lengths come from <assembly>.fai if it exists, otherwise from the FASTA (plain or .gz).
With --agp, every W line becomes a row (status 'draft'), N/U gap lines become gap_after of the
previous row. AGP objects whose name does not look like a chromosome (--chrom-regex, same default
as repeatviz plot) become UNPLACED rows, and contigs not used in the AGP are added as UNPLACED. Edit the file afterwards
(see docs/edits_template.tsv for the column rules).
"""
import argparse
import gzip
import os
import re
import sys

HEADER = "object\tpart\tcontig\tstart\tend\torient\tgap_after\tstatus\tnote"


def lengths(path):
    fai = path if path.endswith(".fai") else path + ".fai"
    if os.path.exists(fai):
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
    out, name, n = {}, None, 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        out[name] = n
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("assembly", help="assembly FASTA or its .fai")
    ap.add_argument("--agp", help="optional AGP to start from")
    ap.add_argument("--chrom-regex", default=r"^(scaf_)?(chr)?(\d+|[XYZWxyzw])$",
                    help="AGP objects matching this are chromosomes [same default as repeatviz plot]")
    a = ap.parse_args()
    L = lengths(a.assembly)
    rows, used = [], set()
    rx = re.compile(a.chrom_regex, re.I)
    if a.agp:
        for line in open(a.agp):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[4] in ("N", "U"):
                if rows and rows[-1][0] == f[0] and rows[-1][0] != "UNPLACED":
                    rows[-1][6] = f[5]
                continue
            if not rx.search(f[0]):
                rows.append(["UNPLACED", "1", f[5], f[6], f[7], f[8], "", "draft", f"AGP object {f[0]}"])
                used.add(f[5])
                continue
            obj, part = f[0], sum(1 for r in rows if r[0] == f[0]) + 1
            rows.append([obj, str(part), f[5], f[6], f[7], f[8], "", "draft", ""])
            used.add(f[5])
    w = sys.stdout.write
    w("# made by make_edits_template.py; see docs/edits_template.tsv for the column rules\n" + HEADER + "\n")
    for r in rows:
        w("\t".join(r) + "\n")
    for c, n in L.items():
        if c not in used:
            w(f"UNPLACED\t1\t{c}\t1\t{n}\t+\t\tdraft\t\n")


if __name__ == "__main__":
    main()
