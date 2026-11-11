"""Synthetic demo assembly for repeatviz (and the tests).

    python3 demo/make_demo.py demo/      ->  demo.fa, demo.edits.tsv, demo.agp, demo.depth.bed

ctgA 400 kb : telomere (TTAGGG)x200 at the start, 60 kb of a 700-bp satellite at the end,
              a 60 kb array (20 copies) of a 3,000-bp unit (beyond --p-long 2000 -> alert) in the middle
ctgB 300 kb : 40 kb of the same 700-bp satellite at its start, a 20 kb segment duplicated
              from ctgA (low-copy / SD), (CA)n microsatellite
ctgC  80 kb : unplaced, mostly a 171-bp satellite
layout: chr1 = ctgA(+) N100 ctgB(+)  ->  the two satellite blocks meet at the junction
depth : 40x, the satellite windows at 120x (collapse)
"""
import random
import sys
from pathlib import Path

random.seed(7)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    return "".join(random.choice("ACGT") if random.random() < rate else c for c in s)


def array(unit, total, rate=0.01):
    return "".join(mutate(unit, rate) for _ in range(total // len(unit) + 1))[:total]


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    sat700, sat171, long3k = rnd(700), rnd(171), rnd(3000)
    sd = rnd(20000)
    a = "TTAGGG" * 200 + rnd(100000) + sd + rnd(80000) + array(long3k, 60000, 0.003) + rnd(78800) + array(sat700, 60000)
    b = array(sat700, 40000) + rnd(100000) + "CA" * 300 + rnd(60000) + mutate(sd, 0.002) + rnd(79400)
    c = rnd(5000) + array(sat171, 70000) + rnd(5000)
    seqs = {"ctgA": a, "ctgB": b, "ctgC": c}
    with open(out / "demo.fa", "w") as o:
        for n, s in seqs.items():
            o.write(f">{n}\n")
            for i in range(0, len(s), 80):
                o.write(s[i:i + 80] + "\n")
    L = {n: len(s) for n, s in seqs.items()}
    with open(out / "demo.edits.tsv", "w") as o:
        o.write("# demo layout\n")
        o.write(f"chr1\t1\tctgA\t1\t{L['ctgA']}\t+\t100\tdraft\t\n")
        o.write(f"chr1\t2\tctgB\t1\t{L['ctgB']}\t+\t\tdraft\t\n")
        o.write(f"UNPLACED\t1\tctgC\t1\t{L['ctgC']}\t+\t\tdraft\t\n")
    with open(out / "demo.agp", "w") as o:
        o.write(f"chr1\t1\t{L['ctgA']}\t1\tW\tctgA\t1\t{L['ctgA']}\t+\n")
        s = L["ctgA"] + 1
        o.write(f"chr1\t{s}\t{s + 99}\t2\tN\t100\tscaffold\tyes\tna\n")
        o.write(f"chr1\t{s + 100}\t{s + 99 + L['ctgB']}\t3\tW\tctgB\t1\t{L['ctgB']}\t+\n")
        o.write(f"ctgC\t1\t{L['ctgC']}\t1\tW\tctgC\t1\t{L['ctgC']}\t+\n")
    sat = {"ctgA": [(L["ctgA"] - 60000, L["ctgA"])], "ctgB": [(0, 40000)], "ctgC": [(5000, 75000)]}
    with open(out / "demo.depth.bed", "w") as o:
        for n, Ln in L.items():
            for s in range(0, Ln, 10000):
                e = min(Ln, s + 10000)
                d = 120 if any(s < se and e > ss for ss, se in sat[n]) else 40
                o.write(f"{n}\t{s}\t{e}\t{d}\n")
    print(f"wrote {out}/demo.fa ({sum(L.values())} bp), demo.edits.tsv, demo.agp, demo.depth.bed")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent)
