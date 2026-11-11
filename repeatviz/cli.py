"""Command-line interface: `repeatviz run` (compute) and `repeatviz plot` (draw)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import __version__, run as run_mod


def plot(args) -> int:
    from . import layout as lay, render, summary
    from .seqio import load_lengths

    rd = Path(args.rundir)
    lengths = load_lengths(rd / "lengths.tsv")
    win = pd.read_csv(rd / "windows.tsv.gz", sep="\t")
    arrays = pd.read_csv(rd / "arrays.tsv", sep="\t")
    if args.agp:
        layout = lay.from_agp(Path(args.agp), lengths, chrom_regex=args.chrom_regex)
    elif args.edits:
        layout = lay.from_edits(Path(args.edits), lengths)
    else:  # no layout: every contig >= --min-chrom-len is its own "chromosome"
        chroms = [lay.Chrom(c, [lay.Part(c, 0, L, "+", 0)], L) for c, L in lengths.items()
                  if L >= args.min_chrom_len]
        layout = lay.Layout(sorted(chroms, key=lambda c: lay.natural_key(c.name)),
                            [(c, L) for c, L in lengths.items() if L < args.min_chrom_len], [])
    has_depth = bool(args.depth)
    if has_depth:
        win, med = summary.add_depth(win, summary.read_depth(Path(args.depth)),
                                     args.collapse, args.low_depth)
        print(f"depth median {med:.1f}x", file=sys.stderr)
    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    al = summary.alerts(win, args.alert_high_copy, args.alert_max_tandem)
    summary.project(al, layout).to_csv(f"{out}.alerts.tsv", sep="\t", index=False)
    summary.project(arrays, layout).to_csv(f"{out}.arrays.layout.tsv", sep="\t", index=False)
    ce = summary.contig_ends(win, arrays, layout, lengths, args.end_flank, args.end_high_copy,
                             args.end_low_copy, args.end_tandem, args.telo_min_bp, args.telo_tip)
    ce.to_csv(f"{out}.contig_ends.tsv", sep="\t", index=False)
    if has_depth:
        win.to_csv(f"{out}.windows_depth.tsv.gz", sep="\t", index=False)
    files = render.render(layout, win, arrays, al, out, args.title or rd.name, has_depth,
                          args.collapse, per_page=args.per_page, width_in=args.width,
                          height_in=args.height, show_alerts=args.alerts,
                          maps=[m.strip() for m in args.maps.split(",") if m.strip()],
                          depth_warn=args.depth_warn)
    print("wrote " + ", ".join(str(f) for f in files)
          + f", {out}.contig_ends.tsv, {out}.alerts.tsv, {out}.arrays.layout.tsv", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="repeatviz", description=(
        "Where are the repeats in a draft assembly, and are they where contigs break? "
        "`run` computes tracks with meryl + ULTRA; `plot` draws them on a chromosome layout."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=f"repeatviz {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    run_mod.add_parser(sub)
    q = sub.add_parser("plot", help="draw tracks on a chromosome layout (pure Python)",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    q.add_argument("rundir", help="output directory of `repeatviz run`")
    q.add_argument("-o", "--out-prefix", required=True)
    g = q.add_mutually_exclusive_group()
    g.add_argument("--agp", help="AGP layout (e.g. synscaf / YaHS output)")
    g.add_argument("--edits", help="curation edits.tsv / edits.draft.tsv")
    q.add_argument("--chrom-regex", default=r"^(scaf_)?(chr)?(\d+|[XYZWxyzw])$",
                   help="AGP objects matching this are chromosomes, others unplaced")
    q.add_argument("--min-chrom-len", type=int, default=10_000_000,
                   help="without a layout, contigs at least this long are drawn as rows")
    q.add_argument("--depth", help="read depth per window (mosdepth regions.bed.gz or bedGraph)")
    q.add_argument("--collapse", type=float, default=1.5, help="rel. depth >= this = collapse")
    q.add_argument("--low-depth", type=float, default=0.3, help="rel. depth < this = low")
    q.add_argument("--alert-high-copy", type=float, default=0.5)
    q.add_argument("--alerts", action="store_true", help="also draw alert markers on the figure (alerts.tsv is always written)")
    q.add_argument("--alert-max-tandem", type=float, default=0.2)
    q.add_argument("--end-flank", type=int, default=200_000, help="bp examined at each contig end")
    q.add_argument("--end-high-copy", type=float, default=0.3)
    q.add_argument("--end-low-copy", type=float, default=0.5)
    q.add_argument("--end-tandem", type=float, default=0.2)
    q.add_argument("--telo-min-bp", type=int, default=300)
    q.add_argument("--telo-tip", type=int, default=20_000, help="telomere must lie this close to the tip")
    q.add_argument("--per-page", type=int, default=15, help="chromosomes per page")
    q.add_argument("--width", type=float, default=11, help="page width (inches; portrait)")
    q.add_argument("--height", type=float, default=17, help="page height (inches)")
    q.add_argument("--maps", default="all,tandem,long,dup",
                   help="which repeat maps to draw, one PDF each: all, tandem (unit <100 bp), "
                        "long (unit >=100 bp), dup (2-9 copy k-mers = segmental duplications)")
    q.add_argument("--depth-warn", type=float, default=2.0,
                   help="yellow warning mark where relative depth >= this")
    q.add_argument("--title", default=None)
    q.set_defaults(func=plot)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
