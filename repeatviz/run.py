"""`repeatviz run`: the compute step (meryl + ULTRA). Plain CLI, no scheduler needed.

Steps (each skipped if its output already exists, so a run can be resumed):
  1. contig lengths                         -> lengths.tsv
  2. meryl count k (asm) ; greater-than     -> asm.meryl, ge2.meryl, ge10.meryl
  3. meryl-lookup -bed-runs                 -> multicopy_ge2.bed, multicopy_ge10.bed
  4. ULTRA short pass (-p P_SHORT)          -> ultra_short.tsv
       tantan (default): only windows tantan flags as tandem (period <= P_SHORT) -> shortpass_regions.bed
       kmer: windows where a k-mer recurs within P_SHORT bp (no extra tool, less sensitive)
       full: whole assembly
  5. long, densely high-copy regions        -> longpass_regions.bed
  6. long pass on those regions, by --long-mode:
       kmer  (default, light) k-mer return-distance period per region -> kmer_long.tsv
             + ULTRA on a few sampled windows of big arrays for the unit -> ultra_sampled.p<P>.tsv
       ultra (detailed, slow) ULTRA -p P_LONG on the whole regions     -> ultra_long.tsv
       none  short pass only
  7. telomere scan (TTAGGG density)       -> telomere_runs.tsv
  8. merge / classify                       -> arrays.tsv, windows.tsv.gz, units.fa
"""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from . import telo, kmerper, seqio, tandem, windows


def log(msg: str) -> None:
    print(f"[repeatviz {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def _cmd(prefix: str) -> list[str]:
    return shlex.split(prefix)


def sh(args: list[str], stdout=None, cwd=None) -> None:
    log("$ " + " ".join(shlex.quote(a) for a in args))
    t = time.time()
    subprocess.run(args, check=True, stdout=stdout, cwd=cwd)
    log(f"  done in {time.time() - t:.0f} s")


def ultra_mem_per_thread(ultra: str, period: int, win_size: int | None, overlap: int | None,
                         probe: Path) -> float:
    """GB per thread reported by `ultra --mem` (1 thread)."""
    args = _cmd(ultra) + ["-p", str(period), "-t", "1", "--mem"]
    if win_size:
        args += ["--win_size", str(win_size)]
    if overlap:
        args += ["--overlap", str(overlap)]
    out = subprocess.run(args + [str(probe)], capture_output=True, text=True).stdout
    m = re.search(r"Total size:\s*([\d.]+)\s*GB", out)
    return float(m.group(1)) if m else 8.0


def run_ultra(ultra: str, fasta: Path, out: Path, period: int, threads: int, memory_gb: float,
              min_length: int, win_size: int | None = None, overlap: int | None = None) -> int:
    per = ultra_mem_per_thread(ultra, period, win_size, overlap, fasta)
    t = max(1, min(threads, int(memory_gb * 0.85 // max(per, 0.01))))
    log(f"ULTRA -p {period}: ~{per:.2f} GB/thread -> {t} thread(s) within {memory_gb} GB")
    args = _cmd(ultra) + ["-p", str(period), "-t", str(t), "--tsv", "--hs", "-c",
                          "--min_length", str(min_length), "-o", str(out)]
    if win_size:
        args += ["--win_size", str(win_size)]
    if overlap:
        args += ["--overlap", str(overlap)]
    logp = str(out)[:-len(".partial")] if str(out).endswith(".partial") else str(out)
    with open(logp + ".log", "w") as lg:
        sh(args + [str(fasta)], stdout=lg)
    return t


P_LADDER = (250, 500, 1000, 2000, 3000, 5000)


def sampled_ultra(ultra: str, fasta: Path, kcalls: list[tandem.Call], out: Path, args) -> None:
    """Refine the unit of big k-mer arrays with ULTRA on one window from the middle of each
    array (max(--sample-win, 4 x period) bp), run with -p just above the period (the -p
    ladder is capped at --p-long).

    Only the consensus is taken from ULTRA (when its period agrees with the k-mer period);
    array extent and period stay from the k-mer pass. Arrays whose period exceeds --p-long
    keep the raw middle copy as their unit.
    """
    pick = [c for c in kcalls if c.length >= args.sample_min_array and c.period <= args.p_long]
    if not pick:
        return
    groups: dict[int, list[tuple[tandem.Call, int, int]]] = {}
    for c in pick:
        p = next((x for x in P_LADDER if x >= 1.2 * c.period), args.p_long)
        p = min(p, args.p_long)
        w = min(c.length, max(args.sample_win, 4 * c.period))
        s = c.start + (c.length - w) // 2
        groups.setdefault(p, []).append((c, s, s + w))
    log(f"sampled ULTRA: {len(pick)} arrays >= {args.sample_min_array} bp, "
        f"{sum(e - s for g in groups.values() for _, s, e in g) / 1e6:.1f} Mb in {len(groups)} -p group(s)")
    todo = {p: g for p, g in groups.items() if not (out / f"ultra_sampled.p{p}.tsv").exists()}
    if todo:   # one pass over the FASTA for all -p groups
        regs = {}
        for p, g in todo.items():
            reg: dict[str, list[tuple[int, int]]] = {}
            for c, s, e in g:
                reg.setdefault(c.contig, []).append((s, e))
            regs[p] = {k: windows.merge_intervals(v) for k, v in reg.items()}
        seqio.extract_regions_multi(fasta, {out / f"ultra_sampled.p{p}.fa": r for p, r in regs.items()})
    for p, g in sorted(groups.items()):
        res = out / f"ultra_sampled.p{p}.tsv"
        if not res.exists():
            rfa = out / f"ultra_sampled.p{p}.fa"
            tmp = _fresh_tmp(res)
            win = max(args.long_win_size, 5 * p)
            run_ultra(ultra, rfa, tmp, p, args.threads, args.memory, args.ultra_min_length,
                      win_size=win, overlap=min(2 * p, win // 2))
            tmp.rename(res)
            rfa.unlink()
        calls = tandem.parse_ultra_tsv(res, seqio.region_name_map)
        by_contig: dict[str, list[tandem.Call]] = {}
        for u in calls:
            by_contig.setdefault(u.contig, []).append(u)
        for c, s, e in g:
            hits = [u for u in by_contig.get(c.contig, ()) if u.end > s and u.start < e]
            if not hits:
                continue
            best = max(hits, key=lambda u: min(u.end, e) - max(u.start, s))
            c.ultra_period = best.period
            if abs(best.period - c.period) <= max(2, 0.03 * c.period):
                c.consensus = best.consensus
                c.method = "kmer+ultra"


def write_units(arrays: list[tandem.Call], path: Path, min_len: int) -> int:
    n = 0
    with open(path, "w") as o:
        for c in sorted(arrays, key=lambda x: -x.length):
            if c.length < min_len or c.cls in ("micro",):
                continue
            o.write(f">{c.contig}:{c.start}-{c.end} period={c.period} array_len={c.length} "
                    f"copies={c.copies:.0f} class={c.cls}\n{c.consensus}\n")
            n += 1
    return n


def run(args) -> int:
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    fasta = seqio.decompress_if_needed(Path(args.assembly).resolve(), out)
    meryl, lookup, ultra = args.meryl_cmd, args.meryl_lookup_cmd, args.ultra_cmd
    t0 = time.time()
    params = {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()
              if k != "func"}
    (out / "run_params.json").write_text(json.dumps(params, indent=1))

    # 1. lengths
    lt = out / "lengths.tsv"
    if not lt.exists():
        log("contig lengths")
        seqio.write_lengths(seqio.read_lengths(fasta), lt)
    lengths = seqio.load_lengths(lt)
    log(f"{len(lengths)} sequences, {sum(lengths.values()) / 1e6:.1f} Mb")

    # 2-3. multi-copy k-mers
    th = {"ge2": 1, "ge10": args.high_copy - 1}
    beds_done = all((out / f"multicopy_{n}.bed").exists() for n in th)
    if not (out / "asm.meryl").exists() and not beds_done:   # beds already made -> k-mer DBs not needed
        tmp = _fresh_tmp(out / "asm.meryl")
        sh(_cmd(meryl) + ["count", f"k={args.k}", f"threads={args.threads}",
                          f"memory={int(args.memory)}", str(fasta), "output", str(tmp)])
        tmp.rename(out / "asm.meryl")
    for name, gt in th.items():
        db, bed = out / f"{name}.meryl", out / f"multicopy_{name}.bed"
        if not db.exists() and not bed.exists():
            tmp = _fresh_tmp(db)
            sh(_cmd(meryl) + ["greater-than", str(gt), str(out / "asm.meryl"),
                              "output", str(tmp), f"threads={args.threads}"])
            tmp.rename(db)
        if not bed.exists():
            tmp = _fresh_tmp(bed)
            sh(_cmd(lookup) + ["-bed-runs", "-sequence", str(fasta), "-mers", str(db),
                               "-output", str(tmp), "-threads", str(args.threads)])
            tmp.rename(bed)
    ge2 = windows.read_bed_runs(out / "multicopy_ge2.bed")
    ge10 = windows.read_bed_runs(out / "multicopy_ge10.bed")

    # 4. short pass (ULTRA -p P_SHORT). prefilter (default): only on windows where a k-mer
    #    recurs within P_SHORT bp (~5% of a mammalian genome) instead of the whole assembly
    us = out / "ultra_short.tsv"
    if not us.exists():
        tmp = _fresh_tmp(us)
        if args.short_mode in ("tantan", "kmer"):
            sb = out / "shortpass_regions.bed"
            if not sb.exists():
                t = time.time()
                if args.short_mode == "tantan":
                    log(f"short-pass prefilter: tantan -w {args.p_short} ({args.threads} threads)")
                    sreg = kmerper.run_tantan_candidates(seqio.iter_fasta(fasta), args.tantan_cmd,
                                                         args.threads, args.p_short)
                else:
                    log(f"short-pass prefilter: k-mer recurrence (k={args.prefilter_k}, "
                        f"<= {args.p_short} bp, {args.threads} threads)")
                    sreg = kmerper.run_short_candidates(seqio.iter_fasta(fasta), args.threads,
                                                        k=args.prefilter_k, maxd=args.p_short)
                stmp = _fresh_tmp(sb)
                with open(stmp, "w") as o:
                    for c, iv in sreg.items():
                        for a, b in iv:
                            o.write(f"{c}\t{a}\t{b}\n")
                stmp.rename(sb)
                log(f"  done in {time.time() - t:.0f} s")
            sreg = {}
            for line in open(sb):
                c, a, b = line.split()[:3]
                sreg.setdefault(c, []).append((int(a), int(b)))
            stot = sum(b - a for iv in sreg.values() for a, b in iv)
            log(f"short-pass windows: {sum(map(len, sreg.values()))} ({stot / 1e6:.1f} Mb, "
                f"{100 * stot / max(1, sum(lengths.values())):.1f}% of the assembly)")
            sfa = out / "shortpass_regions.fa"
            seqio.extract_regions(fasta, sreg, sfa)
            run_ultra(ultra, sfa, tmp, args.p_short, args.threads, args.memory, args.ultra_min_length)
            sfa.unlink()
        else:
            run_ultra(ultra, fasta, tmp, args.p_short, args.threads, args.memory, args.ultra_min_length)
        tmp.rename(us)

    # 5. long, dense high-copy regions
    regions = windows.high_copy_regions(ge10, lengths, merge_gap=args.region_merge_gap,
                                        min_len=args.region_min_len, min_frac=args.region_min_frac,
                                        pad=args.region_pad)
    rbed = out / "longpass_regions.bed"
    with open(rbed, "w") as o:
        for c, iv in regions.items():
            for s, e in iv:
                o.write(f"{c}\t{s}\t{e}\n")
    tot = sum(e - s for iv in regions.values() for s, e in iv)
    log(f"long-pass regions: {sum(map(len, regions.values()))} ({tot / 1e6:.1f} Mb), mode {args.long_mode}")
    short = tandem.parse_ultra_tsv(us, seqio.region_name_map)

    # 6. long pass
    if args.long_mode == "ultra":
        ul = out / "ultra_long.tsv"
        if regions and not ul.exists():
            rfa = out / "longpass_regions.fa"
            seqio.extract_regions(fasta, regions, rfa)
            tmp = _fresh_tmp(ul)
            run_ultra(ultra, rfa, tmp, args.p_long, args.threads, args.memory, args.ultra_min_length,
                      win_size=args.long_win_size, overlap=args.long_overlap)
            tmp.rename(ul)
            rfa.unlink()
        long_ = tandem.parse_ultra_tsv(ul, seqio.region_name_map) if ul.exists() else []
        calls = tandem.drop_inside(short, regions) + long_
        arrays = tandem.merge_calls(calls, tol=args.merge_tol)
    elif args.long_mode == "kmer":
        kl = out / "kmer_long.tsv"
        if not kl.exists():
            log(f"k-mer long pass (k={args.kmer_k}, tile {args.kmer_tile} bp, {args.threads} threads)")
            t = time.time()
            karrs = kmerper.run_regions(seqio.iter_fasta(fasta), regions, args.threads,
                                        k=args.kmer_k, tile=args.kmer_tile, tol=args.kmer_tol,
                                        min_frac=args.kmer_min_frac, min_period=args.p_short + 1)
            tmp = _fresh_tmp(kl)
            kmerper.write_karrays(karrs, tmp)
            tmp.rename(kl)
            log(f"  {len(karrs)} arrays, done in {time.time() - t:.0f} s")
        kcalls = [tandem.from_karray(a) for a in kmerper.read_karrays(kl)]
        if args.sample_ultra:
            sampled_ultra(ultra, fasta, kcalls, out, args)
        spans: dict[str, list[tuple[int, int]]] = {}
        for c in kcalls:
            spans.setdefault(c.contig, []).append((c.start, c.end))
        spans = {c: windows.merge_intervals(v) for c, v in spans.items()}
        arrays = tandem.merge_calls(tandem.drop_covered(short, spans), tol=args.merge_tol) + kcalls
    else:
        arrays = tandem.merge_calls(short, tol=args.merge_tol)

    # 7. telomere by motif density (v0.3): exact TTAGGG / CCCTAA runs, like tidk
    tb = out / "telomere_runs.tsv"
    if not tb.exists():
        t = time.time()
        log(f"telomere scan (TTAGGG/CCCTAA runs: gap <= {args.telo_gap} bp, >= {args.telo_min_len} bp, "
            f"density >= {args.telo_density})")
        runs = telo.scan_fasta(seqio.iter_fasta(fasta), threads=args.threads, max_gap=args.telo_gap,
                               min_len=args.telo_min_len, min_density=args.telo_density)
        tmp = _fresh_tmp(tb)
        telo.write_runs(runs, tmp)
        tmp.rename(tb)
        log(f"  {len(runs)} runs, {sum(r[2] - r[1] for r in runs) / 1e6:.2f} Mb, done in {time.time() - t:.0f} s")
    arrays, relab, added = telo.apply(
        arrays, telo.read_runs(tb), cover=0.5,
        make_call=lambda c, s, e, n, st: tandem.Call(c, s, e, 6, 6.0 * n, telo.MOTIFS[st], float(n),
                                                     method="telo_scan", cls_override="telomere"))
    log(f"telomere: {relab} arrays relabelled by motif density, {added} motif runs added as arrays")

    # 8. merge / classify / tables
    arrays = [a for a in arrays if a.length >= args.min_array]
    n = tandem.write_arrays(arrays, out / "arrays.tsv")
    nu = write_units(arrays, out / "units.fa", args.unit_min_array)
    adf = pd.read_csv(out / "arrays.tsv", sep="\t")
    wdf = windows.build_windows(lengths, args.window, ge2, ge10, adf)
    wdf.to_csv(out / "windows.tsv.gz", sep="\t", index=False)
    log(f"{n} arrays (>= {args.min_array} bp), {nu} units in units.fa; windows.tsv.gz "
        f"({len(wdf)} x {args.window} bp). total {time.time() - t0:.0f} s")
    if not args.keep_meryl:
        for name in ("asm", "ge2", "ge10"):
            subprocess.run(["rm", "-rf", str(out / f"{name}.meryl")])
    return 0


def add_parser(sub) -> None:
    env = os.environ.get
    p = sub.add_parser("run", help="compute tracks (needs meryl, meryl-lookup, ULTRA)",
                       description="Compute multi-copy k-mer and tandem-repeat tracks for an assembly.")
    p.add_argument("assembly", help="assembly FASTA (plain or .gz)")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("-m", "--memory", type=float, default=32, help="memory budget in GB "
                   "(meryl count; ULTRA thread count is chosen to fit)")
    p.add_argument("--window", type=int, default=100_000, help="window size for tracks (bp)")
    p.add_argument("-k", type=int, default=21, help="k-mer size")
    p.add_argument("--high-copy", type=int, default=10, help="k-mer count for 'high-copy'")
    p.add_argument("--p-short", type=int, default=100, help="ULTRA max period, genome-wide pass")
    p.add_argument("--short-mode", choices=("tantan", "kmer", "full"), default="tantan",
                   help="ULTRA short pass only on windows flagged by tantan (default; ~8%% of a "
                        "mammalian genome, ~10x faster), by exact k-mer recurrence (no extra tool, "
                        "misses diverged arrays), or on the whole assembly (full; slow)")
    p.add_argument("--prefilter-k", type=int, default=10, help="k for --short-mode kmer")
    p.add_argument("--long-mode", choices=("kmer", "ultra", "none"), default="kmer",
                   help="long pass on high-copy regions: kmer = light (k-mer return distances, "
                        "minutes; default), ultra = detailed ULTRA -p P_LONG (hours), none = skip")
    p.add_argument("--p-long", type=int, default=2000,
                   help="ULTRA max period: whole regions (--long-mode ultra) / sampled windows (kmer)")
    p.add_argument("--kmer-k", type=int, default=16, help="k for the k-mer long pass (<= 31)")
    p.add_argument("--kmer-tile", type=int, default=5000, help="tile size of the k-mer long pass (bp)")
    p.add_argument("--kmer-tol", type=float, default=0.03, help="relative period tolerance (k-mer pass)")
    p.add_argument("--kmer-min-frac", type=float, default=0.3,
                   help="min share of a tile explained by its period to call it tandem")
    p.add_argument("--no-sample-ultra", dest="sample_ultra", action="store_false",
                   help="kmer mode: skip ULTRA on sampled windows (units = raw middle copy)")
    p.add_argument("--sample-min-array", type=int, default=50000,
                   help="kmer mode: arrays at least this long get a sampled-ULTRA unit")
    p.add_argument("--sample-win", type=int, default=10000,
                   help="kmer mode: sampled window length (at least 4 x period)")
    p.add_argument("--long-win-size", type=int, default=10000,
                   help="ULTRA --win_size for the long pass (smaller = less memory)")
    p.add_argument("--long-overlap", type=int, default=4000,
                   help="ULTRA --overlap for the long pass (>= 2x p-long recommended)")
    p.add_argument("--ultra-min-length", type=int, default=50, help="ULTRA --min_length")
    p.add_argument("--region-min-len", type=int, default=10000)
    p.add_argument("--region-min-frac", type=float, default=0.8)
    p.add_argument("--region-merge-gap", type=int, default=1000)
    p.add_argument("--region-pad", type=int, default=1000)
    p.add_argument("--merge-tol", type=float, default=0.02, help="period tolerance when merging")
    p.add_argument("--min-array", type=int, default=100, help="min array length kept (bp)")
    p.add_argument("--unit-min-array", type=int, default=10000,
                   help="arrays at least this long get their unit written to units.fa")
    p.add_argument("--telo-gap", type=int, default=18,
                   help="telomere scan: chain TTAGGG/CCCTAA hits closer than this (bp; 18 = up to 2 variant units) [18]")
    p.add_argument("--telo-min-len", type=int, default=100, help="telomere scan: minimum run length (bp) [100]")
    p.add_argument("--telo-density", type=float, default=0.5,
                   help="telomere scan: minimum share of the run that is exact motif [0.5]")
    p.add_argument("--keep-meryl", action="store_true", help="keep the meryl databases")
    p.add_argument("--meryl-cmd", default=env("REPEATVIZ_MERYL", "meryl"),
                   help="command (or container prefix) for meryl [env REPEATVIZ_MERYL]")
    p.add_argument("--meryl-lookup-cmd", default=env("REPEATVIZ_MERYL_LOOKUP", "meryl-lookup"),
                   help="[env REPEATVIZ_MERYL_LOOKUP]")
    p.add_argument("--tantan-cmd", default=env("REPEATVIZ_TANTAN", "tantan"), help="[env REPEATVIZ_TANTAN]")
    p.add_argument("--ultra-cmd", default=env("REPEATVIZ_ULTRA", "ultra"), help="[env REPEATVIZ_ULTRA]")
    p.set_defaults(func=run)


def _fresh_tmp(final):
    """Path for writing a step's output before renaming it to `final`, so an interrupted step
    (timeout, kill) leaves no file that a resumed run would mistake for finished output."""
    import shutil
    tmp = final.with_name(final.name + ".partial")
    if tmp.is_dir():
        shutil.rmtree(tmp)
    elif tmp.exists():
        tmp.unlink()
    return tmp
