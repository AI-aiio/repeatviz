from pathlib import Path

import numpy as np
import pandas as pd

from repeatviz import layout as lay, summary, tandem, windows
from repeatviz.seqio import SEP, region_name_map

ULTRA_TSV = """----------------------------
Maximum repeat period: 100
----------------------------
SeqID\tStart\tEnd\tPeriod\tScore\tConsensus\t#copies\t#substitutions\t#insertions\t#deletions\t#Subrepeats\tSubrepeatStarts\tSubrepeatConsensi
t\t1000\t3005\t50\t2612.2\tAAAGCG\t40\t2\t0\t0\t1\t0\tAAAGCG
t\t4000\t4600\t6\t791.7\tAGGGTT\t100\t0\t0\t0\t1\t0\tAGGGTT
r{sep}5000\t0\t25000\t700\t100\tACGT\t35\t0\t0\t0\t1\t0\tACGT
r{sep}5000\t24000\t50000\t703\t100\tACGTA\t37\t0\t0\t0\t1\t0\tACGTA
""".replace("{sep}", SEP)


def test_parse_and_offsets(tmp_path):
    p = tmp_path / "u.tsv"
    p.write_text(ULTRA_TSV)
    calls = tandem.parse_ultra_tsv(p, region_name_map)
    assert [(c.contig, c.start, c.end) for c in calls][:2] == [("t", 1000, 3005), ("t", 4000, 4600)]
    assert calls[2].contig == "r" and calls[2].start == 5000 and calls[3].end == 55000


def test_merge_and_classes(tmp_path):
    p = tmp_path / "u.tsv"
    p.write_text(ULTRA_TSV)
    arr = tandem.merge_calls(tandem.parse_ultra_tsv(p, region_name_map))
    r = [a for a in arr if a.contig == "r"]
    assert len(r) == 1 and (r[0].start, r[0].end) == (5000, 55000) and r[0].members == 2
    cls = {(a.contig, a.start): a.cls for a in arr}
    assert cls[("t", 4000)] == "telomere" and cls[("t", 1000)] == "mini" and r[0].cls == "satellite"


def test_telomere_rotations():
    assert tandem.is_telomere_unit("GGGTTA") and tandem.is_telomere_unit("CTAACC")
    assert tandem.is_telomere_unit("TTAGGGTTAGGG") and not tandem.is_telomere_unit("TTAGGC")


def test_merge_keeps_different_periods_apart():
    a = tandem.Call("c", 0, 1000, 171, 1, "A")
    b = tandem.Call("c", 1000, 2000, 340, 1, "B")
    assert len(tandem.merge_calls([a, b])) == 2


def test_binned_coverage_and_regions():
    cov = windows.binned_coverage([(5, 15), (10, 25)], 30, 10)
    assert cov.tolist() == [5, 10, 5]
    runs = {"c": [(0, 9000), (9500, 20000), (50000, 51000)]}
    reg = windows.high_copy_regions(runs, {"c": 60000}, merge_gap=1000, min_len=10000,
                                    min_frac=0.8, pad=100)
    assert reg == {"c": [(0, 20100)]}


def _layout(tmp_path):
    e = tmp_path / "e.tsv"
    e.write_text("# x\nchr1\t1\tA\t1\t1000\t+\t100\tdecided\t\n"
                 "chr1\t2\tB\t101\t600\t-\t\tdraft\t\nUNPLACED\t1\tC\t1\t50\t+\t\tdecided\t\n")
    return lay.from_edits(e, {"A": 1000, "B": 700, "C": 50, "D": 9})


def test_edits_layout_and_projection(tmp_path):
    L = _layout(tmp_path)
    ch = L.chroms[0]
    assert ch.name == "chr1" and ch.length == 1000 + 100 + 500 and ch.status == "draft"
    b = ch.parts[1]
    assert (b.cstart, b.cend, b.pos) == (100, 600, 1100)
    assert b.to_chrom(600) == 1100 and b.to_chrom(100) == 1600
    assert ("C", 50) in L.unplaced and ("D", 9) in L.unplaced
    df = summary.project(pd.DataFrame({"contig": ["B"], "start": [500], "end": [600]}), L)
    assert (df.chrom[0], df.chrom_start[0], df.chrom_end[0]) == ("chr1", 1100, 1200)


def test_agp_layout(tmp_path):
    a = tmp_path / "x.agp"
    a.write_text("scaf_chr2\t1\t10\t1\tW\tA\t1\t10\t+\nscaf_chr2\t11\t20\t2\tN\t10\tscaffold\tyes\tna\n"
                 "scaf_chr2\t21\t25\t3\tW\tB\t1\t5\t-\nB2\t1\t7\t1\tW\tB2\t1\t7\t+\n")
    L = lay.from_agp(a, {"A": 10, "B": 5, "B2": 7})
    assert [c.name for c in L.chroms] == ["scaf_chr2"] and L.chroms[0].parts[1].pos == 20
    assert L.unplaced == [("B2", 7)]


def test_contig_end_roles(tmp_path):
    L = _layout(tmp_path)
    lengths = {"A": 1000, "B": 700, "C": 50, "D": 9}
    arrays = pd.DataFrame([["A", 900, 1000, 100, 150, 0.7, "satellite", 1, 1, "AC"],
                           ["A", 0, 400, 400, 6, 66, "telomere", 1, 1, "TTAGGG"]],
                          columns=tandem.ARRAY_COLS[:10])   # v0.1 layout still accepted
    win = windows.build_windows(lengths, 100, {}, {"A": [(900, 1000)]}, arrays)
    ce = summary.contig_ends(win, arrays, L, lengths, flank=200, hc=0.3, lc=0.5,
                             tandem_min=0.2, telo_min_bp=300, telo_end=500).set_index(["contig", "end"])
    assert ce.loc[("A", "start"), "role"] == "chr_end" and ce.loc[("A", "start"), "category"] == "telomere"
    assert ce.loc[("A", "end"), "role"] == "junction" and "satellite" in ce.loc[("A", "end"), "category"]
    assert ce.loc[("B", "end"), "role"] == "junction"   # B is '-': its end faces A
    assert ce.loc[("B", "start"), "role"] == "chr_end"
    assert ce.loc[("C", "start"), "role"] == "unplaced"


def test_alerts_and_depth():
    win = pd.DataFrame({"contig": ["c"] * 4, "start": [0, 10, 20, 30], "end": [10, 20, 30, 40],
                        "frac_ge2": [0, 1, 1, 0], "frac_ge10": [0, 0.9, 0.9, 0], "tandem_all": [0, 0, 0, 0]})
    al = summary.alerts(win, 0.5, 0.2)
    assert al[["start", "end"]].values.tolist() == [[10, 30]]
    dep = pd.DataFrame({"contig": ["c"] * 4, "start": [0, 10, 20, 30], "end": [10, 20, 30, 40],
                        "depth": [40.0, 40, 100, 5]})
    w2, med = summary.add_depth(win, dep, 1.5, 0.3)
    assert med == 40 and w2.depth_flag.tolist() == ["", "", "collapse", "low"]


def test_render_smoke(tmp_path):
    from repeatviz import render
    L = _layout(tmp_path)
    lengths = {"A": 1000, "B": 700, "C": 50, "D": 9}
    arrays = pd.DataFrame(columns=tandem.ARRAY_COLS)
    win = windows.build_windows(lengths, 100, {}, {}, arrays)
    files = render.render(L, win, arrays, summary.alerts(win, 0.5, 0.2), tmp_path / "o", "t", False, 1.5)
    assert all(Path(f).exists() for f in files) and sum(str(f).endswith(".pdf") and Path(f).parent.name == "repeatviz_pdf" for f in files) == 4  # all / tandem / long / dup in repeatviz_pdf/


# ---------------------------------------------------------------- v0.2 k-mer long pass
def _rnd(n, rng):
    return "".join(rng.choice(list("ACGT"), n))


def _array(unit, total, rate, rng):
    out = []
    while sum(map(len, out)) < total:
        u = np.array(list(unit))
        m = rng.random(len(u)) < rate
        u[m] = rng.choice(list("ACGT"), m.sum())
        out.append("".join(u))
    return "".join(out)[:total]


def test_kmer_period_long_unit_and_satellite():
    from repeatviz import kmerper
    rng = np.random.default_rng(1)
    seq = (_rnd(20000, rng) + _array(_rnd(3000, rng), 60000, 0.003, rng) + _rnd(20000, rng)
           + _array(_rnd(700, rng), 30000, 0.01, rng) + _rnd(10000, rng))
    arr = kmerper.analyse_region("c", 1000, seq)
    assert [a.period for a in arr] == [3000, 700]
    assert abs(arr[0].start - 21000) < 100 and abs(arr[0].end - 81000) < 100
    assert arr[0].tandem_frac > 0.8 and len(arr[0].unit) == 3000


def test_kmer_no_period_limit_and_random_is_empty():
    from repeatviz import kmerper
    rng = np.random.default_rng(2)
    seq = _rnd(10000, rng) + _array(_rnd(40000, rng), 200000, 0.002, rng) + _rnd(10000, rng)
    arr = kmerper.analyse_region("c", 0, seq)
    assert len(arr) == 1 and arr[0].period == 40000
    assert kmerper.analyse_region("r", 0, _rnd(50000, rng)) == []


def test_kmer_short_units_left_to_short_pass_and_n_runs():
    from repeatviz import kmerper
    rng = np.random.default_rng(3)
    seq = _array(_rnd(30, rng), 20000, 0.01, rng) + "N" * 5000 + _array(_rnd(500, rng), 20000, 0.01, rng)
    arr = kmerper.analyse_region("c", 0, seq)
    assert [a.period for a in arr] == [500] and arr[0].start >= 25000


def test_kmer_hor_reports_monomer():
    from repeatviz import kmerper
    rng = np.random.default_rng(4)
    mono = _rnd(171, rng)
    hor = "".join(_array(mono, 171, 0.15, rng) for _ in range(6))       # 6 diverged monomers
    arr = kmerper.analyse_region("c", 0, _array(hor, 100000, 0.005, rng))
    assert arr and arr[0].period == 1026


def test_arrays_tsv_v02_columns(tmp_path):
    c = tandem.Call("c", 0, 6000, 3000, 1, "ACGT", 2, method="kmer", tandem_frac=0.9, period2=171)
    tandem.write_arrays([c], tmp_path / "a.tsv")
    df = pd.read_csv(tmp_path / "a.tsv", sep="\t")
    assert list(df.columns) == tandem.ARRAY_COLS and df["method"][0] == "kmer" and df["period2"][0] == 171
    assert df.columns[6] == "class"


def test_short_prefilter_keeps_arrays_drops_random():
    from repeatviz import kmerper
    rng = np.random.default_rng(5)
    parts = [_rnd(30000, rng), _array("AC", 200, 0.02, rng), _rnd(30000, rng),
             _array(_rnd(60, rng), 2000, 0.05, rng), _rnd(30000, rng),
             _array("TTAGGG", 3000, 0.01, rng), _rnd(10000, rng)]
    seq = "".join(parts)
    iv = kmerper.short_candidates(seq, chunk=25000)   # small chunk: arrays cross chunk cuts
    pos, arrs = 0, []
    for i, p in enumerate(parts):
        if i % 2:
            arrs.append((pos, pos + len(p)))
        pos += len(p)
    for a, b in arrs:
        cov = sum(max(0, min(e, b) - max(s, a)) for s, e in iv)
        assert cov >= 0.95 * (b - a), (a, b, iv)
    assert sum(e - s for s, e in iv) < 0.2 * len(seq)


def test_short_prefilter_parallel_matches_serial():
    from repeatviz import kmerper
    rng = np.random.default_rng(6)
    seq = _rnd(50000, rng) + _array("CAG", 900, 0.02, rng) + _rnd(50000, rng)
    one = kmerper.run_short_candidates([("x", seq)], threads=1, piece=40000)
    two = kmerper.run_short_candidates([("x", seq)], threads=2, piece=40000)
    assert one == two and one["x"] == kmerper.short_candidates(seq)


def test_drop_covered_keeps_calls_running_beyond_kmer_arrays():
    calls = [tandem.Call("c", 0, 1000, 66, 10, "A"), tandem.Call("c", 2000, 2100, 30, 5, "A")]
    spans = {"c": [(100, 500), (600, 700)]}
    assert [c.start for c in tandem.drop_covered(calls, spans)] == [0, 2000]
    assert tandem.drop_covered(calls, {"c": [(0, 950), (2000, 2100)]}) == []


def test_merge_runs():
    from repeatviz import kmerper
    iv = [(100, 130), (140, 200), (1000, 1010), (5000, 5100)]
    assert kmerper.merge_runs(iv, 6000) == [(0, 300), (4900, 5200)]


def test_tantan_prefilter_if_installed():
    import shutil
    import pytest
    from repeatviz import kmerper
    if not shutil.which("tantan"):
        pytest.skip("tantan not installed")
    rng = np.random.default_rng(7)
    seq = _rnd(20000, rng) + _array(_rnd(45, rng), 3000, 0.12, rng) + _rnd(20000, rng)   # diverged
    iv = kmerper.run_tantan_candidates([("x", seq)])["x"]
    assert sum(max(0, min(e, 23000) - max(s, 20000)) for s, e in iv) > 0.9 * 3000


def test_telo_scan_variants_and_density():
    from repeatviz import telo, tandem
    rnd = "ACGTTGCA" * 40
    # 400 bp of TTAGGG with a TTGGGG variant every 5th unit (ULTRA may call this period 18 / variant consensus)
    tel = "".join("TTGGGG" if i % 5 == 4 else "TTAGGG" for i in range(67))
    seq = rnd + tel + rnd + "CCCTAA" * 30 + rnd + "TTAGGG" * 3 + rnd
    runs = telo.scan_seq("c1", seq)
    assert [(r[4], r[2] - r[1] >= 390) for r in runs][:1] == [("+", True)]
    assert any(r[4] == "-" and r[2] - r[1] == 180 for r in runs)
    assert not any(r[2] - r[1] < 100 for r in runs)          # 3 copies = not a run
    # an ULTRA-style call with a period-18 consensus over the run gets relabelled
    a = tandem.Call("c1", len(rnd), len(rnd) + len(tel), 18, 100.0, "AGGGTTAGGGTTAGGGTT", 22.0)
    assert a.cls != "telomere"
    by = {"c1": [(s, e, n, st) for _, s, e, n, st in runs]}
    arrays, relab, added = telo.apply([a], by, make_call=lambda c, s, e, n, st: tandem.Call(
        c, s, e, 6, 6.0 * n, telo.MOTIFS[st], float(n), method="telo_scan", cls_override="telomere"))
    assert a.cls == "telomere" and relab == 1 and added == 1   # the CCCTAA run is added
