"""Layout-aware summaries: depth track, contig-end report, alerts, arrays on chromosomes."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from .layout import Layout, contig_index
from .tandem import CLASSES


# ---------------------------------------------------------------- depth
def read_depth(path: Path) -> pd.DataFrame:
    """mosdepth *.regions.bed.gz or a bedGraph: contig start end depth."""
    op = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.split()
            rows.append((f[0], int(f[1]), int(f[2]), float(f[-1])))
    return pd.DataFrame(rows, columns=["contig", "start", "end", "depth"])


def add_depth(win: pd.DataFrame, depth: pd.DataFrame, high: float, low: float) -> tuple[pd.DataFrame, float]:
    """Length-weighted mean depth per window; rel_depth = depth / genome median."""
    d = depth[depth["depth"] > 0]
    med = float(np.median(np.repeat(d["depth"].values, ((d["end"] - d["start"]) // 1000).clip(1))))
    out = []
    for c, g in win.groupby("contig", sort=False):
        dd = depth[depth["contig"] == c]
        if dd.empty:
            out.append(np.full(len(g), np.nan))
            continue
        ds, de, dv = dd["start"].values, dd["end"].values, dd["depth"].values
        vals = []
        for s, e in zip(g["start"].values, g["end"].values):
            lo = np.searchsorted(de, s, side="right")
            hi = np.searchsorted(ds, e, side="left")
            ov = np.minimum(de[lo:hi], e) - np.maximum(ds[lo:hi], s)
            vals.append(float((dv[lo:hi] * ov).sum() / ov.sum()) if ov.sum() > 0 else np.nan)
        out.append(np.array(vals))
    win = win.copy()
    win["depth"] = np.concatenate(out) if out else []
    win["rel_depth"] = (win["depth"] / med).round(3)
    win["depth_flag"] = np.where(win["rel_depth"] >= high, "collapse",
                                 np.where(win["rel_depth"] < low, "low", ""))
    return win, med


# ---------------------------------------------------------------- alerts
def alerts(win: pd.DataFrame, hc: float, max_tandem: float) -> pd.DataFrame:
    """Windows that are mostly high-copy k-mers but carry little tandem annotation.
    With the k-mer long pass (no period limit) these are mostly segmental-duplication /
    interspersed-repeat clusters or diverged arrays the pass could not resolve; with
    --long-mode ultra/none they also include tandem units longer than --p-long.
    Adjacent windows are merged."""
    m = win[(win["frac_ge10"] >= hc) & (win["tandem_all"] < max_tandem)]
    rows = []
    for c, g in m.groupby("contig", sort=False):
        cur = None
        for r in g.itertuples(index=False):
            if cur and r.start == cur["end"]:
                cur["end"] = r.end
                cur["n"] += 1
                cur["hc"] += r.frac_ge10
                cur["td"] += r.tandem_all
            else:
                if cur:
                    rows.append(cur)
                cur = {"contig": c, "start": r.start, "end": r.end, "n": 1,
                       "hc": r.frac_ge10, "td": r.tandem_all}
        if cur:
            rows.append(cur)
    df = pd.DataFrame(rows, columns=["contig", "start", "end", "n", "hc", "td"])
    if df.empty:
        return pd.DataFrame(columns=["contig", "start", "end", "length", "mean_frac_high_copy",
                                     "mean_frac_tandem", "alert"])
    df["length"] = df["end"] - df["start"]
    df["mean_frac_high_copy"] = (df["hc"] / df["n"]).round(3)
    df["mean_frac_tandem"] = (df["td"] / df["n"]).round(3)
    df["alert"] = "high-copy, little tandem: SD/interspersed cluster or unresolved array"
    return df[["contig", "start", "end", "length", "mean_frac_high_copy", "mean_frac_tandem", "alert"]]


# ---------------------------------------------------------------- projection
def project(df: pd.DataFrame, layout: Layout) -> pd.DataFrame:
    """Add chrom / chrom_start / chrom_end for rows (contig, start, end) placed in the layout."""
    idx = contig_index(layout)
    ch, cs, ce = [], [], []
    for r in df.itertuples(index=False):
        hit = None
        for c, p in idx.get(r.contig, ()):
            if r.end > p.cstart and r.start < p.cend:
                a = p.to_chrom(max(r.start, p.cstart))
                b = p.to_chrom(min(r.end, p.cend))
                hit = (c.name, min(a, b), max(a, b))
                break
        ch.append(hit[0] if hit else "")
        cs.append(hit[1] if hit else -1)
        ce.append(hit[2] if hit else -1)
    df = df.copy()
    df["chrom"], df["chrom_start"], df["chrom_end"] = ch, cs, ce
    return df


# ---------------------------------------------------------------- contig ends
def contig_ends(win: pd.DataFrame, arrays: pd.DataFrame, layout: Layout, lengths: dict[str, int],
                flank: int, hc: float, lc: float, tandem_min: float, telo_min_bp: int,
                telo_end: int) -> pd.DataFrame:
    """One row per contig end: what repeats sit within `flank` bp of it.

    role: chr_end (outer end of a chromosome), junction (faces another contig),
    unplaced. category joins every signal that applies; 'no repeat' when none
    does (break likely not repeat-driven: haplotype divergence, coverage, ...).
    """
    idx = contig_index(layout)
    role = {}
    for ch in layout.chroms:
        for i, p in enumerate(ch.parts):
            first_side = "start" if p.orient == "+" else "end"   # side facing chromosome start
            last_side = "end" if p.orient == "+" else "start"
            role[(p.contig, first_side)] = (ch.name, "chr_end" if i == 0 else "junction")
            role[(p.contig, last_side)] = (ch.name, "chr_end" if i == len(ch.parts) - 1 else "junction")
    has_depth = "rel_depth" in win.columns
    arr_by = {c: g for c, g in arrays.groupby("contig")}
    win_by = {c: g for c, g in win.groupby("contig")}
    rows = []
    for c, L in lengths.items():
        a = arr_by.get(c, arrays.iloc[0:0])
        w = win_by.get(c)
        for side in ("start", "end"):
            s, e = (0, min(flank, L)) if side == "start" else (max(0, L - flank), L)
            ov = np.minimum(a["end"], e) - np.maximum(a["start"], s)
            sel = a[ov > 0].assign(ov=ov[ov > 0])
            bp = {k: int(sel.loc[sel["class"] == k, "ov"].sum()) for k in CLASSES}
            top = sel.sort_values("ov", ascending=False).head(1)
            tes, tee = (0, min(telo_end, L)) if side == "start" else (max(0, L - telo_end), L)
            tov = np.minimum(a["end"], tee) - np.maximum(a["start"], tes)
            telo_tip = int(tov[(tov > 0) & (a["class"] == "telomere")].sum())
            ww = w[(w["end"] > s) & (w["start"] < e)] if w is not None else None
            f10 = float(ww["frac_ge10"].mean()) if ww is not None and len(ww) else 0.0
            f2 = float(ww["frac_ge2"].mean()) if ww is not None and len(ww) else 0.0
            dflag = ""
            rd = np.nan
            if has_depth and ww is not None and len(ww):
                rd = float(ww["rel_depth"].median())
                fl = set(ww["depth_flag"].dropna()) - {""}
                dflag = ",".join(sorted(fl))
            span = e - s
            cats = []
            if telo_tip >= telo_min_bp:
                cats.append("telomere")
            if (bp["satellite"] + bp["macro"]) / span >= tandem_min:
                cats.append("satellite")
            if (bp["micro"] + bp["mini"]) / span >= tandem_min:
                cats.append("short-tandem")
            if f10 >= hc:
                cats.append("high-copy")
            elif f2 >= lc:
                cats.append("low-copy/SD")
            if "collapse" in dflag:
                cats.append("collapse")
            chrom, r = role.get((c, side), ("", "unplaced" if c not in idx else "internal"))
            rows.append({
                "contig": c, "end": side, "contig_len": L, "chrom": chrom, "role": r,
                "flank_bp": span,
                **{f"bp_{k}": bp[k] for k in CLASSES},
                "telomere_at_tip_bp": telo_tip,
                "top_period": int(top["period"].iloc[0]) if len(top) else "",
                "top_class": top["class"].iloc[0] if len(top) else "",
                "top_unit": (top["consensus"].iloc[0][:60] if len(top) else ""),
                "frac_high_copy": round(f10, 3), "frac_low_copy": round(f2, 3),
                "rel_depth_median": round(rd, 2) if rd == rd else "",
                "depth_flags": dflag,
                "category": "+".join(cats) if cats else "no repeat",
            })
    return pd.DataFrame(rows)
