"""matplotlib rendering: one PDF per repeat map (all / short-unit tandem / long-unit tandem).

Each chromosome = a contig bar coloured as a repeat heatmap (blue = repeat-poor, grey = genome
median, red = repeat-rich) with a read-depth wave of the same height underneath.
Portrait pages; contigs are separated by a white gap with thin black outlines.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from .layout import Layout, Part, display_name  # noqa: E402

# ---- repeat maps: key -> (title, legend description, score function, min. mid-point) ----
MAPS = {
    "all": ("All repeats",
            "share of window in high-copy k-mers (>=10x) or tandem arrays",
            lambda r: np.maximum(r["frac_ge10"].values, r["tandem_all"].values), 0.0),
    "tandem": ("Tandem repeats - short units (<100 bp, incl. telomere)",
               "share of window in tandem arrays with unit <100 bp",
               lambda r: (r["tandem_telomere"] + r["tandem_micro"] + r["tandem_mini"]).values, 0.02),
    "long": ("Tandem repeats - long units (>=100 bp, satellites)",
             "share of window in tandem arrays with unit >=100 bp",
             lambda r: (r["tandem_satellite"] + r["tandem_macro"]).values, 0.02),
    "dup": ("Duplications (2-9 copies, e.g. segmental dup.)",
            "share of window in k-mers occurring 2-9x in the assembly",
            lambda r: np.clip(r["frac_ge2"].values - r["frac_ge10"].values, 0, 1), 0.0),
}
MAP_ORDER = ("all", "tandem", "long", "dup")

REPEAT_CMAP = LinearSegmentedColormap.from_list(
    "repeat", ["#2166ac", "#92c5de", "#efefef", "#ef8a62", "#67001f"])
SMOOTH = 5  # windows in the display moving average (100 kb windows -> 500 kb)
TELO_COLOR = "#ffd400"   # yellow (black outline): distinct from the red/blue heatmap, also for colour-blind readers
HIGH_COLOR = "#c0392b"
DEPTH_LINE = "#333333"
DEPTH_COLLAPSE = "#d6604d"
DEPTH_LOW = "#4393c3"
WARN_COLOR = "#ffd400"
DEPTH_TOP = 3.0  # wave top = 3x median

# sizes (points / axis units)
FS_TITLE, FS_SUB, FS_CHR, FS_CTG, FS_TICK, FS_LEG = 22, 14, 15, 10, 13, 12
BAR_H, DEP_H, GAP, LAB_H, ROW_SPACE = 0.30, 0.30, 0.04, 0.16, 0.22
ROW_H = LAB_H + BAR_H + GAP + DEP_H + ROW_SPACE


def _part_rows(win_by: dict, p: Part) -> pd.DataFrame | None:
    g = win_by.get(p.contig)
    if g is None:
        return None
    g = g[(g["end"] > p.cstart) & (g["start"] < p.cend)]
    return g.iloc[::-1] if p.orient == "-" else g


def _draw_strip(ax, x0, x1, y, h, vals):
    """vals: 0-1 positions on REPEAT_CMAP. Drawn as vector quads (pcolormesh): an embedded
    1-row image gets resampled by the PDF backend/viewers into diagonal streaks."""
    if len(vals) == 0:
        return
    xe = np.linspace(x0, x1, len(vals) + 1)
    ax.pcolormesh(xe, [y, y + h], np.asarray(vals, dtype=float)[None, :], cmap=REPEAT_CMAP,
                  vmin=0, vmax=1, shading="flat", linewidth=0, antialiased=False, zorder=2)


def _norm(score: np.ndarray, mid: float) -> np.ndarray:
    """0 -> 0 (blue), mid -> 0.5 (grey), 1 -> 1 (dark red)."""
    x = np.where(score <= mid, 0.5 * score / max(mid, 1e-6),
                 0.5 + 0.5 * (score - mid) / max(1 - mid, 1e-6))
    return np.clip(x, 0, 1)


def _flag(ax, x, y, h, w, left=False):
    """telomere flag: pole + red pennant pointing away from the chromosome body."""
    ax.plot([x, x], [y, y + h], color="black", lw=1.0, zorder=7, clip_on=False, solid_capstyle="butt")
    d = -w if left else w
    ax.add_patch(plt.Polygon([(x, y + h), (x + d, y + h - h * 0.27), (x, y + h - h * 0.54)],
                             closed=True, fc=TELO_COLOR, ec="black", lw=0.5, zorder=8, clip_on=False))


def _runs(mask: np.ndarray, max_gap: int = 2) -> list[tuple[int, int]]:
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        return []
    out, s, e = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - e <= max_gap + 1:
            e = i
        else:
            out.append((s, e))
            s = e = i
    out.append((s, e))
    return out


def draw_object(ax, y, name, parts, length, win_by, arrays_by, alerts_by, has_depth, collapse,
                warn, score_fn, mid, xspan, status="", show_alerts=False, label_contigs=True,
                depth_label=True):
    """y = top of the row. Returns nothing; row height is ROW_H."""
    yb = y - LAB_H - BAR_H                      # bar bottom
    gap_w = xspan * 0.0018                      # white gap between contigs
    flag_w = xspan * 0.012
    peaks = []                                  # (x, depth) of high-depth runs, merged per row below
    for i, p in enumerate(parts):
        x0, x1 = p.pos, p.pos + p.length
        lx0 = x0 + (gap_w / 2 if i else 0)
        lx1 = x1 - (gap_w / 2 if i < len(parts) - 1 else 0)
        rows = _part_rows(win_by, p)
        if rows is not None and not rows.empty:
            sc = pd.Series(score_fn(rows)).rolling(SMOOTH, center=True, min_periods=1).mean().values
            _draw_strip(ax, lx0, lx1, yb, BAR_H, _norm(sc, mid))
        else:
            ax.add_patch(Rectangle((lx0, yb), lx1 - lx0, BAR_H, fc="#dddddd", lw=0, zorder=2))
        ax.add_patch(Rectangle((lx0, yb), lx1 - lx0, BAR_H, fc="none", ec="black", lw=0.7, zorder=4))
        if label_contigs:
            ya = yb + BAR_H + LAB_H * 0.45
            ax.annotate("", xy=(lx0, ya), xytext=(lx1, ya), zorder=4,
                        arrowprops=dict(arrowstyle="<|-|>", color="#555555", lw=0.8,
                                        shrinkA=0, shrinkB=0, mutation_scale=7))
            if p.length > xspan * 0.05:
                ax.text((lx0 + lx1) / 2, ya,
                        f"{p.contig.replace('ptg0000', '').replace('ptg000', '')}({p.orient})",
                        ha="center", va="center", fontsize=FS_CTG, color="#333333", zorder=5,
                        bbox=dict(fc="white", ec="none", pad=1.2))
        # depth wave (same height as the bar)
        if has_depth and rows is not None and not rows.empty:
            y0 = yb - GAP - DEP_H
            rd = rows["rel_depth"].values.astype(float)
            rdc = np.clip(rd, 0, DEPTH_TOP)
            xs = np.linspace(lx0, lx1, len(rd) + 1)
            xs = (xs[:-1] + xs[1:]) / 2
            yy = y0 + rdc / DEPTH_TOP * DEP_H
            base = y0 + DEP_H / DEPTH_TOP
            ax.add_patch(Rectangle((lx0, y0), lx1 - lx0, DEP_H, fc="#fafafa", ec="#bbbbbb", lw=0.4, zorder=1))
            ax.plot([lx0, lx1], [base, base], color="#999999", lw=0.5, ls=(0, (2, 2)), zorder=1)
            ax.fill_between(xs, base, yy, where=rdc >= collapse, color=DEPTH_COLLAPSE, lw=0, zorder=2,
                            interpolate=True)
            ax.fill_between(xs, base, yy, where=rdc <= 0.5, color=DEPTH_LOW, lw=0, zorder=2, interpolate=True)
            ax.plot(xs, yy, color=DEPTH_LINE, lw=0.6, zorder=3)
            absd = rows["depth"].values.astype(float) if "depth" in rows else rd
            for s, e in _runs(np.nan_to_num(rd) >= warn):
                k = s + int(np.nanargmax(absd[s:e + 1]))
                peaks.append((xs[k], absd[k]))
        # telomere flags
        a = arrays_by.get(p.contig)
        if a is not None:
            t = a[(a["class"] == "telomere") & (a["length"] >= 300)
                  & (a["end"] > p.cstart) & (a["start"] < p.cend)]
            for r in t.itertuples(index=False):
                xx = p.to_chrom(min(max((r.start + r.end) // 2, p.cstart), p.cend))
                _flag(ax, xx, yb + BAR_H * 0.5, BAR_H * 0.5 + LAB_H + 0.08, flag_w, left=xx < length / 2)
        al = alerts_by.get(p.contig) if show_alerts else None
        if al is not None:
            al = al[(al["end"] > p.cstart) & (al["start"] < p.cend)]
            for r in al.itertuples(index=False):
                xx = p.to_chrom(min(max((r.start + r.end) // 2, p.cstart), p.cend))
                ax.plot(xx, yb + BAR_H + 0.015, marker="v", ms=4, color="#ff7f00", mec="none", zorder=6,
                        clip_on=False)
    # high-depth marks: merge peaks closer than 4 % of the page width (e.g. both sides of a junction)
    if has_depth and peaks:
        peaks.sort()
        merged = [list(peaks[0])]
        for x, d in peaks[1:]:
            if x - merged[-1][0] < xspan * 0.04:
                if d > merged[-1][1]:
                    merged[-1] = [x, d]
            else:
                merged.append([x, d])
        y0 = yb - GAP - DEP_H
        for x, d in merged:
            ax.plot(x, y0 + DEP_H, marker="v", ms=8, color=HIGH_COLOR, mec="white", mew=0.5, zorder=10,
                    clip_on=False)
            ha = "left" if x < xspan * 0.06 else ("right" if x > length - xspan * 0.06 else "center")
            dx = xspan * 0.008 * (1 if ha == "left" else (-1 if ha == "right" else 0))
            ax.text(x + dx, y0 + DEP_H * 0.72, f"high: max {d:.0f}x", ha=ha, va="top",
                    fontsize=FS_CTG - 2, color=HIGH_COLOR, zorder=10, clip_on=False,
                    bbox=dict(fc="white", ec="none", alpha=1.0, pad=0.6))
    ax.text(-xspan * 0.015, yb + BAR_H / 2, display_name(name), ha="right", va="center",
            fontsize=FS_CHR if label_contigs else FS_CHR - 3, color="#c0392b" if status == "draft" else "black")
    if has_depth and depth_label:
        ax.text(-xspan * 0.015, yb - GAP - DEP_H / 2, "depth", ha="right", va="center",
                fontsize=FS_CTG, color="dimgrey")


def _legend(fig, map_key, mid, has_depth, collapse, warn, show_alerts):
    title, desc, _, _ = MAPS[map_key]
    # big colour bar
    cax = fig.add_axes([0.20, 0.088, 0.60, 0.016])
    grad = np.linspace(0, 1, 512)[None, :]
    cax.imshow(REPEAT_CMAP(grad), aspect="auto", extent=(0, 1, 0, 1))
    cax.set_yticks([])
    cax.set_xticks([0, 0.5, 1])
    cax.set_xticklabels([f"repeat-poor\n(0%)", f"genome median\n({mid:.0%})", f"repeat-rich\n(100%)"],
                        fontsize=FS_LEG)
    for s in cax.spines.values():
        s.set_linewidth(0.6)
    cax.set_title(f"repeat heatmap: {desc}", fontsize=FS_LEG, pad=6)
    h = [plt.Line2D([], [], marker=">", ms=12, ls="", color=TELO_COLOR, mec="black",
                    label="telomere (>=300 bp array)"),
         plt.Line2D([], [], color="#555555", lw=1, marker="<", ms=5, label="contig (name between arrows)")]
    if has_depth:
        h += [plt.Line2D([], [], color=DEPTH_LINE, lw=1.2, label=f"read depth (dotted = median, top {DEPTH_TOP:g}x)"),
              Patch(fc=DEPTH_COLLAPSE, label=f"depth >= {collapse:g}x median"),
              Patch(fc=DEPTH_LOW, label="depth <= 0.5x median"),
              plt.Line2D([], [], ls="", marker="v", ms=10, color=HIGH_COLOR,
                         label=f"very high depth (>= {warn:g}x median), label = peak")]
    if show_alerts:
        h += [plt.Line2D([], [], marker="v", ls="", ms=7, color="#ff7f00", label="alert: high-copy, little tandem")]
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=FS_LEG,
               frameon=False, handlelength=1.6, columnspacing=1.5)


def render(layout: Layout, win: pd.DataFrame, arrays: pd.DataFrame, alerts: pd.DataFrame,
           out_prefix: Path, title: str, has_depth: bool, high: float,
           per_page: int = 15, width_in: float = 11, height_in: float = 17, max_unplaced: int = 60,
           unplaced_per_page: int = 30, show_alerts: bool = False, maps=MAP_ORDER,
           depth_warn: float = 2.0) -> list[Path]:
    win_by = {c: g for c, g in win.groupby("contig", sort=False)}
    arrays_by = {c: g for c, g in arrays.groupby("contig", sort=False)}
    alerts_by = {c: g for c, g in alerts.groupby("contig", sort=False)}
    objects = [(c.name, c.parts, c.length, c.status) for c in layout.chroms]
    unp = sorted(layout.unplaced, key=lambda x: -x[1])
    unp_objs = [(f"{c}", [Part(c, 0, L, "+", 0)], L, "") for c, L in unp[:max_unplaced]]
    pages = [("chromosomes", objects[i:i + per_page]) for i in range(0, len(objects), per_page)]
    pages += [("UNPLACED", unp_objs[i:i + unplaced_per_page]) for i in range(0, len(unp_objs), unplaced_per_page)]
    written = []
    # v0.3: all PDFs together in <prefix dir>/repeatviz_pdf/ (PNGs and tables stay next to the prefix)
    pdf_dir = Path(out_prefix).parent / "repeatviz_pdf"
    pdf_dir.mkdir(parents=True, exist_ok=True)
    for mk in maps:
        mtitle, _, score_fn, mid_floor = MAPS[mk]
        mid = max(float(np.nanmedian(score_fn(win))) if len(win) else 0.3, mid_floor)
        pdf_path = pdf_dir / f"{Path(out_prefix).name}.{mk}.pdf"
        with PdfPages(pdf_path) as pdf:
            for pi, (kind, objs) in enumerate(pages, 1):
                is_unp = kind == "UNPLACED"
                n_rows = per_page if not is_unp else unplaced_per_page
                maxlen = max(o[2] for o in objs)
                fig = plt.figure(figsize=(width_in, height_in))
                left = 0.11 if not is_unp else 0.20
                ax = fig.add_axes([left, 0.155, 0.97 - left, 0.75])
                y = 0.0
                for name, parts, length, status in objs:
                    draw_object(ax, y, name, parts, length, win_by, arrays_by, alerts_by, has_depth,
                                high, depth_warn, score_fn, mid, maxlen, status, show_alerts,
                                label_contigs=not is_unp, depth_label=not is_unp)
                    y -= ROW_H
                ax.set_xlim(0, maxlen * 1.03)
                ax.set_ylim(-ROW_H * n_rows, 0.05)
                ax.set_yticks([])
                lg = np.log10(max(maxlen, 10))
                step = 10 ** int(lg) / (2 if maxlen / 10 ** int(lg) < 3 else 1)
                ticks = np.arange(0, maxlen + step, step)
                ax.set_xticks(ticks)
                ax.set_xticklabels([f"{t / 1e6:g}" for t in ticks], fontsize=FS_TICK)
                ax.set_xlabel("Mb", fontsize=FS_TICK)
                for s in ("left", "right", "top"):
                    ax.spines[s].set_visible(False)
                fig.text(0.03, 0.975, mtitle, fontsize=FS_TITLE, fontweight="bold", va="top")
                sub = (f"UNPLACED (longest {len(unp_objs)} of {len(unp)}, "
                       f"{sum(L for _, L in unp) / 1e6:.1f} Mb total)" if is_unp else "chromosomes")
                fig.text(0.03, 0.948, f"{title}  |  {sub}  |  page {pi}/{len(pages)}",
                         fontsize=FS_SUB, va="top", color="#444444")
                _legend(fig, mk, mid, has_depth, high, depth_warn, show_alerts)
                pdf.savefig(fig)
                png = Path(f"{out_prefix}.{mk}.p{pi}.png")
                fig.savefig(png, dpi=110)
                written.append(png)
                plt.close(fig)
        written.insert(0, pdf_path)
    return written
