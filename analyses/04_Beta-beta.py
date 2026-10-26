#!/usr/bin/env python3
"""
Beta-beta comparison of clinical PD GWAS signals and pathological LBD GWAS effects
- Project: Genetic determinants of Lewy body pathology
- Version: Python 3; numpy, pandas, scipy, matplotlib, openpyxl, adjustText (optional)
- Last Updated: 16-SEP-2026

Overview
- For every clinical PD independent signal found in the LBD summary statistics,
  alleles are harmonised and both betas are oriented to the PD risk-increasing allele.
- Reports sign concordance (binomial test), OLS slope, and York (errors-in-both-variables,
  Deming-type) slope overall and excluding high-leverage loci; draws the beta-beta plot.

Usage
    python 04_Beta-beta.py --loci PD_independent_signals.xlsx --lbd LBD_meta.filtered.tab \\
        --table-out loci_matched.tsv --above-out loci_above_line.tsv --out beta_beta.png
    python 04_Beta-beta.py --from-table loci_matched.tsv --view zoom --genes SNCA GBA1 TMEM175 HLA
  --loci  .xlsx with CHR, BP, Nearest Gene, Known or Novel, Reference allele,
          Alternative allele, Effect allele, and Beta/SE/P columns for --pd-set
  --lbd   tab-separated with CHR, BP, Allele1 (effect allele), Allele2, Effect, StdErr, P-value
  Both inputs must be on the same genome build. See --help for all options.

"""
import argparse

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# --pd-set → (beta, SE, P) column names in the PD loci table
PD_SETS = {
    "case-control": ("Beta, case-control", "SE, case-control", "P, case-control"),
    "all":          ("Beta, all studies", "SE, all studies", "P, all studies"),
    "biobank":      ("Beta, biobank", "SE, biobank", "P, biobank"),
}
DEFAULT_LEVERAGE_GENES = ["GBA1", "SNCA", "TMEM175"]
PLOT_COLS = ["CHR", "BP", "GENE", "STATUS", "EA", "OA", "BETA_pd", "SE_pd", "P_pd",
             "BETA_lbd_al", "SE_lbd", "P_lbd"]


def york_regression(x, y, sx, sy, tol=1e-12, max_iter=200):
    """York (2004) errors-in-both-variables regression with per-point standard
    errors (sx, sy) and uncorrelated errors. Returns slope, intercept, slope SE."""
    x, y, sx, sy = map(np.asarray, (x, y, sx, sy))
    wx, wy = 1.0 / sx**2, 1.0 / sy**2
    b = np.polyfit(x, y, 1)[0]                     # OLS starting slope
    for _ in range(max_iter):
        W = (wx * wy) / (wx + b * b * wy)
        Wsum = W.sum()
        x_bar = (W * x).sum() / Wsum
        y_bar = (W * y).sum() / Wsum
        U, V = x - x_bar, y - y_bar
        beta = W * (U / wy + b * V / wx)
        b_new = (W * beta * V).sum() / (W * beta * U).sum()
        if abs(b_new - b) < tol * (abs(b_new) + 1e-30):
            b = b_new
            break
        b = b_new
    # final quantities at the converged slope
    W = (wx * wy) / (wx + b * b * wy)
    Wsum = W.sum()
    x_bar = (W * x).sum() / Wsum
    y_bar = (W * y).sum() / Wsum
    U, V = x - x_bar, y - y_bar
    beta = W * (U / wy + b * V / wx)
    a = y_bar - b * x_bar
    x_adj = x_bar + beta
    x_adj_bar = (W * x_adj).sum() / Wsum
    u = x_adj - x_adj_bar
    slope_se = np.sqrt(1.0 / (W * u * u).sum())
    return b, a, slope_se


def load_pd_loci(path, pd_set):
    bcol, scol, pcol = PD_SETS[pd_set]
    df = pd.read_excel(path)
    df = df.rename(columns={
        "Nearest Gene": "GENE", "Known or Novel": "STATUS",
        "Reference allele": "REF", "Alternative allele": "ALT",
        "Effect allele": "EA",
        bcol: "BETA_pd", scol: "SE_pd", pcol: "P_pd"})
    df["CHR"] = df["CHR"].astype(str).str.replace("chr", "", regex=False)
    df["BP"] = df["BP"].astype(int)
    for c in ["EA", "REF", "ALT"]:
        df[c] = df[c].astype(str).str.upper()
    df["OA"] = np.where(df["EA"] == df["REF"], df["ALT"], df["REF"])
    return df[["CHR", "BP", "GENE", "STATUS", "EA", "OA", "BETA_pd", "SE_pd", "P_pd"]]


def load_lbd(path):
    cols = ["CHR", "BP", "Allele1", "Allele2", "Effect", "StdErr", "P-value"]
    df = pd.read_csv(path, sep="\t", usecols=cols, dtype={"CHR": str})
    df = df.rename(columns={"Allele1": "A1", "Allele2": "A2",
                            "Effect": "BETA_lbd", "StdErr": "SE_lbd", "P-value": "P_lbd"})
    df["CHR"] = df["CHR"].str.replace("chr", "", regex=False)
    df["A1"] = df["A1"].str.upper()
    df["A2"] = df["A2"].str.upper()
    return df


def harmonise(pdf, lbd):
    """Merge PD loci onto the LBD GWAS, align alleles, and orient every variant to
    the PD risk-increasing allele. Returns the plotted table."""
    print(f"PD independent signals: {len(pdf)}")
    m = pdf.merge(lbd, on=["CHR", "BP"], how="inner")
    print(f"Matched in LBD by CHR:BP: {len(m)} / {len(pdf)}")

    same = (m.EA == m.A1) & (m.OA == m.A2)
    swap = (m.EA == m.A2) & (m.OA == m.A1)
    m = m[same | swap].copy()
    is_same = (m.EA == m.A1) & (m.OA == m.A2)
    m["BETA_lbd_al"] = np.where(is_same, m.BETA_lbd, -m.BETA_lbd)
    print(f"Allele-harmonised loci plotted: {len(m)}")
    if len(m) < 2:
        raise SystemExit("Too few harmonised loci; check builds/alleles.")

    # orient every variant to the PD risk-increasing allele (PD beta positive)
    flip = m.BETA_pd < 0
    m.loc[flip, ["BETA_pd", "BETA_lbd_al"]] *= -1
    m.loc[flip, ["EA", "OA"]] = m.loc[flip, ["OA", "EA"]].values
    print(f"Flipped to PD risk allele: {int(flip.sum())} / {len(m)}")
    return m


def load_plot_table(path):
    """Read a harmonised table previously written by --table-out."""
    m = pd.read_csv(path, sep="\t", dtype={"CHR": str})
    missing = [c for c in PLOT_COLS if c not in m.columns]
    if missing:
        raise SystemExit(f"{path} is missing columns: {', '.join(missing)}")
    print(f"Loaded {len(m)} harmonised loci from {path}")
    return m


def view_limits(x, y, args):
    """Axis limits for the requested view. Returns (xlo, xhi, ylo, yhi)."""
    if args.view == "zoom":
        hi = args.zoom_lim if args.zoom_lim is not None else \
            float(np.quantile(np.concatenate([x, y]), args.zoom_quantile))
        lo_y = min(0.0, y.min())
        pad = 0.10 * (hi - lo_y)
        return -pad, hi + pad, lo_y - pad, hi + pad
  
    padx = 0.07 * (x.max() - min(0.0, x.min()))
    pady = 0.07 * (y.max() - y.min())
    pad = max(padx, pady)
    # the rightmost point usually carries a label, so give that edge extra room
    return min(0.0, x.min()) - pad, x.max() + 1.8 * pad, y.min() - pad, y.max() + pad


def figure_size(xlo, xhi, ylo, yhi, args):
    """Size the canvas to the data box so an equal-aspect panel does not letterbox."""
    if args.figsize:
        return tuple(args.figsize)
    xr, yr = xhi - xlo, yhi - ylo
    base, margin = args.panel_size, 1.3
    if yr >= xr:
        h, w = base, base * xr / yr
    else:
        w, h = base, base * yr / xr
    return (np.clip(w, 3.2, 12) + margin, np.clip(h, 3.2, 12) + margin)


def gene_match(series, tokens):
    """True where a gene name matches a requested token: exact, or the token as the
    prefix of a hyphenated family name so 'HLA' picks up HLA-DRB1."""
    up = series.str.upper()
    hit = pd.Series(False, index=series.index)
    for t in tokens:
        t = t.strip().upper()
        hit |= (up == t) | up.str.startswith(t + "-")
    return hit


def draw_points(ax, m, s, edge_lw, alpha, legend=True):
    colours = {"Known": "steelblue", "Novel": "darkorange"}
    # groupby drops missing keys, which would silently remove those loci from the plot
    for status, sub in m.assign(STATUS=m.STATUS.fillna("Unknown")).groupby("STATUS"):
        ax.scatter(sub.BETA_pd, sub.BETA_lbd_al, s=s,
                   facecolor=colours.get(status, "grey"), edgecolor="black",
                   linewidth=edge_lw, alpha=alpha,
                   label=status if legend else None, zorder=3)
    return colours


def draw_offscale(ax, off, named, colours, xlo, xhi, ylo, yhi):
    """Draw points outside the view as outward triangles on the frame, so a zoomed
    panel still shows that they exist. `named` is the set of row labels to annotate.
    Returns (x, y, text, align) label items, placed later with the in-view labels."""
    xr, yr = xhi - xlo, yhi - ylo
    items = []
    for idx, r in off.iterrows():
        right, left = r.BETA_pd > xhi, r.BETA_pd < xlo
        top, bottom = r.BETA_lbd_al > yhi, r.BETA_lbd_al < ylo
        mk = "^" if top else "v" if bottom else ">" if right else "<"
        # nudge the marker just inside the frame so it is not clipped in half
        cx = np.clip(r.BETA_pd, xlo + 0.015 * xr, xhi - 0.015 * xr)
        cy = np.clip(r.BETA_lbd_al, ylo + 0.015 * yr, yhi - 0.015 * yr)
        ax.scatter([cx], [cy], marker=mk, s=70, zorder=6,
                   facecolor=colours.get(r.STATUS, "grey"), edgecolor="black", linewidth=0.5)
        if idx in named:
            # the text states the true coordinates, since the marker sits on the frame;
            # align it towards the inside of the panel so it is not cut off
            align = dict(ha="right" if right else "left" if left else "center",
                         va="top" if top else "bottom" if bottom else "center")
            items.append((cx, cy, f"{r.GENE}\n({r.BETA_pd:.2f}, {r.BETA_lbd_al:.2f})", align))
    return items


def pick_labels(m, genes, lbd_p, tag=""):
    """Rows to label: the best-P_lbd signal of each requested gene, or of every gene
    passing the LBD significance threshold when no gene list is given."""
    if genes:
        sel = m[gene_match(m.GENE, genes)]
        out = sel.loc[sel.groupby("GENE")["P_lbd"].idxmin()] if len(sel) else sel
        print(f"{tag}labelling requested genes: {', '.join(sorted(out.GENE))}")
        absent = [g for g in genes if not gene_match(m.GENE, [g]).any()]
        if absent:
            print(f"{tag}  not among the loci shown in this panel: {', '.join(absent)}")
    else:
        thresh = lbd_p if lbd_p is not None else 0.05 / len(m)
        sig = m[m.P_lbd < thresh]
        out = sig.loc[sig.groupby("GENE")["P_lbd"].idxmin()]
        print(f"{tag}LBD-significant loci labelled (P_lbd < {thresh:.2e}): {len(out)}")
    return out


def label_items(rows):
    """(x, y, gene, align) label items for rows drawn at their true coordinates."""
    return [(r.BETA_pd, r.BETA_lbd_al, r.GENE, {}) for _, r in rows.iterrows()]


def place_labels(ax, items, fontsize, lw=0.5):
    """Draw italic labels at their (x, y) anchors, then let adjustText spread them."""
    texts = [ax.text(x, y, text, fontsize=fontsize, fontstyle="italic",
                     fontweight="bold", zorder=7, **align)
             for x, y, text, align in items]
    if not texts:
        return
    try:
        from adjustText import adjust_text
    except ImportError:
        print("Tip: pip install adjustText for non-overlapping labels")
        return
    kw = dict(ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=lw))
    # min_arrow_len=0 keeps a leader line on every moved label, so a gene name is
    # never left floating next to the wrong point; fall back for older adjustText
    for extra in (dict(ensure_inside_axes=True, min_arrow_len=0),
                  dict(ensure_inside_axes=True), {}):
        try:
            adjust_text(texts, **kw, **extra)
            return
        except TypeError:
            continue


def parse_args():
    ap = argparse.ArgumentParser(
        description="Beta-beta plot of clinical PD GWAS independent signals against a "
                    "pathological LBD GWAS.",
        epilog="See the module docstring for input formats and examples.")
    src = ap.add_argument_group("input (either --loci and --lbd, or --from-table)")
    src.add_argument("--loci", help="PD independent-signal table (.xlsx)")
    src.add_argument("--lbd", help="LBD GWAS summary statistics (tab-separated)")
    src.add_argument("--from-table", help="replot from a TSV previously written by --table-out "
                                          "(skips re-reading the LBD summary statistics)")

    io = ap.add_argument_group("output")
    io.add_argument("--out", default="beta_beta_loci.png",
                    help="figure path; format from the extension (default beta_beta_loci.png)")
    io.add_argument("--table-out", default=None, help="TSV of all plotted loci")
    io.add_argument("--above-out", default=None, help="TSV of loci above the y = x line")

    ana = ap.add_argument_group("analysis")
    ana.add_argument("--pd-set", choices=list(PD_SETS), default="case-control",
                     help="which PD beta/SE/P columns to use (default case-control)")
    ana.add_argument("--leverage-genes", nargs="+", default=DEFAULT_LEVERAGE_GENES,
                     metavar="GENE",
                     help="genes excluded in the sensitivity York/Deming fit "
                          f"(default {' '.join(DEFAULT_LEVERAGE_GENES)})")
    ana.add_argument("--lbd-p", type=float, default=None,
                     help="label loci with P_lbd below this (default Bonferroni 0.05/n)")
    ana.add_argument("--genes", nargs="+", default=None, metavar="GENE",
                     help="label only these genes, ignoring --lbd-p. A bare family name "
                          "matches its members, e.g. HLA matches HLA-DRB1.")
    ana.add_argument("--no-fit-line", action="store_true",
                     help="do not draw the York/Deming fit line")
    ana.add_argument("--label1", default=None,
                     help="x-axis cohort label (default 'clinical PD GWAS (<pd-set>)')")
    ana.add_argument("--label2", default="pathological LBD GWAS", help="y-axis cohort label")

    view = ap.add_argument_group("view / zoom")
    view.add_argument("--view", choices=["full", "zoom", "inset"], default="full",
                      help="full: all data. zoom: magnify the core and mark off-scale loci "
                           "on the frame. inset: full view plus a magnified inset.")
    view.add_argument("--zoom-lim", type=float, default=None,
                      help="upper beta shown in the zoom view (default: --zoom-quantile of the data)")
    view.add_argument("--zoom-quantile", type=float, default=0.95,
                      help="quantile of the pooled PD and LBD betas that sets the zoom box "
                           "(default 0.95)")
    view.add_argument("--inset-lim", type=float, default=None,
                      help="upper beta shown in the inset (default: 0.90 quantile)")
    view.add_argument("--inset-genes", nargs="+", default=None, metavar="GENE",
                      help="genes to label inside the inset (default none; pass the same "
                           "list as --genes to mirror the main panel)")
    view.add_argument("--inset-label-size", type=float, default=None,
                      help="font size for inset gene labels (default --label-size minus 2)")
    view.add_argument("--inset-loc", default="upper left", choices=["lower right", "upper left"],
                      help="corner of the main panel holding the inset (default upper left)")
    view.add_argument("--inset-rect", type=float, nargs=4, default=None,
                      metavar=("X", "Y", "W", "H"),
                      help="exact inset box in axes fractions, e.g. --inset-rect 0.55 0.04 0.42 0.42")
    view.add_argument("--panel-size", type=float, default=7.0,
                      help="length of the longer panel edge in inches (default 7)")
    view.add_argument("--figsize", type=float, nargs=2, default=None, metavar=("W", "H"),
                      help="override the whole figure size in inches, e.g. --figsize 8 10")
    view.add_argument("--legend-loc", default="upper left",
                      help="matplotlib legend location (default 'upper left')")

    style = ap.add_argument_group("point / error-bar styling")
    style.add_argument("--point-size", type=float, default=40)
    style.add_argument("--point-alpha", type=float, default=0.9)
    style.add_argument("--label-size", type=float, default=9)
    style.add_argument("--errorbars", choices=["ci", "se", "none"], default="ci",
                       help="ci: 1.96*SE whiskers (default). se: 1*SE. none: hide them, "
                            "which uncrowds a dense core.")
    style.add_argument("--errorbar-alpha", type=float, default=1.0)
    return ap.parse_args()


def main():
    args = parse_args()

    if args.from_table:
        m = load_plot_table(args.from_table)
    elif args.loci and args.lbd:
        m = harmonise(load_pd_loci(args.loci, args.pd_set), load_lbd(args.lbd))
    else:
        raise SystemExit("Provide either --loci and --lbd, or --from-table.")

    label1 = args.label1 or f"clinical PD GWAS ({args.pd_set})"

    # ---- statistics (always on every harmonised locus) ----------------------
    x, y = m.BETA_pd.values, m.BETA_lbd_al.values
    ols = stats.linregress(x, y)
    conc = int(np.sum(np.sign(x) == np.sign(y)))
    binom = stats.binomtest(conc, len(m), 0.5).pvalue
    print(f"n={len(m)} | concordant {conc}/{len(m)} ({100*conc/len(m):.0f}%), "
          f"binomial P={binom:.2e}")
    print(f"OLS slope={ols.slope:.3f}, r={ols.rvalue:.3f}, P={ols.pvalue:.2e}")

    d_slope, d_int, d_se = york_regression(x, y, m.SE_pd.values, m.SE_lbd.values)
    lo, hi = d_slope - 1.96 * d_se, d_slope + 1.96 * d_se
    print(f"Deming slope={d_slope:.3f} (95% CI {lo:.3f}-{hi:.3f})")
    leverage = {g.upper() for g in args.leverage_genes}
    keep = ~m.GENE.str.upper().isin(leverage)
    ds2, di2, dse2 = york_regression(x[keep], y[keep], m.SE_pd.values[keep], m.SE_lbd.values[keep])
    print(f"Deming slope excl. {'/'.join(sorted(leverage))}: {ds2:.3f} "
          f"(95% CI {ds2-1.96*dse2:.3f}-{ds2+1.96*dse2:.3f})")

    above = m[m.BETA_lbd_al > m.BETA_pd].copy()
    above["delta_LBD_minus_PD"] = above.BETA_lbd_al - above.BETA_pd
    above = above.sort_values("delta_LBD_minus_PD", ascending=False)
    print(f"Loci above y=x (larger effect in pathological LBD): {len(above)}")

    # ---- figure ------------------------------------------------------------
    xlo, xhi, ylo, yhi = view_limits(x, y, args)
    fig, ax = plt.subplots(figsize=figure_size(xlo, xhi, ylo, yhi, args))
    # fix the view first: adjustText keeps labels inside the axes limits current when
    # it runs, so autoscaled limits here would push labels outside the final frame
    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi); ax.set_aspect("equal")
    diag =[min(xlo, ylo), max(xhi, yhi)]
    ax.plot(diag, diag, color="grey", ls=":", lw=0.9, zorder=1, label="y = x")
    ax.axhline(0, color="grey", lw=0.4); ax.axvline(0, color="grey", lw=0.4)

    if args.errorbars != "none":
        k = 1.96 if args.errorbars == "ci" else 1.0
        ax.errorbar(x, y, yerr=k * m.SE_lbd.values, xerr=k * m.SE_pd.values,
                    fmt="none", elinewidth=0.5, ecolor="lightgrey",
                    alpha=args.errorbar_alpha, zorder=2)
    colours = draw_points(ax, m, args.point_size, 0.4, args.point_alpha)

    if not args.no_fit_line:
        xs = np.array([xlo, xhi])
        ax.plot(xs, d_int + d_slope * xs, color="firebrick", lw=1.3, ls="--",
                zorder=3, label=f"Deming (slope {d_slope:.2f})")

    # the inset is drawn before the gene labels so adjustText can route them around it
    if args.view == "inset":
        ihi = args.inset_lim if args.inset_lim is not None else \
            float(np.quantile(np.concatenate([x, y]), 0.90))
        ilo_y = min(0.0, y.min()) - 0.05 * ihi
        rect = args.inset_rect or ([0.54, 0.05, 0.43, 0.43]
                                   if args.inset_loc == "lower right"
                                   else [0.09, 0.50, 0.43, 0.41])
        axin = ax.inset_axes(rect, zorder=8)
        axin.plot([ilo_y, ihi], [ilo_y, ihi], color="grey", ls=":", lw=0.8, zorder=1)
        axin.axhline(0, color="grey", lw=0.4); axin.axvline(0, color="grey", lw=0.4)
        draw_points(axin, m, args.point_size * 0.55, 0.3, args.point_alpha, legend=False)
        if not args.no_fit_line:
            xs = np.array([-0.05 * ihi, ihi])
            axin.plot(xs, d_int + d_slope * xs, color="firebrick", lw=1.0, ls="--", zorder=3)
        axin.set_xlim(-0.05 * ihi, ihi); axin.set_ylim(ilo_y, ihi)
        axin.set_aspect("equal")
        # inward, sparse ticks so they do not collide with the main axis labels
        axin.tick_params(labelsize=6, direction="in", pad=2, length=3)
        axin.xaxis.set_major_locator(MaxNLocator(4))
        axin.yaxis.set_major_locator(MaxNLocator(4))
        axin.patch.set_facecolor("white"); axin.patch.set_alpha(1.0)
        ax.indicate_inset_zoom(axin, edgecolor="black", lw=0.8, alpha=0.7)
        # the inset is opaque, so warn about any locus it hides in the main panel
        rx0, rx1 = xlo + rect[0] * (xhi - xlo), xlo + (rect[0] + rect[2]) * (xhi - xlo)
        ry0, ry1 = ylo + rect[1] * (yhi - ylo), ylo + (rect[1] + rect[3]) * (yhi - ylo)
        hidden = m[m.BETA_pd.between(rx0, rx1) & m.BETA_lbd_al.between(ry0, ry1) &
                   ((m.BETA_pd > ihi) | (m.BETA_lbd_al > ihi))]
        if len(hidden):
            print(f"WARNING: the inset covers {len(hidden)} locus/loci in the main panel: "
                  f"{', '.join(sorted(hidden.GENE.unique()))}. Move it with --inset-rect.")
        if args.inset_genes:
            # only loci inside the magnified box can carry a label there
            vis = m[m.BETA_pd.between(-0.05 * ihi, ihi) &
                    m.BETA_lbd_al.between(ilo_y, ihi)]
            place_labels(axin, label_items(pick_labels(vis, args.inset_genes, args.lbd_p,
                                                       tag="inset: ")),
                         args.inset_label_size or max(args.label_size - 2, 5), lw=0.4)

    to_label = pick_labels(m, args.genes, args.lbd_p, tag="main: ")
    inside = to_label[(to_label.BETA_pd.between(xlo, xhi)) &
                      (to_label.BETA_lbd_al.between(ylo, yhi))]

    labels = label_items(inside)
    off = m[(m.BETA_pd > xhi) | (m.BETA_pd < xlo) |
            (m.BETA_lbd_al > yhi) | (m.BETA_lbd_al < ylo)]
    if len(off):
        # every off-scale locus gets a marker; a label goes only on a selected signal
        # (from pick_labels) that fell outside the view, so labelling follows the same
        # rule as in-view points and no gene is labelled twice
        named = set(to_label.index.intersection(off.index))
        print(f"Off-scale in this view, marked on the frame: "
              f"{', '.join(sorted(off.GENE))}"
              + (f" (named: {', '.join(sorted(off.loc[list(named), 'GENE']))})"
                 if named else " (none named)"))
        labels += draw_offscale(ax, off, named, colours, xlo, xhi, ylo, yhi)
        # legend entry explaining the triangles
        ax.scatter([], [], marker="^", s=70, facecolor="none", edgecolor="black",
                   linewidth=0.5, label=f"off-scale locus (n={len(off)})")
    place_labels(ax, labels, args.label_size)

    # the axes are different phenotypes: a cross-disorder effect comparison,
    # not a discovery/replication design
    ax.set_xlabel(f"Beta, {label1}")
    ax.set_ylabel(f"Beta, {args.label2}")
    ax.legend(fontsize=8, loc=args.legend_loc)

    fig.tight_layout(); fig.savefig(args.out, dpi=300)
    print(f"Saved {args.out}")

    if args.table_out:
        m.sort_values("P_pd")[PLOT_COLS].to_csv(args.table_out, sep="\t", index=False)
        print(f"Saved loci table to {args.table_out}")
    if args.above_out:
        above[PLOT_COLS + ["delta_LBD_minus_PD"]].to_csv(args.above_out, sep="\t", index=False)
        print(f"Saved above-line table to {args.above_out}")


if __name__ == "__main__":
    main()
