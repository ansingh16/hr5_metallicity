#!/usr/bin/env python3
"""
Phase 2 paper figures for TNG-Cluster metallicity gradient study.

  Fig 6: Inner/outer gradient decomposition
  Fig 7: Environmental dependence (cluster-centric distance)
  Fig 8: Velocity dispersion and stellar age correlations
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import Normalize
from matplotlib import cm
from scipy import stats as sp_stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

STYLE_PATH = os.path.join(os.path.dirname(__file__), "..", "paper.mplstyle")
plt.style.use(STYLE_PATH)

BASEDIR = os.path.join(os.path.dirname(__file__), "..")
CATDIR = os.path.join(BASEDIR, cfg.CATALOG_DIR)
PLOTDIR = os.path.join(BASEDIR, cfg.PLOT_DIR, "paper")

MASS_BINS = [
    (3.162e9,  1e10,    r"$10^{9.5}\!-\!10^{10}$",   "#0072B2"),
    (1e10,     3.162e10, r"$10^{10}\!-\!10^{10.5}$",  "#D55E00"),
    (3.162e10, 1e11,    r"$10^{10.5}\!-\!10^{11}$",   "#009E73"),
    (1e11,     3.162e11, r"$10^{11}\!-\!10^{11.5}$",  "#E69F00"),
    (3.162e11, 1e13,    r"$>10^{11.5}$",              "#56B4E9"),
]

WONG = ["#0072B2", "#D55E00", "#009E73", "#E69F00",
        "#56B4E9", "#CC79A7", "#F0E442"]


def bootstrap_median_ci(values, n_boot=1000, ci=68):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 5:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(42)
    medians = np.array([
        np.median(rng.choice(values, size=len(values), replace=True))
        for _ in range(n_boot)
    ])
    lo = np.percentile(medians, (100 - ci) / 2)
    hi = np.percentile(medians, 100 - (100 - ci) / 2)
    return np.median(values), lo, hi


def running_median(x, y, n_bins=8, min_count=10):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < min_count:
        return np.array([]), np.array([]), np.array([]), np.array([])
    bins = np.linspace(np.percentile(x, 2), np.percentile(x, 98), n_bins + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    meds, lo, hi = [], [], []
    for i in range(len(bins) - 1):
        m = (x >= bins[i]) & (x < bins[i + 1])
        if m.sum() >= min_count:
            meds.append(np.median(y[m]))
            lo.append(np.percentile(y[m], 16))
            hi.append(np.percentile(y[m], 84))
        else:
            meds.append(np.nan)
            lo.append(np.nan)
            hi.append(np.nan)
    return centers, np.array(meds), np.array(lo), np.array(hi)


def add_minor_ticks(ax):
    if ax.get_xscale() == "linear":
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    if ax.get_yscale() == "linear":
        ax.yaxis.set_minor_locator(AutoMinorLocator())


def savefig(fig, name):
    os.makedirs(PLOTDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(PLOTDIR, f"{name}.{ext}"),
                    bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {name}")


def load_data():
    evo_file = os.path.join(CATDIR, "gradients_evolution.csv")
    df_all = pd.read_csv(evo_file)
    df_z0 = df_all[np.isclose(df_all["redshift"], 0.0, atol=0.01)].copy()

    ref_file = os.path.join(CATDIR, "elliptical_catalog_snap099.csv")
    if os.path.exists(ref_file):
        ref = pd.read_csv(ref_file)
        for col in ["SFR", "sSFR"]:
            if col not in df_z0.columns and col in ref.columns:
                merge = ref[["subhalo_id", col]]
                df_z0 = df_z0.merge(merge, on="subhalo_id", how="left")

    print(f"Loaded {len(df_all)} rows ({df_all['redshift'].nunique()} "
          f"redshifts), {len(df_z0)} at z=0")
    return df_all, df_z0


# ═══════════════════════════════════════════════════════════════════════════
# Figure 6: Inner/outer gradient decomposition
# ═══════════════════════════════════════════════════════════════════════════

def fig6_inner_outer(df_z0):
    """Inner vs outer gradients: radial breaks reveal in-situ vs accreted."""
    from scipy.stats import gaussian_kde

    fig, axes = plt.subplots(2, 2, figsize=(7.014, 5.5))

    # --- Top row: KDE density scatter + contours ---
    for ax, inner_col, outer_col, label in [
        (axes[0, 0], "slope_inner_FeH_star", "slope_outer_FeH_star",
         r"$\nabla[\mathrm{Fe/H}]_\star$"),
        (axes[0, 1], "slope_inner_OFe_star", "slope_outer_OFe_star",
         r"$\nabla[\mathrm{O/Fe}]_\star$"),
    ]:
        mask = df_z0[inner_col].notna() & df_z0[outer_col].notna()
        x = df_z0.loc[mask, outer_col].values
        y = df_z0.loc[mask, inner_col].values

        xy = np.vstack([x, y])
        kde = gaussian_kde(xy)
        density = kde(xy)

        idx = density.argsort()
        x_s, y_s, d_s = x[idx], y[idx], density[idx]

        ax.scatter(x_s, y_s, c=d_s, cmap="Greys", s=8, alpha=0.8,
                   rasterized=True, edgecolors="none", zorder=1)

        xg = np.linspace(x.min(), x.max(), 80)
        yg = np.linspace(y.min(), y.max(), 80)
        Xg, Yg = np.meshgrid(xg, yg)
        Zg = kde(np.vstack([Xg.ravel(), Yg.ravel()])).reshape(Xg.shape)
        levels = np.percentile(density, [30, 68, 95])
        ax.contour(Xg, Yg, Zg, levels=sorted(levels),
                   colors="r", linewidths=0.6, alpha=0.8, zorder=2)

        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "r--", lw=0.8, alpha=0.6, zorder=3)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        ax.set_xlabel(f"Outer {label}" + r" [dex/$R_e$]")
        ax.set_ylabel(f"Inner {label}" + r" [dex/$R_e$]")
        add_minor_ticks(ax)

    # --- Bottom row: grouped bars inner/outer by mass bin ---
    for ax, inner_col, outer_col, ylabel in [
        (axes[1, 0], "slope_inner_FeH_star", "slope_outer_FeH_star",
         r"Median $\nabla[\mathrm{Fe/H}]_\star$"),
        (axes[1, 1], "slope_inner_OFe_star", "slope_outer_OFe_star",
         r"Median $\nabla[\mathrm{O/Fe}]_\star$"),
    ]:
        x_pos = np.arange(len(MASS_BINS))
        width = 0.35

        for offset, col, lbl, fc in [
            (-width / 2, inner_col, "Inner", WONG[0]),
            (+width / 2, outer_col, "Outer", WONG[1]),
        ]:
            meds, errs_lo, errs_hi = [], [], []
            for mlo, mhi, _, _ in MASS_BINS:
                sel = df_z0[(df_z0["M_star_Msun"] >= mlo)
                            & (df_z0["M_star_Msun"] < mhi)]
                med, lo, hi = bootstrap_median_ci(sel[col].values)
                meds.append(med)
                errs_lo.append(med - lo)
                errs_hi.append(hi - med)
            meds = np.array(meds)
            ax.bar(x_pos + offset, meds, width, label=lbl,
                   color=fc, alpha=0.8, edgecolor="0.3", linewidth=0.4)
            ax.errorbar(x_pos + offset, meds,
                        yerr=[errs_lo, errs_hi],
                        fmt="none", ecolor="0.2", elinewidth=0.7,
                        capsize=2, capthick=0.5)

        ax.set_xticks(x_pos)
        ax.set_xticklabels([b[2] for b in MASS_BINS])
        ax.set_xlabel(r"$M_\star\;[M_\odot]$")
        ax.set_ylabel(ylabel + r" [dex/$R_e$]")
        if ax is axes[1, 1]:
            ax.legend(loc="upper right")
        ax.axhline(0, color="0.5", lw=0.4, ls=":")
        add_minor_ticks(ax)

    for i, letter in enumerate(["a", "b", "c", "d"]):
        ax = axes.flat[i]
        ax.text(0.96, 0.96 if i < 2 else 0.04, f"({letter})",
                transform=ax.transAxes, fontweight="bold",
                va="top" if i < 2 else "bottom", ha="right")

    fig.subplots_adjust(hspace=0.35, wspace=0.3)
    savefig(fig, "phase2_fig06_inner_outer")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 7: Environmental dependence (R/R200)
# ═══════════════════════════════════════════════════════════════════════════

def fig7_environment(df_z0):
    """Cluster-centric distance effects on gradients."""
    fig, axes = plt.subplots(2, 2, figsize=(7.014, 5.5))

    # --- Top row: gradient vs R/R200 scatter + running median ---
    for ax, col, ylabel, ylim in [
        (axes[0, 0], "slope_FeH_star",
         r"$\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]",
         (-0.7, 0.15)),
        (axes[0, 1], "slope_OFe_star",
         r"$\nabla[\mathrm{O/Fe}]_\star$ [dex/$R_e$]",
         (-0.15, 0.25)),
    ]:
        mask = df_z0[col].notna() & df_z0["R_over_R200"].notna()
        x = df_z0.loc[mask, "R_over_R200"].values
        y = df_z0.loc[mask, col].values

        ax.scatter(x, y, s=2, alpha=0.25, color="0.5", rasterized=True,
                   edgecolors="none")

        # Running median
        cx, my, lo, hi = running_median(x, y, n_bins=10, min_count=15)
        valid = np.isfinite(my)
        ax.plot(cx[valid], my[valid], color=WONG[0], lw=1.5, zorder=5)
        ax.fill_between(cx[valid], lo[valid], hi[valid],
                        color=WONG[0], alpha=0.2, zorder=4)

        ax.set_xlabel(r"$R/R_{200}$")
        ax.set_ylabel(ylabel)
        ax.set_ylim(ylim)
        ax.set_xlim(-0.1, 5.0)
        ax.axhline(0, color="0.5", lw=0.4, ls=":")
        add_minor_ticks(ax)

    # --- Bottom left: boxplot by radial bin ---
    radial_bins = [(0, 0.5, r"$<0.5$"),
                   (0.5, 1.0, r"$0.5\!-\!1$"),
                   (1.0, 2.0, r"$1\!-\!2$"),
                   (2.0, 5.0, r"$>2$")]

    ax = axes[1, 0]
    box_data = []
    box_labels = []
    box_n = []
    for rlo, rhi, lbl in radial_bins:
        sel = df_z0[(df_z0["R_over_R200"] >= rlo)
                    & (df_z0["R_over_R200"] < rhi)]
        vals = sel["slope_FeH_star"].dropna().values
        box_data.append(vals)
        box_labels.append(lbl)
        box_n.append(len(vals))

    bp = ax.boxplot(box_data, positions=range(len(box_data)),
                    widths=0.5, patch_artist=True,
                    showfliers=False, medianprops=dict(color="k", lw=1.0))
    for patch, color in zip(bp["boxes"],
                            [WONG[0], WONG[1], WONG[2], WONG[3]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    ax.set_xticks(range(len(box_labels)))
    ax.set_xticklabels(box_labels)
    ax.set_xlabel(r"$R/R_{200}$ bin")
    ax.set_ylabel(r"$\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]")
    ax.axhline(0, color="0.5", lw=0.4, ls=":")
    for i, n in enumerate(box_n):
        ax.text(i, ax.get_ylim()[1] * 0.95, f"$N$={n}",
                ha="center", color="0.4")
    add_minor_ticks(ax)

    # --- Bottom right: inner/outer [Fe/H] gradient by radial bin ---
    ax = axes[1, 1]
    x_pos = np.arange(len(radial_bins))
    width = 0.35

    for offset, col, lbl, fc in [
        (-width / 2, "slope_inner_FeH_star", "Inner", WONG[0]),
        (+width / 2, "slope_outer_FeH_star", "Outer", WONG[1]),
    ]:
        meds, errs_lo, errs_hi = [], [], []
        for rlo, rhi, _ in radial_bins:
            sel = df_z0[(df_z0["R_over_R200"] >= rlo)
                        & (df_z0["R_over_R200"] < rhi)]
            med, lo, hi = bootstrap_median_ci(sel[col].values)
            meds.append(med)
            errs_lo.append(med - lo)
            errs_hi.append(hi - med)
        meds = np.array(meds)
        ax.bar(x_pos + offset, meds, width, label=lbl,
               color=fc, alpha=0.8, edgecolor="0.3", linewidth=0.4)
        ax.errorbar(x_pos + offset, meds,
                    yerr=[errs_lo, errs_hi],
                    fmt="none", ecolor="0.2", elinewidth=0.7,
                    capsize=2, capthick=0.5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([b[2] for b in radial_bins])
    ax.set_xlabel(r"$R/R_{200}$ bin")
    ax.set_ylabel(r"Median $\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]")
    ax.legend(loc="lower left")
    ax.axhline(0, color="0.5", lw=0.4, ls=":")
    add_minor_ticks(ax)

    for i, letter in enumerate(["a", "b", "c", "d"]):
        ax = axes.flat[i]
        ax.text(0.04, 0.96, f"({letter})",
                transform=ax.transAxes, fontweight="bold",
                va="top", ha="left")

    fig.subplots_adjust(hspace=0.35, wspace=0.3)
    savefig(fig, "phase2_fig07_environment")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 8: Velocity dispersion and stellar age
# ═══════════════════════════════════════════════════════════════════════════

def fig8_sigma_age(df_z0):
    """Gradient correlations with velocity dispersion and stellar age."""
    fig, axes = plt.subplots(2, 2, figsize=(7.014, 5.5))

    panels = [
        (axes[0, 0], "vel_disp_kms", "slope_FeH_star",
         r"$\sigma_\star$ [km s$^{-1}$]",
         r"$\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]",
         (-0.7, 0.15), (30, 500)),
        (axes[0, 1], "vel_disp_kms", "slope_OFe_star",
         r"$\sigma_\star$ [km s$^{-1}$]",
         r"$\nabla[\mathrm{O/Fe}]_\star$ [dex/$R_e$]",
         (-0.15, 0.25), (30, 500)),
        (axes[1, 0], "age_star_global", "slope_FeH_star",
         r"$\langle\mathrm{Age}\rangle_\star$ [Gyr]",
         r"$\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]",
         (-0.7, 0.15), None),
        (axes[1, 1], "age_star_global", "slope_OFe_star",
         r"$\langle\mathrm{Age}\rangle_\star$ [Gyr]",
         r"$\nabla[\mathrm{O/Fe}]_\star$ [dex/$R_e$]",
         (-0.15, 0.25), None),
    ]

    for ax, xcol, ycol, xlabel, ylabel, ylim, xlim in panels:
        mask = df_z0[xcol].notna() & df_z0[ycol].notna()
        x = df_z0.loc[mask, xcol].values
        y = df_z0.loc[mask, ycol].values

        ax.scatter(x, y, s=2, alpha=0.25, color="0.5", rasterized=True,
                   edgecolors="none")

        # Running median
        cx, my, lo, hi = running_median(x, y, n_bins=8, min_count=15)
        valid = np.isfinite(my)
        ax.plot(cx[valid], my[valid], color=WONG[1], lw=1.5, zorder=5)
        ax.fill_between(cx[valid], lo[valid], hi[valid],
                        color=WONG[1], alpha=0.2, zorder=4)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_ylim(ylim)
        if xlim is not None:
            ax.set_xlim(xlim)
            ax.set_xscale("log")
        ax.axhline(0, color="0.5", lw=0.4, ls=":")
        add_minor_ticks(ax)

    for i, letter in enumerate(["a", "b", "c", "d"]):
        ax = axes.flat[i]
        ax.text(0.04, 0.96, f"({letter})",
                transform=ax.transAxes, fontweight="bold",
                va="top", ha="left")

    fig.subplots_adjust(hspace=0.35, wspace=0.3)
    savefig(fig, "phase2_fig08_sigma_age")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    df_all, df_z0 = load_data()

    print("\n--- Figure 6: Inner/outer decomposition ---")
    fig6_inner_outer(df_z0)

    print("\n--- Figure 7: Environmental dependence ---")
    fig7_environment(df_z0)

    print("\n--- Figure 8: Velocity dispersion & age ---")
    fig8_sigma_age(df_z0)

    # Print summary statistics for paper text
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS FOR PAPER")
    print("=" * 60)

    print("\n--- Inner/outer decomposition ---")
    for name, icol, ocol in [
        ("[Fe/H]", "slope_inner_FeH_star", "slope_outer_FeH_star"),
        ("[O/Fe]", "slope_inner_OFe_star", "slope_outer_OFe_star"),
    ]:
        inner = df_z0[icol].dropna().values
        outer = df_z0[ocol].dropna().values
        n = min(len(inner), len(outer))
        stat, pval = sp_stats.wilcoxon(inner[:n], outer[:n])
        print(f"  {name}: inner median={np.median(inner):.3f}, "
              f"outer median={np.median(outer):.3f}, "
              f"ratio={np.median(inner)/np.median(outer):.2f}, "
              f"Wilcoxon p={pval:.1e}")

    print("\n--- Environmental correlations ---")
    for col, name in [("slope_FeH_star", "[Fe/H]"),
                      ("slope_OFe_star", "[O/Fe]")]:
        mask = df_z0[col].notna() & df_z0["R_over_R200"].notna()
        rho, p = sp_stats.spearmanr(
            df_z0.loc[mask, "R_over_R200"], df_z0.loc[mask, col])
        print(f"  R/R200 vs {name}: rho={rho:.3f}, p={p:.1e}")

    r_split = 0.5
    inner_cl = df_z0[df_z0["R_over_R200"] < r_split]
    outer_cl = df_z0[df_z0["R_over_R200"] >= r_split]
    print(f"  Inner cluster (R<{r_split}R200): N={len(inner_cl)}, "
          f"median [Fe/H] grad={np.nanmedian(inner_cl.slope_FeH_star):.3f}")
    print(f"  Outer cluster (R>={r_split}R200): N={len(outer_cl)}, "
          f"median [Fe/H] grad={np.nanmedian(outer_cl.slope_FeH_star):.3f}")
    stat, p = sp_stats.mannwhitneyu(
        inner_cl["slope_FeH_star"].dropna(),
        outer_cl["slope_FeH_star"].dropna())
    print(f"  Mann-Whitney U: p={p:.1e}")

    print("\n--- Velocity dispersion & age ---")
    for xcol, xname in [("vel_disp_kms", "sigma"),
                        ("age_star_global", "age")]:
        for ycol, yname in [("slope_FeH_star", "[Fe/H]"),
                            ("slope_OFe_star", "[O/Fe]")]:
            mask = df_z0[xcol].notna() & df_z0[ycol].notna()
            rho, p = sp_stats.spearmanr(
                df_z0.loc[mask, xcol], df_z0.loc[mask, ycol])
            print(f"  {xname} vs {yname}: rho={rho:.3f}, p={p:.1e}")

    print("\nDone.")
