#!/usr/bin/env python3
"""
Phase 3 paper figures for TNG-Cluster metallicity gradient study.

  Fig 9:  CNO nucleosynthetic channel diagnostics
  Fig 10: Comparison with observational constraints
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
    print(f"Loaded {len(df_all)} rows ({df_all['redshift'].nunique()} "
          f"redshifts), {len(df_z0)} at z=0")
    return df_all, df_z0


# ═══════════════════════════════════════════════════════════════════════════
# Figure 9: CNO nucleosynthetic channel diagnostics
# ═══════════════════════════════════════════════════════════════════════════

def fig9_cno_diagnostics(df_z0, df_all):
    """CNO abundance ratio gradients as nucleosynthetic probes."""
    fig, axes = plt.subplots(2, 2, figsize=(7.014, 5.5))

    cno_ratios = [
        ("slope_NO_star", r"$\nabla[\mathrm{N/O}]_\star$", WONG[0]),
        ("slope_CO_star", r"$\nabla[\mathrm{C/O}]_\star$", WONG[1]),
        ("slope_CN_star", r"$\nabla[\mathrm{C/N}]_\star$", WONG[2]),
    ]

    # --- (a) Overlaid histograms of CNO gradient distributions ---
    ax = axes[0, 0]
    bin_edges = np.linspace(-0.12, 0.04, 40)
    for col, label, color in cno_ratios:
        vals = df_z0[col].dropna().values
        med = np.median(vals)
        ax.hist(vals, bins=bin_edges, histtype="stepfilled", alpha=0.35,
                color=color, edgecolor=color, linewidth=0.8, label=label)
        ax.axvline(med, color=color, ls="--", lw=1.0, alpha=0.8)

    ax.set_xlabel(r"Gradient slope [dex/$R_e$]")
    ax.set_ylabel("Count")
    ax.legend(loc="upper left")
    add_minor_ticks(ax)

    # --- (b) [C/N] gradient vs [Fe/H] gradient ---
    ax = axes[0, 1]
    mask = df_z0["slope_CN_star"].notna() & df_z0["slope_FeH_star"].notna()
    x = df_z0.loc[mask, "slope_FeH_star"].values
    y = df_z0.loc[mask, "slope_CN_star"].values

    ax.scatter(x, y, s=2, alpha=0.25, color="0.5", rasterized=True,
               edgecolors="none")

    cx, my, lo, hi = running_median(x, y, n_bins=8, min_count=15)
    valid = np.isfinite(my)
    ax.plot(cx[valid], my[valid], color=WONG[1], lw=1.5, zorder=5)
    ax.fill_between(cx[valid], lo[valid], hi[valid],
                    color=WONG[1], alpha=0.2, zorder=4)

    ax.set_xlabel(r"$\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]")
    ax.set_ylabel(r"$\nabla[\mathrm{C/N}]_\star$ [dex/$R_e$]")
    ax.axhline(0, color="0.5", lw=0.4, ls=":")
    ax.axvline(0, color="0.5", lw=0.4, ls=":")
    add_minor_ticks(ax)

    # --- (c) CNO gradients by mass bin (grouped bars) ---
    ax = axes[1, 0]
    x_pos = np.arange(len(MASS_BINS))
    n_ratios = len(cno_ratios)
    total_width = 0.75
    bar_width = total_width / n_ratios

    for j, (col, label, color) in enumerate(cno_ratios):
        meds, errs_lo, errs_hi = [], [], []
        for mlo, mhi, _, _ in MASS_BINS:
            sel = df_z0[(df_z0["M_star_Msun"] >= mlo)
                        & (df_z0["M_star_Msun"] < mhi)]
            med, lo, hi = bootstrap_median_ci(sel[col].values)
            meds.append(med)
            errs_lo.append(med - lo)
            errs_hi.append(hi - med)
        meds = np.array(meds)
        offset = (j - (n_ratios - 1) / 2) * bar_width
        ax.bar(x_pos + offset, meds, bar_width, label=label,
               color=color, alpha=0.8, edgecolor="0.3", linewidth=0.4)
        ax.errorbar(x_pos + offset, meds,
                    yerr=[errs_lo, errs_hi],
                    fmt="none", ecolor="0.2", elinewidth=0.7,
                    capsize=2, capthick=0.5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([b[2] for b in MASS_BINS])
    ax.set_xlabel(r"$M_\star$ bin [$\mathrm{M}_\odot$]")
    ax.set_ylabel(r"Median gradient [dex/$R_e$]")
    ax.legend(loc="lower right", ncol=1)
    ax.axhline(0, color="0.5", lw=0.4, ls=":")
    add_minor_ticks(ax)

    # --- (d) Redshift evolution of CNO gradient medians ---
    ax = axes[1, 1]
    snaps_with_cno = df_all.groupby("snap_num").apply(
        lambda g: g["slope_NO_star"].notna().sum() > 10
    )
    valid_snaps = snaps_with_cno[snaps_with_cno].index

    for col, label, color in cno_ratios:
        zz, med_vals, lo_vals, hi_vals = [], [], [], []
        for snap in sorted(valid_snaps, reverse=True):
            sub = df_all[df_all["snap_num"] == snap]
            z_val = sub["redshift"].iloc[0]
            vals = sub[col].dropna().values
            if len(vals) < 10:
                continue
            med, lo, hi = bootstrap_median_ci(vals)
            zz.append(z_val)
            med_vals.append(med)
            lo_vals.append(lo)
            hi_vals.append(hi)

        zz = np.array(zz)
        med_vals = np.array(med_vals)
        lo_vals = np.array(lo_vals)
        hi_vals = np.array(hi_vals)
        order = np.argsort(zz)
        zz, med_vals, lo_vals, hi_vals = (
            zz[order], med_vals[order], lo_vals[order], hi_vals[order])

        ax.plot(zz, med_vals, "o-", color=color, ms=4, lw=1.2, label=label)
        ax.fill_between(zz, lo_vals, hi_vals, color=color, alpha=0.15)

    ax.set_xlabel("Redshift $z$")
    ax.set_ylabel(r"Median gradient [dex/$R_e$]")
    ax.axhline(0, color="0.5", lw=0.4, ls=":")
    ax.legend(loc="lower left")
    ax.set_xlim(-0.1, max(zz) + 0.15 if len(zz) else 1.0)
    add_minor_ticks(ax)

    for i, letter in enumerate(["a", "b", "c", "d"]):
        ax = axes.flat[i]
        ax.text(0.04, 0.96, f"({letter})",
                transform=ax.transAxes, fontweight="bold",
                va="top", ha="left")

    fig.subplots_adjust(hspace=0.35, wspace=0.30)
    savefig(fig, "phase3_fig09_cno_diagnostics")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 10: Comparison with observational constraints
# ═══════════════════════════════════════════════════════════════════════════

# Verified literature values for [Z/H] or [Fe/H] gradients.
# Parikh+21 Table 2: ETG [Z/H] within 1 Re, dex/Re (linear).
# Santucci+20: centrals D[Z/H]/Dlog(R/Re)=-0.25 -> ~-0.10 dex/Re (linear).
# Cheng+24: intrinsic (PSF-corrected) D[Fe/H]/Dlog(R/Re)=-0.151 -> ~-0.07 dex/Re.
# Cheng+25: JWST-SUSPENSE 8 galaxies, flat [Fe/H] gradients.

OBS_MASS = {
    r"Parikh+21 (MaNGA)": {
        "marker": "s", "color": "#CC79A7", "ms": 7, "zorder": 10,
        "logM": [10.05, 10.45, 10.70, 10.90, 11.15],
        "grad": [-0.13, -0.06, -0.09, -0.19, -0.17],
        "err":  [0.03,  0.02,  0.04,  0.03,  0.04],
    },
}

OBS_REDSHIFT = [
    {"label": r"Parikh+21 (MaNGA)",  "z": 0.03, "grad": -0.11, "err": 0.05,
     "marker": "s", "color": "#CC79A7", "ms": 8},
    {"label": r"Santucci+20 (SAMI)", "z": 0.05, "grad": -0.10, "err": 0.01,
     "marker": "D", "color": "#56B4E9", "ms": 7},
    {"label": r"Cheng+24 (LEGA-C)", "z": 0.80, "grad": -0.07, "err": 0.02,
     "marker": "^", "color": "#009E73", "ms": 8},
    {"label": r"Cheng+25 (SUSPENSE)", "z": 1.70, "grad": 0.00, "err": 0.10,
     "marker": "p", "color": "#E69F00", "ms": 8},
]


def fig10_obs_comparison(df_z0, df_all):
    """Comparison of TNG-Cluster predictions with observational data."""
    fig, axes = plt.subplots(1, 2, figsize=(7.014, 3.2))

    # --- (a) [Fe/H] gradient vs stellar mass at z=0 ---
    ax = axes[0]
    mask = df_z0["slope_FeH_star"].notna() & df_z0["M_star_Msun"].notna()
    logm = np.log10(df_z0.loc[mask, "M_star_Msun"].values)
    grad = df_z0.loc[mask, "slope_FeH_star"].values

    ax.scatter(logm, grad, s=2, alpha=0.15, color="0.6", rasterized=True,
               edgecolors="none", label="TNG-Cluster ($z=0$)")

    cx, my, lo, hi = running_median(logm, grad, n_bins=10, min_count=15)
    valid = np.isfinite(my)
    ax.plot(cx[valid], my[valid], color="0.15", lw=2.0, zorder=6,
            label="TNG-Cluster median")
    ax.fill_between(cx[valid], lo[valid], hi[valid],
                    color="0.15", alpha=0.15, zorder=5)

    for name, obs in OBS_MASS.items():
        ax.errorbar(obs["logM"], obs["grad"], yerr=obs["err"],
                    fmt=obs["marker"], color=obs["color"], ms=obs["ms"],
                    capsize=3, capthick=1.0, elinewidth=1.0,
                    markeredgecolor="0.2", markeredgewidth=0.5,
                    zorder=obs["zorder"], label=name)

    ax.set_xlabel(r"$\log\,(M_\star / \mathrm{M}_\odot)$")
    ax.set_ylabel(r"$\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]")
    ax.set_xlim(9.2, 12.2)
    ax.set_ylim(-0.7, 0.15)
    ax.axhline(0, color="0.5", lw=0.4, ls=":")
    ax.legend(loc="upper left")
    add_minor_ticks(ax)

    # Annotation about methodology caveat
    ax.text(0.96, 0.04,
            r"Obs: $\nabla[\mathrm{Z/H}]$",
            transform=ax.transAxes, va="bottom", ha="right",
            color="0.4", style="italic")

    # --- (b) Median [Fe/H] gradient vs redshift ---
    ax = axes[1]

    # TNG-Cluster track
    zz, med_vals, lo_vals, hi_vals = [], [], [], []
    for snap in sorted(df_all["snap_num"].unique(), reverse=True):
        sub = df_all[df_all["snap_num"] == snap]
        z_val = sub["redshift"].iloc[0]
        vals = sub["slope_FeH_star"].dropna().values
        if len(vals) < 20:
            continue
        med, lo, hi = bootstrap_median_ci(vals)
        zz.append(z_val)
        med_vals.append(med)
        lo_vals.append(lo)
        hi_vals.append(hi)

    zz = np.array(zz)
    med_vals = np.array(med_vals)
    lo_vals = np.array(lo_vals)
    hi_vals = np.array(hi_vals)
    order = np.argsort(zz)
    zz, med_vals, lo_vals, hi_vals = (
        zz[order], med_vals[order], lo_vals[order], hi_vals[order])

    ax.plot(zz, med_vals, "o-", color="0.15", ms=4, lw=1.5,
            label="TNG-Cluster", zorder=5)
    ax.fill_between(zz, lo_vals, hi_vals, color="0.15", alpha=0.15,
                    zorder=4)

    for obs in OBS_REDSHIFT:
        ax.errorbar(obs["z"], obs["grad"], yerr=obs["err"],
                    fmt=obs["marker"], color=obs["color"], ms=obs["ms"],
                    capsize=3, capthick=1.0, elinewidth=1.0,
                    markeredgecolor="0.2", markeredgewidth=0.5,
                    zorder=10, label=obs["label"])

    ax.set_xlabel("Redshift $z$")
    ax.set_ylabel(r"Median $\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]")
    ax.set_xlim(-0.1, max(zz) + 0.15 if len(zz) else 1.0)
    ax.set_ylim(-0.45, 0.15)
    ax.axhline(0, color="0.5", lw=0.4, ls=":")
    ax.legend(loc="lower left")
    add_minor_ticks(ax)

    for i, letter in enumerate(["a", "b"]):
        ax = axes[i]
        ax.text(0.04, 0.96, f"({letter})",
                transform=ax.transAxes, fontweight="bold",
                va="top", ha="left")

    fig.subplots_adjust(wspace=0.30)
    savefig(fig, "phase3_fig10_obs_comparison")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    df_all, df_z0 = load_data()

    print("\n--- Figure 9: CNO diagnostics ---")
    fig9_cno_diagnostics(df_z0, df_all)

    print("\n--- Figure 10: Observational comparison ---")
    fig10_obs_comparison(df_z0, df_all)

    # ── Summary statistics ──
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS FOR PAPER")
    print("=" * 60)

    print("\n--- CNO gradient distributions at z=0 ---")
    for col, name in [("slope_NO_star", "[N/O]"),
                      ("slope_CO_star", "[C/O]"),
                      ("slope_CN_star", "[C/N]"),
                      ("slope_MgFe_star", "[Mg/Fe]")]:
        v = df_z0[col].dropna().values
        pct_neg = 100 * np.mean(v < 0)
        print(f"  {name}: median={np.median(v):.4f}, "
              f"16-84%=[{np.percentile(v,16):.4f}, {np.percentile(v,84):.4f}], "
              f"{pct_neg:.1f}% negative")

    print("\n--- CNO mass correlations ---")
    for col, name in [("slope_NO_star", "[N/O]"),
                      ("slope_CO_star", "[C/O]"),
                      ("slope_CN_star", "[C/N]"),
                      ("slope_MgFe_star", "[Mg/Fe]")]:
        mask = df_z0[col].notna() & df_z0["M_star_Msun"].notna()
        rho, p = sp_stats.spearmanr(
            np.log10(df_z0.loc[mask, "M_star_Msun"]), df_z0.loc[mask, col])
        print(f"  log(M*) vs {name}: rho={rho:.3f}, p={p:.2e}")

    print("\n--- [C/N] vs [Fe/H] gradient correlation ---")
    mask = df_z0["slope_CN_star"].notna() & df_z0["slope_FeH_star"].notna()
    rho, p = sp_stats.spearmanr(
        df_z0.loc[mask, "slope_FeH_star"], df_z0.loc[mask, "slope_CN_star"])
    print(f"  rho={rho:.3f}, p={p:.2e}")

    print("\n--- CNO inner/outer comparison ---")
    for name, icol, ocol in [
        ("[N/O]", "slope_inner_NO_star", "slope_outer_NO_star"),
        ("[C/O]", "slope_inner_CO_star", "slope_outer_CO_star"),
        ("[C/N]", "slope_inner_CN_star", "slope_outer_CN_star"),
    ]:
        inner = df_z0[icol].dropna().values
        outer = df_z0[ocol].dropna().values
        print(f"  {name}: inner med={np.median(inner):.4f}, "
              f"outer med={np.median(outer):.4f}, "
              f"ratio={np.median(inner)/np.median(outer):.2f}")

    print("\nDone.")
