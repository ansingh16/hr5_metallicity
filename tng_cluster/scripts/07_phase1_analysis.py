#!/usr/bin/env python3
"""
Phase 1 Analysis: Core results for TNG-Cluster metallicity gradient paper.

Analyses implemented:
  1. Gradient distributions at z=0 (histograms + summary statistics)
  2. Mass dependence of gradients (slope vs M_star)
  4. Redshift evolution of gradient slopes (all snapshots)
  7+8. Alpha-element gradients + [O/Fe] vs [Fe/H] diagnostic
  9. Gas depletion, SFR profiles, quenching

Outputs:
  catalogs/gradients_evolution.csv        — rebuilt combined catalog
  plots/phase1_fig01_z0_distributions.pdf — gradient distributions at z=0
  plots/phase1_fig02_mass_dependence.pdf  — gradient vs stellar mass
  plots/phase1_fig03_redshift_evolution.pdf — gradient evolution z=0-8
  plots/phase1_fig04_sample_size.pdf      — sample completeness
  plots/phase1_fig05_alpha_gradients.pdf  — [O/Fe], [Mg/Fe] gradient analysis
  plots/phase1_fig06_ofe_feh_tracks.pdf   — [O/Fe] vs [Fe/H] diagnostic
  plots/phase1_fig07_gas_quenching.pdf    — gas depletion + SFR evolution
  tables/phase1_table1_z0_stats.csv       — summary statistics at z=0
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy import stats as sp_stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

STYLE_PATH = os.path.join(os.path.dirname(__file__), "..", "paper.mplstyle")
plt.style.use(STYLE_PATH)

BASEDIR = os.path.join(os.path.dirname(__file__), "..")
CATDIR = os.path.join(BASEDIR, cfg.CATALOG_DIR)
PLOTDIR = os.path.join(BASEDIR, cfg.PLOT_DIR)
TABLEDIR = os.path.join(BASEDIR, "tables")

MASS_BINS = [
    (2e9,  1e10, r"$10^{9.3}\!-\!10^{10}$",           "#0072B2"),
    (1e10, 5e10, r"$10^{10}\!-\!10^{10.7}$",           "#D55E00"),
    (5e10, 2e11, r"$10^{10.7}\!-\!10^{11.3}$",         "#009E73"),
    (2e11, 1e13, r"$>10^{11.3}$",                       "#E69F00"),
]

WONG_COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00",
               "#56B4E9", "#CC79A7", "#F0E442"]


# ═══════════════════════════════════════════════════════════════════════════
# Utilities
# ═══════════════════════════════════════════════════════════════════════════

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
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())


# ═══════════════════════════════════════════════════════════════════════════
# Step 0: Rebuild combined catalog
# ═══════════════════════════════════════════════════════════════════════════

def rebuild_combined_catalog():
    """Combine all per-snapshot gradient catalogs into one file."""
    import glob

    files = sorted(glob.glob(os.path.join(CATDIR, "gradients_snap*.csv")))
    frames = []
    for f in files:
        snap = int(os.path.basename(f).replace("gradients_snap", "").replace(".csv", ""))
        df = pd.read_csv(f)
        if "snap_num" not in df.columns:
            df["snap_num"] = snap
        if "redshift" not in df.columns and snap in cfg.SNAPSHOTS:
            df["redshift"] = cfg.SNAPSHOTS[snap]
        frames.append(df)
        print(f"  snap{snap:03d} (z={cfg.SNAPSHOTS.get(snap, '?'):.2f}): "
              f"{len(df)} galaxies, {len(df.columns)} cols")

    combined = pd.concat(frames, ignore_index=True, sort=False)

    # For z=0 rows, merge environment info from elliptical catalog
    ref_file = os.path.join(CATDIR, "elliptical_catalog_snap099.csv")
    if os.path.exists(ref_file):
        ref = pd.read_csv(ref_file)
        env_cols = ["vel_disp_kms", "SFR", "sSFR", "cluster_M200_Msun",
                    "cluster_R200_pkpc"]
        ref_merge = ref[["subhalo_id"] + env_cols].copy()

        z0_mask = combined["redshift"] == 0.0
        for col in env_cols:
            if col not in combined.columns:
                combined[col] = np.nan

        z0_ids = combined.loc[z0_mask, "subhalo_id"].values
        ref_dict = ref_merge.set_index("subhalo_id").to_dict("index")
        for idx in combined.index[z0_mask]:
            sid = combined.at[idx, "subhalo_id"]
            if sid in ref_dict:
                for col in env_cols:
                    combined.at[idx, col] = ref_dict[sid].get(col, np.nan)

        # Compute R/R200 for z=0 if not present
        if "R_over_R200" not in combined.columns:
            combined["R_over_R200"] = np.nan

    out_file = os.path.join(CATDIR, "gradients_evolution.csv")
    combined.to_csv(out_file, index=False)
    print(f"\n  Combined: {len(combined)} rows, "
          f"{len(combined.columns)} columns -> {out_file}")
    return combined


# ═══════════════════════════════════════════════════════════════════════════
# Analysis 1: Gradient distributions at z=0
# ═══════════════════════════════════════════════════════════════════════════

def analysis1_z0_distributions(df_z0):
    """Multi-panel histograms + summary statistics table."""
    print("\n" + "=" * 60)
    print("Analysis 1: Gradient distributions at z=0")
    print("=" * 60)

    slope_cols = [
        ("slope_FeH_star",   r"$\nabla$[Fe/H]$_\star$",   "Stellar"),
        ("slope_Z_star",     r"$\nabla Z_\star$",          "Stellar"),
        ("slope_OFe_star",   r"$\nabla$[O/Fe]$_\star$",    "Stellar"),
        ("slope_MgFe_star",  r"$\nabla$[Mg/Fe]$_\star$",   "Stellar"),
        ("slope_FeH_gas",    r"$\nabla$[Fe/H]$_\mathrm{gas}$", "Gas"),
        ("slope_Z_gas",      r"$\nabla Z_\mathrm{gas}$",       "Gas"),
        ("slope_NO_star",    r"$\nabla$[N/O]$_\star$",     "Stellar"),
        ("slope_CO_star",    r"$\nabla$[C/O]$_\star$",     "Stellar"),
        ("slope_CN_star",    r"$\nabla$[C/N]$_\star$",     "Stellar"),
        ("slope_SFR_profile", r"$\nabla$SFR",              "Gas"),
    ]

    # --- Figure 1: Histograms ---
    fig, axes = plt.subplots(2, 5, figsize=(7.014, 4.0))
    axes_flat = axes.flatten()

    stats_rows = []
    for i, (col, label, phase) in enumerate(slope_cols):
        ax = axes_flat[i]
        vals = df_z0[col].dropna().values
        n_valid = len(vals)

        if n_valid > 5:
            med = np.median(vals)
            p16, p84 = np.percentile(vals, [16, 84])
            std = np.std(vals)
            frac_neg = np.mean(vals < 0) * 100
            frac_pos = np.mean(vals > 0) * 100

            ax.hist(vals, bins=30, color=WONG_COLORS[0] if phase == "Stellar"
                    else WONG_COLORS[1], alpha=0.7, edgecolor="white",
                    linewidth=0.3)
            ax.axvline(med, color="k", ls="-", lw=1.0)
            ax.axvline(0, color="grey", ls="--", lw=0.6)
            ax.set_title(label, fontsize=7, pad=2)

            stat_text = (f"med={med:.3f}\n"
                         f"N={n_valid}")
            ax.text(0.97, 0.95, stat_text, transform=ax.transAxes,
                    fontsize=5.5, va="top", ha="right",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec="none", alpha=0.8))

            stats_rows.append({
                "Gradient": col.replace("slope_", ""),
                "Phase": phase,
                "N": n_valid,
                "Median": f"{med:.4f}",
                "16th": f"{p16:.4f}",
                "84th": f"{p84:.4f}",
                "Std": f"{std:.4f}",
                "Neg%": f"{frac_neg:.1f}",
                "Pos%": f"{frac_pos:.1f}",
            })
        else:
            ax.set_visible(False)

    for ax in axes_flat:
        if ax.get_visible():
            ax.set_xlabel(r"Slope [dex/$R_e$]", fontsize=6)
            ax.set_ylabel("Count", fontsize=6)
            add_minor_ticks(ax)

    fig.suptitle(f"TNG-Cluster: Metallicity Gradient Distributions at z=0 "
                 f"(N={len(df_z0)})", fontsize=9, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig01_z0_distributions.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig01_z0_distributions.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved phase1_fig01_z0_distributions.pdf")

    # --- Table 1: Summary statistics ---
    stats_df = pd.DataFrame(stats_rows)
    stats_file = os.path.join(TABLEDIR, "phase1_table1_z0_stats.csv")
    stats_df.to_csv(stats_file, index=False)
    print(f"  Saved {stats_file}")
    print(stats_df.to_string(index=False))

    return stats_df


# ═══════════════════════════════════════════════════════════════════════════
# Analysis 2: Mass dependence of gradients
# ═══════════════════════════════════════════════════════════════════════════

def analysis2_mass_dependence(df_z0, df_all):
    """Gradient slope vs stellar mass at z=0, z=0.5, z=1.0, z=2.0."""
    print("\n" + "=" * 60)
    print("Analysis 2: Mass dependence of gradients")
    print("=" * 60)

    slope_keys = [
        ("slope_FeH_star",  r"$\nabla$[Fe/H]$_\star$"),
        ("slope_OFe_star",  r"$\nabla$[O/Fe]$_\star$"),
        ("slope_MgFe_star", r"$\nabla$[Mg/Fe]$_\star$"),
        ("slope_Z_star",    r"$\nabla Z_\star$"),
        ("slope_Z_gas",     r"$\nabla Z_\mathrm{gas}$"),
    ]

    z_panels = [0.0, 0.5, 1.0, 2.0]
    n_rows = len(z_panels)
    n_cols = len(slope_keys)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.014, 7.5),
                             sharex=True)

    for row, z_sel in enumerate(z_panels):
        sub = df_all[np.isclose(df_all["redshift"], z_sel, atol=0.05)]
        if len(sub) < 10:
            for col in range(n_cols):
                axes[row, col].text(0.5, 0.5, f"z={z_sel}\nN<10",
                                    transform=axes[row, col].transAxes,
                                    ha="center", va="center", fontsize=8)
            continue

        logm = np.log10(sub["M_star_Msun"].values)

        for col, (slope_col, ylabel) in enumerate(slope_keys):
            ax = axes[row, col]
            y = sub[slope_col].values
            valid = np.isfinite(y)

            if valid.sum() > 10:
                ax.scatter(logm[valid], y[valid], s=1.5, alpha=0.15,
                           color=WONG_COLORS[0], rasterized=True)

                c, med, lo, hi = running_median(logm, y, n_bins=7,
                                                min_count=8)
                good = np.isfinite(med)
                if good.any():
                    ax.fill_between(c[good], lo[good], hi[good],
                                    alpha=0.2, color=WONG_COLORS[1])
                    ax.plot(c[good], med[good], "o-",
                            color=WONG_COLORS[1], ms=3, lw=1.2)

                # Spearman correlation
                xv, yv = logm[valid], y[valid]
                rho, pval = sp_stats.spearmanr(xv, yv)
                ax.text(0.03, 0.05,
                        f"$\\rho$={rho:.2f}\nN={valid.sum()}",
                        transform=ax.transAxes, fontsize=5.5, va="bottom")

            ax.axhline(0, color="grey", ls="--", lw=0.6)

            if row == 0:
                ax.set_title(ylabel, fontsize=7)
            if row == n_rows - 1:
                ax.set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
            if col == 0:
                ax.set_ylabel(f"z={z_sel:.1f}\nSlope [dex/$R_e$]",
                              fontsize=7)
            add_minor_ticks(ax)

    fig.suptitle("Gradient Slope vs Stellar Mass", fontsize=10, y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig02_mass_dependence.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig02_mass_dependence.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved phase1_fig02_mass_dependence.pdf")


# ═══════════════════════════════════════════════════════════════════════════
# Analysis 4: Redshift evolution of gradient slopes
# ═══════════════════════════════════════════════════════════════════════════

def analysis4_redshift_evolution(df_all):
    """Gradient evolution with redshift, by mass bin."""
    print("\n" + "=" * 60)
    print("Analysis 4: Redshift evolution of gradient slopes")
    print("=" * 60)

    slope_keys = [
        ("slope_FeH_star",    r"$\nabla$[Fe/H]$_\star$"),
        ("slope_OFe_star",    r"$\nabla$[O/Fe]$_\star$"),
        ("slope_MgFe_star",   r"$\nabla$[Mg/Fe]$_\star$"),
        ("slope_Z_star",      r"$\nabla Z_\star$"),
        ("slope_FeH_gas",     r"$\nabla$[Fe/H]$_\mathrm{gas}$"),
        ("slope_Z_gas",       r"$\nabla Z_\mathrm{gas}$"),
    ]

    redshifts = sorted(df_all["redshift"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(7.014, 5.0))
    axes_flat = axes.flatten()

    for i, (slope_col, ylabel) in enumerate(slope_keys):
        ax = axes_flat[i]

        # All galaxies combined
        meds_all, lo_all, hi_all, zs_all = [], [], [], []
        for z in redshifts:
            sub = df_all.loc[np.isclose(df_all["redshift"], z, atol=0.01),
                             slope_col]
            med, lo, hi = bootstrap_median_ci(sub)
            if np.isfinite(med):
                meds_all.append(med)
                lo_all.append(lo)
                hi_all.append(hi)
                zs_all.append(z)

        if zs_all:
            ax.fill_between(zs_all, lo_all, hi_all, alpha=0.12, color="grey")
            ax.plot(zs_all, meds_all, "s-", color="grey", ms=3, lw=0.8,
                    alpha=0.5, label="All")

        # By mass bin
        for mlo, mhi, mlabel, mcolor in MASS_BINS:
            meds, los, his, zs = [], [], [], []
            for z in redshifts:
                mask = (np.isclose(df_all["redshift"], z, atol=0.01) &
                        (df_all["M_star_Msun"] >= mlo) &
                        (df_all["M_star_Msun"] < mhi))
                sub = df_all.loc[mask, slope_col]
                med, lo, hi = bootstrap_median_ci(sub)
                if np.isfinite(med):
                    meds.append(med); los.append(lo); his.append(hi)
                    zs.append(z)
            if zs:
                ax.fill_between(zs, los, his, alpha=0.12, color=mcolor)
                ax.plot(zs, meds, "o-", color=mcolor, ms=3, label=mlabel)

        ax.axhline(0, color="grey", ls="--", lw=0.6)
        ax.set_xlabel("Redshift")
        ax.set_ylabel(ylabel + r" [dex/$R_e$]", fontsize=7)
        ax.set_xlim(-0.2, max(redshifts) + 0.3)
        add_minor_ticks(ax)

        if i == 0:
            ax.legend(fontsize=5, loc="lower left", ncol=1)

    fig.suptitle("TNG-Cluster: Metallicity Gradient Evolution", fontsize=10,
                 y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig03_redshift_evolution.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig03_redshift_evolution.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved phase1_fig03_redshift_evolution.pdf")

    # --- Sample size figure ---
    fig, ax = plt.subplots(1, 1, figsize=(3.321, 2.5))
    for mlo, mhi, mlabel, mcolor in MASS_BINS:
        counts, zs = [], []
        for z in redshifts:
            mask = (np.isclose(df_all["redshift"], z, atol=0.01) &
                    (df_all["M_star_Msun"] >= mlo) &
                    (df_all["M_star_Msun"] < mhi))
            counts.append(mask.sum())
            zs.append(z)
        ax.plot(zs, counts, "o-", color=mcolor, ms=3, label=mlabel)

    # Total
    total_counts = [np.isclose(df_all["redshift"], z, atol=0.01).sum()
                    for z in redshifts]
    ax.plot(redshifts, total_counts, "ks-", ms=4, lw=1.2, label="Total")

    ax.set_xlabel("Redshift")
    ax.set_ylabel("N galaxies")
    ax.legend(fontsize=5.5, loc="upper right")
    ax.set_yscale("log")
    add_minor_ticks(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig04_sample_size.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig04_sample_size.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved phase1_fig04_sample_size.pdf")


# ═══════════════════════════════════════════════════════════════════════════
# Analysis 7+8: Alpha-element gradients + [O/Fe] vs [Fe/H]
# ═══════════════════════════════════════════════════════════════════════════

def analysis7_8_alpha_elements(df_z0, df_all):
    """Alpha-element gradient analysis + chemical evolution tracks."""
    print("\n" + "=" * 60)
    print("Analysis 7+8: Alpha-element gradients + [O/Fe] vs [Fe/H]")
    print("=" * 60)

    # --- Figure 5: Alpha-element gradient properties ---
    fig, axes = plt.subplots(2, 3, figsize=(7.014, 4.5))

    # Row 1: [O/Fe] and [Mg/Fe] gradient vs mass at z=0
    for col_idx, (scol, label) in enumerate([
        ("slope_OFe_star", r"$\nabla$[O/Fe]$_\star$"),
        ("slope_MgFe_star", r"$\nabla$[Mg/Fe]$_\star$"),
    ]):
        ax = axes[0, col_idx]
        logm = np.log10(df_z0["M_star_Msun"].values)
        y = df_z0[scol].values
        valid = np.isfinite(y)
        ax.scatter(logm[valid], y[valid], s=1.5, alpha=0.15,
                   color=WONG_COLORS[0], rasterized=True)
        c, med, lo, hi = running_median(logm, y, n_bins=7, min_count=8)
        good = np.isfinite(med)
        if good.any():
            ax.fill_between(c[good], lo[good], hi[good], alpha=0.2,
                            color=WONG_COLORS[1])
            ax.plot(c[good], med[good], "o-", color=WONG_COLORS[1], ms=3)
        ax.axhline(0, color="grey", ls="--", lw=0.6)
        ax.set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
        ax.set_ylabel(label + r" [dex/$R_e$]")
        ax.set_title(f"{label} vs mass (z=0)", fontsize=7)
        add_minor_ticks(ax)

    # Row 1, col 3: [O/Fe] gradient evolution with z
    ax = axes[0, 2]
    redshifts = sorted(df_all["redshift"].unique())
    for mlo, mhi, mlabel, mcolor in MASS_BINS:
        meds, zs = [], []
        for z in redshifts:
            mask = (np.isclose(df_all["redshift"], z, atol=0.01) &
                    (df_all["M_star_Msun"] >= mlo) &
                    (df_all["M_star_Msun"] < mhi))
            sub = df_all.loc[mask, "slope_OFe_star"]
            med, _, _ = bootstrap_median_ci(sub)
            if np.isfinite(med):
                meds.append(med); zs.append(z)
        if zs:
            ax.plot(zs, meds, "o-", color=mcolor, ms=3, label=mlabel)
    ax.axhline(0, color="grey", ls="--", lw=0.6)
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"$\nabla$[O/Fe]$_\star$ [dex/$R_e$]")
    ax.set_title("[O/Fe] gradient evolution", fontsize=7)
    ax.legend(fontsize=5, loc="best")
    add_minor_ticks(ax)

    # Row 2: Global [O/Fe] vs mass, [O/Fe] vs [Fe/H] at z=0, [Mg/Fe] vs [Fe/H]
    ax = axes[1, 0]
    logm = np.log10(df_z0["M_star_Msun"].values)
    ofe = df_z0["OFe_star_global"].values
    valid = np.isfinite(ofe)
    ax.scatter(logm[valid], ofe[valid], s=1.5, alpha=0.15,
               color=WONG_COLORS[2], rasterized=True)
    c, med, lo, hi = running_median(logm, ofe, n_bins=7, min_count=8)
    good = np.isfinite(med)
    if good.any():
        ax.plot(c[good], med[good], "o-", color=WONG_COLORS[1], ms=3)
    ax.set_xlabel(r"$\log_{10}(M_\star/M_\odot)$")
    ax.set_ylabel(r"[O/Fe]$_\star$ (global)")
    ax.set_title("Global [O/Fe] vs mass (z=0)", fontsize=7)
    add_minor_ticks(ax)

    # [O/Fe] vs [Fe/H] at z=0 colored by mass
    ax = axes[1, 1]
    feh = df_z0["FeH_star_global"].values
    ofe = df_z0["OFe_star_global"].values
    valid = np.isfinite(feh) & np.isfinite(ofe)
    sc = ax.scatter(feh[valid], ofe[valid], c=logm[valid], s=2, alpha=0.4,
                    cmap="viridis", vmin=9.5, vmax=12, rasterized=True)
    plt.colorbar(sc, ax=ax, label=r"$\log M_\star$", pad=0.02)
    ax.set_xlabel(r"[Fe/H]$_\star$ (global)")
    ax.set_ylabel(r"[O/Fe]$_\star$ (global)")
    ax.set_title("[O/Fe] vs [Fe/H] at z=0", fontsize=7)
    add_minor_ticks(ax)

    # [Mg/Fe] vs [Fe/H] at z=0 colored by mass
    ax = axes[1, 2]
    mgfe = df_z0["MgFe_star_global"].values
    valid_mg = np.isfinite(feh) & np.isfinite(mgfe)
    sc2 = ax.scatter(feh[valid_mg], mgfe[valid_mg], c=logm[valid_mg], s=2,
                     alpha=0.4, cmap="viridis", vmin=9.5, vmax=12,
                     rasterized=True)
    plt.colorbar(sc2, ax=ax, label=r"$\log M_\star$", pad=0.02)
    ax.set_xlabel(r"[Fe/H]$_\star$ (global)")
    ax.set_ylabel(r"[Mg/Fe]$_\star$ (global)")
    ax.set_title("[Mg/Fe] vs [Fe/H] at z=0", fontsize=7)
    add_minor_ticks(ax)

    fig.suptitle("Alpha-Element Abundance Properties", fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig05_alpha_gradients.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig05_alpha_gradients.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved phase1_fig05_alpha_gradients.pdf")

    # --- Figure 6: [O/Fe] vs [Fe/H] evolutionary tracks ---
    z_select = [z for z in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.01]
                if z in df_all["redshift"].unique()
                or any(np.isclose(df_all["redshift"].unique(), z, atol=0.05))]
    z_available = sorted(df_all["redshift"].unique())
    z_select_actual = []
    for zt in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.01]:
        diffs = np.abs(np.array(z_available) - zt)
        closest = z_available[np.argmin(diffs)]
        if diffs.min() < 0.15 and closest not in z_select_actual:
            z_select_actual.append(closest)

    n_panels = len(z_select_actual)
    ncols = min(4, n_panels)
    nrows = (n_panels + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(7.014, 2.8 * nrows),
                             squeeze=False)
    axes_flat = axes.flatten()

    for i, z in enumerate(z_select_actual):
        ax = axes_flat[i]
        sub = df_all[np.isclose(df_all["redshift"], z, atol=0.05)]
        feh_g = sub["FeH_star_global"].values
        ofe_g = sub["OFe_star_global"].values
        valid = np.isfinite(feh_g) & np.isfinite(ofe_g)

        if valid.sum() > 5:
            ax.scatter(feh_g[valid], ofe_g[valid], s=1, alpha=0.15,
                       color="grey", rasterized=True)
            for mlo, mhi, mlabel, mcolor in MASS_BINS:
                mm = valid & (sub["M_star_Msun"].values >= mlo) & \
                     (sub["M_star_Msun"].values < mhi)
                if mm.sum() >= 10:
                    feh_bins = np.linspace(np.percentile(feh_g[mm], 5),
                                           np.percentile(feh_g[mm], 95), 7)
                    c, med, lo, hi = running_median(feh_g[mm], ofe_g[mm],
                                                    n_bins=6, min_count=5)
                    good = np.isfinite(med)
                    if good.any():
                        ax.plot(c[good], med[good], "o-", ms=2.5, lw=1,
                                color=mcolor,
                                label=mlabel if i == 0 else None)

        ax.set_title(f"z={z:.1f} (N={len(sub)})", fontsize=7)
        ax.set_xlabel(r"[Fe/H]$_\star$")
        ax.set_ylabel(r"[O/Fe]$_\star$")
        add_minor_ticks(ax)

    # Hide unused panels
    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].set_visible(False)

    if n_panels > 0:
        axes_flat[0].legend(fontsize=5, loc="best")

    fig.suptitle("[O/Fe] vs [Fe/H] Across Redshift", fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig06_ofe_feh_tracks.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig06_ofe_feh_tracks.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved phase1_fig06_ofe_feh_tracks.pdf")


# ═══════════════════════════════════════════════════════════════════════════
# Analysis 9: Gas depletion, SFR profiles, quenching
# ═══════════════════════════════════════════════════════════════════════════

def analysis9_gas_quenching(df_z0, df_all):
    """Gas fraction evolution, SFR gradients, Z_gas vs Z_star comparison."""
    print("\n" + "=" * 60)
    print("Analysis 9: Gas depletion, SFR profiles, quenching")
    print("=" * 60)

    redshifts = sorted(df_all["redshift"].unique())

    fig, axes = plt.subplots(2, 3, figsize=(7.014, 4.5))

    # --- Panel (0,0): Fraction of gas-bearing ellipticals vs z ---
    ax = axes[0, 0]
    gas_fracs, gas_zs = [], []
    for z in redshifts:
        sub = df_all[np.isclose(df_all["redshift"], z, atol=0.01)]
        if len(sub) < 5:
            continue
        has_gas = sub["slope_Z_gas"].notna().sum()
        gas_fracs.append(has_gas / len(sub) * 100)
        gas_zs.append(z)

    ax.plot(gas_zs, gas_fracs, "o-", color=WONG_COLORS[0], ms=4)
    ax.set_xlabel("Redshift")
    ax.set_ylabel("Gas-bearing fraction [%]")
    ax.set_title("Galaxies with measurable gas", fontsize=7)
    ax.set_ylim(0, 105)
    add_minor_ticks(ax)

    # --- Panel (0,1): SFR gradient evolution ---
    ax = axes[0, 1]
    for mlo, mhi, mlabel, mcolor in MASS_BINS:
        meds, zs = [], []
        for z in redshifts:
            mask = (np.isclose(df_all["redshift"], z, atol=0.01) &
                    (df_all["M_star_Msun"] >= mlo) &
                    (df_all["M_star_Msun"] < mhi))
            sub = df_all.loc[mask, "slope_SFR_profile"]
            med, _, _ = bootstrap_median_ci(sub)
            if np.isfinite(med):
                meds.append(med); zs.append(z)
        if zs:
            ax.plot(zs, meds, "o-", color=mcolor, ms=3, label=mlabel)
    ax.axhline(0, color="grey", ls="--", lw=0.6)
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"$\nabla$SFR [dex/$R_e$]")
    ax.set_title("SFR gradient evolution", fontsize=7)
    ax.legend(fontsize=5, loc="best")
    add_minor_ticks(ax)

    # --- Panel (0,2): Z_gas gradient evolution ---
    ax = axes[0, 2]
    for mlo, mhi, mlabel, mcolor in MASS_BINS:
        meds, zs = [], []
        for z in redshifts:
            mask = (np.isclose(df_all["redshift"], z, atol=0.01) &
                    (df_all["M_star_Msun"] >= mlo) &
                    (df_all["M_star_Msun"] < mhi))
            sub = df_all.loc[mask, "slope_Z_gas"]
            med, _, _ = bootstrap_median_ci(sub)
            if np.isfinite(med):
                meds.append(med); zs.append(z)
        if zs:
            ax.plot(zs, meds, "o-", color=mcolor, ms=3, label=mlabel)
    ax.axhline(0, color="grey", ls="--", lw=0.6)
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"$\nabla Z_\mathrm{gas}$ [dex/$R_e$]")
    ax.set_title("Gas metallicity gradient evolution", fontsize=7)
    add_minor_ticks(ax)

    # --- Panel (1,0): Z_gas vs Z_star gradient at z=0 ---
    ax = axes[1, 0]
    x = df_z0["slope_Z_star"].values
    y = df_z0["slope_Z_gas"].values
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() > 5:
        ax.scatter(x[valid], y[valid], s=2, alpha=0.2, color=WONG_COLORS[0],
                   rasterized=True)
        lims = [min(x[valid].min(), y[valid].min()) - 0.1,
                max(x[valid].max(), y[valid].max()) + 0.1]
        ax.plot(lims, lims, "k--", lw=0.6, alpha=0.5)
        rho, _ = sp_stats.spearmanr(x[valid], y[valid])
        ax.text(0.05, 0.95, f"$\\rho$={rho:.2f}", transform=ax.transAxes,
                fontsize=6, va="top")
    ax.set_xlabel(r"$\nabla Z_\star$ [dex/$R_e$]")
    ax.set_ylabel(r"$\nabla Z_\mathrm{gas}$ [dex/$R_e$]")
    ax.set_title(r"$Z_\mathrm{gas}$ vs $Z_\star$ gradient (z=0)", fontsize=7)
    add_minor_ticks(ax)

    # --- Panel (1,1): Z_star and Z_gas gradient evolution (median) ---
    ax = axes[1, 1]
    for slope_col, label, color, ls in [
        ("slope_Z_star", r"$\nabla Z_\star$", WONG_COLORS[0], "-"),
        ("slope_Z_gas", r"$\nabla Z_\mathrm{gas}$", WONG_COLORS[1], "--"),
        ("slope_FeH_star", r"$\nabla$[Fe/H]$_\star$", WONG_COLORS[2], "-"),
    ]:
        meds, zs = [], []
        for z in redshifts:
            sub = df_all.loc[np.isclose(df_all["redshift"], z, atol=0.01),
                             slope_col]
            med, _, _ = bootstrap_median_ci(sub)
            if np.isfinite(med):
                meds.append(med); zs.append(z)
        if zs:
            ax.plot(zs, meds, "o" + ls, color=color, ms=3, label=label)
    ax.axhline(0, color="grey", ls="--", lw=0.6)
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"Median slope [dex/$R_e$]")
    ax.set_title("Stellar vs gas gradient evolution", fontsize=7)
    ax.legend(fontsize=5, loc="best")
    add_minor_ticks(ax)

    # --- Panel (1,2): Gradient slope vs sSFR at z=0 ---
    ax = axes[1, 2]
    if "sSFR" in df_z0.columns:
        ssfr = df_z0["sSFR"].values
        feh_slope = df_z0["slope_FeH_star"].values
        valid = np.isfinite(ssfr) & np.isfinite(feh_slope) & (ssfr > 0)
        if valid.sum() > 10:
            log_ssfr = np.log10(ssfr[valid])
            ax.scatter(log_ssfr, feh_slope[valid], s=1.5, alpha=0.15,
                       color=WONG_COLORS[0], rasterized=True)
            c, med, lo, hi = running_median(log_ssfr, feh_slope[valid],
                                            n_bins=6, min_count=5)
            good = np.isfinite(med)
            if good.any():
                ax.fill_between(c[good], lo[good], hi[good], alpha=0.2,
                                color=WONG_COLORS[1])
                ax.plot(c[good], med[good], "o-", color=WONG_COLORS[1], ms=3)
            rho, _ = sp_stats.spearmanr(log_ssfr, feh_slope[valid])
            ax.text(0.05, 0.95, f"$\\rho$={rho:.2f}", transform=ax.transAxes,
                    fontsize=6, va="top")
    ax.axhline(0, color="grey", ls="--", lw=0.6)
    ax.set_xlabel(r"$\log_{10}$ sSFR [yr$^{-1}$]")
    ax.set_ylabel(r"$\nabla$[Fe/H]$_\star$ [dex/$R_e$]")
    ax.set_title("[Fe/H] gradient vs sSFR (z=0)", fontsize=7)
    add_minor_ticks(ax)

    fig.suptitle("Gas Depletion, SFR, and Quenching", fontsize=10, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig07_gas_quenching.pdf"),
                bbox_inches="tight")
    fig.savefig(os.path.join(PLOTDIR, "phase1_fig07_gas_quenching.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved phase1_fig07_gas_quenching.pdf")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(PLOTDIR, exist_ok=True)
    os.makedirs(TABLEDIR, exist_ok=True)

    print("Rebuilding combined catalog...")
    df_all = rebuild_combined_catalog()

    df_z0 = df_all[np.isclose(df_all["redshift"], 0.0, atol=0.01)].copy()
    print(f"\nz=0 sample: {len(df_z0)} galaxies")
    print(f"Full sample: {len(df_all)} rows across "
          f"{df_all['redshift'].nunique()} redshifts")

    analysis1_z0_distributions(df_z0)
    analysis2_mass_dependence(df_z0, df_all)
    analysis4_redshift_evolution(df_all)
    analysis7_8_alpha_elements(df_z0, df_all)
    analysis9_gas_quenching(df_z0, df_all)

    print("\n" + "=" * 60)
    print("Phase 1 analysis complete!")
    print("=" * 60)
    print(f"\nPlots saved to: {PLOTDIR}/")
    print(f"Tables saved to: {TABLEDIR}/")

    # Summary of what needs re-processing for Phase 2
    old_snaps = [25, 33, 50, 67, 78, 84, 91]
    print(f"\nNote: Snapshots {old_snaps} have old-format catalogs (30 cols).")
    print("Re-run pipeline on these for inner/outer slopes, CNO, and "
          "environment info needed in Phase 2.")


if __name__ == "__main__":
    main()
