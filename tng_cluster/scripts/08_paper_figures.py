#!/usr/bin/env python3
"""
Paper figures for TNG-Cluster metallicity gradient study.

Produces 5 focused, story-driven figures (no grids) in plots/paper/.
Each figure communicates one clear result:

  Fig 1: z=0 gradient distributions ([Fe/H] + [O/Fe])
  Fig 2: Mass dependence of [Fe/H] gradient at z=0
  Fig 3: Gradient evolution with redshift by mass bin
  Fig 4: [O/Fe] vs [Fe/H] chemical enrichment diagnostic
  Fig 5: Gas depletion and quenching connection
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
# Figure 1: z=0 gradient distributions
# ═══════════════════════════════════════════════════════════════════════════

def fig1_z0_distributions(df_z0):
    """[Fe/H] and [O/Fe] gradient histograms — the headline result."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.014, 2.8))

    for ax, col, label, color in [
        (ax1, "slope_FeH_star", r"$\nabla[\mathrm{Fe/H}]_\star$", WONG[0]),
        (ax2, "slope_OFe_star", r"$\nabla[\mathrm{O/Fe}]_\star$", WONG[2]),
    ]:
        vals = df_z0[col].dropna().values
        med = np.median(vals)
        p16, p84 = np.percentile(vals, [16, 84])
        frac = np.mean(vals < 0) * 100 if "FeH" in col else np.mean(vals > 0) * 100
        sign_label = "negative" if "FeH" in col else "positive"

        ax.hist(vals, bins=30, color=color, alpha=0.75, edgecolor="white",
                linewidth=0.3)
        ax.axvline(med, color="k", ls="-", lw=1.0, zorder=5)
        ax.axvline(0, color="grey", ls="--", lw=0.6)
        ax.axvspan(p16, p84, alpha=0.08, color="k")

        ax.text(0.97, 0.95,
                f"med = {med:+.3f}\n"
                f"[{p16:+.3f}, {p84:+.3f}]\n"
                f"{frac:.1f}% {sign_label}",
                transform=ax.transAxes, va="top", ha="right",
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec="none", alpha=0.85))

        ax.set_xlabel(r"Slope [dex/$R_e$]")
        ax.set_ylabel("Count")
        ax.set_title(label)
        add_minor_ticks(ax)

    fig.text(0.5, -0.02, f"N = {len(df_z0)} ellipticals at z = 0",
             ha="center", style="italic")
    fig.tight_layout()
    savefig(fig, "fig1_z0_distributions")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2: Mass dependence
# ═══════════════════════════════════════════════════════════════════════════

def fig2_mass_dependence(df_z0):
    """[Fe/H] gradient vs stellar mass at z=0."""
    fig, ax = plt.subplots(figsize=(3.321, 2.5))

    logm = np.log10(df_z0["M_star_Msun"].values)
    y = df_z0["slope_FeH_star"].values
    valid = np.isfinite(y) & np.isfinite(logm)

    for mlo, mhi, mlabel, mcolor in MASS_BINS:
        lm = np.log10(mlo)
        lh = np.log10(mhi)
        mm = valid & (logm >= lm) & (logm < lh)
        ax.scatter(logm[mm], y[mm], s=2, alpha=0.2, color=mcolor,
                   rasterized=True, label=mlabel)

    c, med, lo, hi = running_median(logm[valid], y[valid], n_bins=8,
                                     min_count=10)
    good = np.isfinite(med)
    if good.any():
        ax.fill_between(c[good], lo[good], hi[good], alpha=0.15, color="k")
        ax.plot(c[good], med[good], "o-", color="k", ms=3, lw=1.2,
                zorder=5, label="Running median")

    ax.axhline(0, color="grey", ls="--", lw=0.6)

    ax.set_xlabel(r"$\log_{10}(M_\star\,/\,M_\odot)$")
    ax.set_ylabel(r"$\nabla[\mathrm{Fe/H}]_\star$ [dex/$R_e$]")
    ax.legend(loc="upper right", ncol=2, handletextpad=0.3,
              columnspacing=0.8)
    add_minor_ticks(ax)
    fig.tight_layout()
    savefig(fig, "fig2_mass_dependence")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3: Gradient evolution with redshift
# ═══════════════════════════════════════════════════════════════════════════

def fig3_redshift_evolution(df_all):
    """[Fe/H] and [O/Fe] gradient evolution by mass bin — the central result."""
    fig, axes = plt.subplots(2, 2, figsize=(7.014, 4.5),
                             height_ratios=[3, 1], sharex=True)

    redshifts = sorted(df_all["redshift"].unique())
    z_max = max(redshifts) + 0.1

    for col_idx, (slope_col, ylabel) in enumerate([
        ("slope_FeH_star", r"$\nabla[\mathrm{Fe/H}]_\star$"),
        ("slope_OFe_star", r"$\nabla[\mathrm{O/Fe}]_\star$"),
    ]):
        ax = axes[0, col_idx]
        ax_n = axes[1, col_idx]

        # All galaxies (grey background)
        meds_all, lo_all, hi_all, zs_all = [], [], [], []
        for z in redshifts:
            if z > z_max:
                continue
            sub = df_all.loc[np.isclose(df_all["redshift"], z, atol=0.01),
                             slope_col]
            med, lo, hi = bootstrap_median_ci(sub)
            if np.isfinite(med):
                meds_all.append(med)
                lo_all.append(lo)
                hi_all.append(hi)
                zs_all.append(z)
        if zs_all:
            ax.fill_between(zs_all, lo_all, hi_all, alpha=0.1, color="grey")
            ax.plot(zs_all, meds_all, "s-", color="grey", ms=2.5, lw=0.7,
                    alpha=0.5, label="All", zorder=1)

        for mlo, mhi, mlabel, mcolor in MASS_BINS:
            meds, los, his, zs = [], [], [], []
            counts_z, counts_n = [], []
            for z in redshifts:
                if z > z_max:
                    continue
                mask = (np.isclose(df_all["redshift"], z, atol=0.01) &
                        (df_all["M_star_Msun"] >= mlo) &
                        (df_all["M_star_Msun"] < mhi))
                n = mask.sum()
                counts_z.append(z)
                counts_n.append(n)
                sub = df_all.loc[mask, slope_col]
                med, lo, hi = bootstrap_median_ci(sub)
                if np.isfinite(med):
                    meds.append(med); los.append(lo); his.append(hi)
                    zs.append(z)
            if zs:
                ax.fill_between(zs, los, his, alpha=0.1, color=mcolor)
                ax.plot(zs, meds, "o-", color=mcolor, ms=3, label=mlabel,
                        zorder=3)
            ax_n.plot(counts_z, counts_n, "o-", color=mcolor, ms=2, lw=0.8)

        ax.axhline(0, color="grey", ls="--", lw=0.6)
        ax.set_ylabel(ylabel + r" [dex/$R_e$]")
        ax.set_xlim(-0.3, z_max)
        add_minor_ticks(ax)

        # Sample size panel
        total_c = [np.isclose(df_all["redshift"], z, atol=0.01).sum()
                    for z in redshifts if z <= z_max]
        total_z = [z for z in redshifts if z <= z_max]
        ax_n.plot(total_z, total_c, "s-", color="grey", ms=2, lw=0.8)
        ax_n.set_yscale("log")
        ax_n.set_ylabel("N")
        ax_n.set_xlabel("Redshift")
        ax_n.set_xlim(-0.3, z_max)
        ax_n.set_ylim(bottom=0.8)

    axes[0, 0].set_ylim(-0.38, 0.05)
    axes[0, 0].legend(loc="upper right", ncol=1)
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.08)
    savefig(fig, "fig3_gradient_evolution")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 4: [O/Fe] vs [Fe/H] chemical diagnostic
# ═══════════════════════════════════════════════════════════════════════════

def fig4_ofe_feh_diagnostic(df_z0, df_all):
    """[O/Fe] vs [Fe/H]: (a) mass-binned at z=0, (b) redshift evolution."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.014, 3.0),
                                    sharey=True)

    # ── Panel (a): mass-binned medians at z=0 ──
    logm = np.log10(df_z0["M_star_Msun"].values)
    feh = df_z0["FeH_star_global"].values
    ofe = df_z0["OFe_star_global"].values
    valid = np.isfinite(feh) & np.isfinite(ofe) & np.isfinite(logm)

    for mlo, mhi, mlabel, mcolor in MASS_BINS:
        lm, lh = np.log10(mlo), np.log10(mhi)
        mm = valid & (logm >= lm) & (logm < lh)
        if mm.sum() >= 15:
            c, med, lo, hi = running_median(feh[mm], ofe[mm], n_bins=6,
                                             min_count=5)
            good = np.isfinite(med)
            if good.any():
                ax1.fill_between(c[good], lo[good], hi[good],
                                 color=mcolor, alpha=0.18)
                ax1.plot(c[good], med[good], "o-", ms=3, lw=1.4,
                         color=mcolor, label=mlabel, zorder=4)

    ax1.set_xlabel(r"$[\mathrm{Fe/H}]_\star$ (global)")
    ax1.set_ylabel(r"$[\mathrm{O/Fe}]_\star$ (global)")
    ax1.legend(loc="upper right", handletextpad=0.3,
               title=r"$M_\star\;[M_\odot]$")
    ax1.text(0.05, 0.05, "(a) $z=0$ mass bins",
             transform=ax1.transAxes, va="bottom")
    add_minor_ticks(ax1)

    # ── Panel (b): redshift evolution track ──
    target_z = sorted(df_all["redshift"].unique())
    cmap = plt.cm.plasma_r
    norm = plt.Normalize(vmin=0, vmax=max(target_z) + 0.1)

    med_fehs, med_ofes, plot_zs = [], [], []
    for z in target_z:
        sub = df_all[np.isclose(df_all["redshift"], z, atol=0.05)]
        feh_z = sub["FeH_star_global"].dropna()
        ofe_z = sub["OFe_star_global"].dropna()
        common = sub.loc[feh_z.index.intersection(ofe_z.index)]
        if len(common) < 5:
            continue
        med_fehs.append(common["FeH_star_global"].median())
        med_ofes.append(common["OFe_star_global"].median())
        plot_zs.append(z)

    med_fehs = np.array(med_fehs)
    med_ofes = np.array(med_ofes)
    plot_zs = np.array(plot_zs)

    ax2.plot(med_fehs, med_ofes, "-", color="0.5", lw=1.0, zorder=2)
    for i, z in enumerate(plot_zs):
        ax2.scatter(med_fehs[i], med_ofes[i], s=55, c=[cmap(norm(z))],
                    marker="*", edgecolors="k", linewidths=0.4, zorder=5,
                    label=f"$z={z:g}$")

    ax2.set_xlabel(r"$[\mathrm{Fe/H}]_\star$ (global)")
    ax2.legend(loc="upper right", handletextpad=0.3,
               title="Redshift", ncol=2)
    ax2.text(0.05, 0.05, "(b) Redshift evolution",
             transform=ax2.transAxes, va="bottom")
    add_minor_ticks(ax2)

    fig.tight_layout()
    savefig(fig, "fig4_ofe_feh_diagnostic")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 5: Gas depletion and quenching
# ═══════════════════════════════════════════════════════════════════════════

def fig5_gas_quenching(df_z0, df_all):
    """Gas-bearing fraction vs z + stellar vs gas gradient evolution."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.014, 2.8))

    redshifts = sorted(df_all["redshift"].unique())

    # Left: gas-bearing fraction
    gas_fracs, gas_zs = [], []
    for z in redshifts:
        sub = df_all[np.isclose(df_all["redshift"], z, atol=0.01)]
        if len(sub) < 5:
            continue
        has_gas = sub["slope_Z_gas"].notna().sum()
        gas_fracs.append(has_gas / len(sub) * 100)
        gas_zs.append(z)

    ax1.plot(gas_zs, gas_fracs, "o-", color=WONG[0], ms=4)
    ax1.fill_between(gas_zs, 0, gas_fracs, alpha=0.1, color=WONG[0])
    ax1.set_xlabel("Redshift")
    ax1.set_ylabel("Gas-bearing fraction [%]")
    ax1.set_ylim(0, 105)
    ax1.set_xlim(-0.3, max(redshifts) + 0.5)

    # Annotate key transitions
    for z_ann, frac_ann in zip(gas_zs, gas_fracs):
        if z_ann in (0.0, 0.3, 0.7):
            ax1.annotate(f"{frac_ann:.0f}%", (z_ann, frac_ann),
                         textcoords="offset points", xytext=(6, 4),
                         color=WONG[0])

    add_minor_ticks(ax1)

    # Right: stellar vs gas gradient evolution (medians)
    tracks = [
        ("slope_FeH_star", r"$\nabla[\mathrm{Fe/H}]_\star$", WONG[0], "-"),
        ("slope_Z_star",   r"$\nabla Z_\star$",               WONG[2], "-"),
        ("slope_Z_gas",    r"$\nabla Z_\mathrm{gas}$",        WONG[1], "--"),
    ]
    for slope_col, label, color, ls in tracks:
        meds, los, his, zs = [], [], [], []
        for z in redshifts:
            sub = df_all.loc[np.isclose(df_all["redshift"], z, atol=0.01),
                             slope_col]
            med, lo, hi = bootstrap_median_ci(sub)
            if np.isfinite(med):
                meds.append(med); los.append(lo); his.append(hi)
                zs.append(z)
        if zs:
            ax2.fill_between(zs, los, his, alpha=0.08, color=color)
            ax2.plot(zs, meds, "o" + ls, color=color, ms=3, label=label)

    ax2.axhline(0, color="grey", ls="--", lw=0.6)
    ax2.set_xlabel("Redshift")
    ax2.set_ylabel(r"Median slope [dex/$R_e$]")
    ax2.set_xlim(-0.3, max(redshifts) + 0.5)
    ax2.legend(loc="lower left")
    add_minor_ticks(ax2)

    fig.tight_layout()
    savefig(fig, "fig5_gas_quenching")


# ═══════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(PLOTDIR, exist_ok=True)
    df_all, df_z0 = load_data()

    fig1_z0_distributions(df_z0)
    fig2_mass_dependence(df_z0)
    fig3_redshift_evolution(df_all)
    fig4_ofe_feh_diagnostic(df_z0, df_all)
    fig5_gas_quenching(df_z0, df_all)

    print(f"\nAll paper figures saved to {PLOTDIR}/")


if __name__ == "__main__":
    main()
