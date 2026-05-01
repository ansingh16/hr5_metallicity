#!/usr/bin/env python3
"""
Generate redshift evolution plots for metallicity gradients.

Figure 1: 2x3 gradient evolution by mass bin
Figure 2: 1x3 mass-gradient relation at different redshifts
Figure 3: 2x4 binned [O/Fe] vs [Fe/H] across redshift
Figure 4: 2x3 inner vs outer gradient evolution
Figure 5: 1x3 CNO gradient evolution
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

style_path = os.path.join(os.path.dirname(__file__), "..", "paper.mplstyle")
plt.style.use(style_path)

MASS_BINS = [
    (2e9, 1e10, r"$2{\times}10^9 < M_* < 10^{10}$", "C0"),
    (1e10, 5e10, r"$10^{10} < M_* < 5{\times}10^{10}$", "C1"),
    (5e10, 2e11, r"$5{\times}10^{10} < M_* < 2{\times}10^{11}$", "C2"),
    (2e11, 1e13, r"$M_* > 2{\times}10^{11}$", "C3"),
]


def bootstrap_median_ci(values, n_boot=1000, ci=68):
    values = values.dropna().values
    if len(values) < 5:
        return np.nan, np.nan, np.nan
    medians = np.array([
        np.median(np.random.choice(values, size=len(values), replace=True))
        for _ in range(n_boot)
    ])
    lo = np.percentile(medians, (100 - ci) / 2)
    hi = np.percentile(medians, 100 - (100 - ci) / 2)
    return np.median(values), lo, hi


def plot_gradient_evolution(df, slope_col, ylabel, title, ax, show_mass_bins=True):
    redshifts = sorted(df["redshift"].unique())
    if show_mass_bins:
        for mlo, mhi, label, color in MASS_BINS:
            meds, los, his, zs = [], [], [], []
            for z in redshifts:
                mask = (df["redshift"] == z) & (df["M_star_Msun"] >= mlo) & (df["M_star_Msun"] < mhi)
                sub = df.loc[mask, slope_col]
                if len(sub) < 5:
                    continue
                med, lo, hi = bootstrap_median_ci(sub)
                meds.append(med); los.append(lo); his.append(hi); zs.append(z)
            if zs:
                ax.fill_between(zs, los, his, alpha=0.15, color=color)
                ax.plot(zs, meds, "o-", color=color, ms=4, label=label)
    else:
        meds, los, his = [], [], []
        for z in redshifts:
            sub = df.loc[df["redshift"] == z, slope_col]
            med, lo, hi = bootstrap_median_ci(sub)
            meds.append(med); los.append(lo); his.append(hi)
        ax.fill_between(redshifts, los, his, alpha=0.2, color="k")
        ax.plot(redshifts, meds, "ko-", ms=5)
    ax.axhline(0, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("Redshift")
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def binned_median(x, y, bins, min_count=5):
    centers = 0.5 * (bins[:-1] + bins[1:])
    meds, lo, hi = [], [], []
    for i in range(len(bins) - 1):
        m = (x >= bins[i]) & (x < bins[i + 1]) & np.isfinite(y)
        n = m.sum()
        if n >= min_count:
            meds.append(np.median(y[m]))
            lo.append(np.percentile(y[m], 16))
            hi.append(np.percentile(y[m], 84))
        else:
            meds.append(np.nan); lo.append(np.nan); hi.append(np.nan)
    return centers, np.array(meds), np.array(lo), np.array(hi)


def main():
    csv_path = os.path.join(os.path.dirname(__file__), "..", cfg.CATALOG_DIR,
                            "gradients_evolution.csv")
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} rows, redshifts: {sorted(df['redshift'].unique())}")

    outdir = os.path.join(os.path.dirname(__file__), "..", cfg.PLOT_DIR)
    os.makedirs(outdir, exist_ok=True)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 1: 2x3 gradient evolution by mass bin
    # ══════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(2, 3, figsize=(7.014, 5.0))

    plot_gradient_evolution(df, "slope_FeH_star",
                           r"$\nabla$[Fe/H]$_*$ [dex/$R_e$]",
                           "Stellar [Fe/H] gradient", axes[0, 0])
    axes[0, 0].legend(loc="lower left", fontsize=5.5)

    plot_gradient_evolution(df, "slope_OFe_star",
                           r"$\nabla$[O/Fe]$_*$ [dex/$R_e$]",
                           "Stellar [O/Fe] gradient", axes[0, 1])

    plot_gradient_evolution(df, "slope_MgFe_star",
                           r"$\nabla$[Mg/Fe]$_*$ [dex/$R_e$]",
                           "Stellar [Mg/Fe] gradient", axes[0, 2])

    # Row 2: global properties + sample size
    for mlo, mhi, label, color in MASS_BINS:
        meds, zs = [], []
        for z in sorted(df["redshift"].unique()):
            mask = (df["redshift"] == z) & (df["M_star_Msun"] >= mlo) & (df["M_star_Msun"] < mhi)
            sub = df.loc[mask, "FeH_star_global"].dropna()
            if len(sub) >= 5:
                meds.append(np.median(sub)); zs.append(z)
        if zs:
            axes[1, 0].plot(zs, meds, "o-", color=color, ms=4, label=label)
    axes[1, 0].set_xlabel("Redshift")
    axes[1, 0].set_ylabel(r"[Fe/H]$_*$ (global)")
    axes[1, 0].set_title("Global [Fe/H] evolution")
    axes[1, 0].legend(loc="lower left", fontsize=5.5)

    for mlo, mhi, label, color in MASS_BINS:
        meds, zs = [], []
        for z in sorted(df["redshift"].unique()):
            mask = (df["redshift"] == z) & (df["M_star_Msun"] >= mlo) & (df["M_star_Msun"] < mhi)
            sub = df.loc[mask, "OFe_star_global"].dropna()
            if len(sub) >= 5:
                meds.append(np.median(sub)); zs.append(z)
        if zs:
            axes[1, 1].plot(zs, meds, "o-", color=color, ms=4)
    axes[1, 1].set_xlabel("Redshift")
    axes[1, 1].set_ylabel(r"[O/Fe]$_*$ (global)")
    axes[1, 1].set_title("Global [O/Fe] evolution")

    redshifts = sorted(df["redshift"].unique())
    for mlo, mhi, label, color in MASS_BINS:
        counts = [((df["redshift"] == z) & (df["M_star_Msun"] >= mlo) & (df["M_star_Msun"] < mhi)).sum()
                  for z in redshifts]
        axes[1, 2].plot(redshifts, counts, "o-", color=color, ms=4, label=label)
    axes[1, 2].set_xlabel("Redshift")
    axes[1, 2].set_ylabel("N galaxies")
    axes[1, 2].set_title("Sample size")
    axes[1, 2].legend(loc="upper right", fontsize=5.5)

    fig.suptitle("TNG-Cluster: Metallicity Gradient Evolution (z=0\u20133)", fontsize=12, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(outdir, "evolution_gradients.png"),
                dpi=300, bbox_inches="tight")
    print("Saved evolution_gradients.png")
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 2: mass-gradient at different redshifts
    # ══════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(7.014, 2.6))
    z_select = [0.0, 0.5, 1.0, 2.0, 3.0]
    cmap = plt.cm.coolwarm_r
    colors_z = {z: cmap(i / (len(z_select) - 1)) for i, z in enumerate(z_select)}

    for z in z_select:
        sub = df[df["redshift"] == z]
        if len(sub) < 10:
            continue
        logm = np.log10(sub["M_star_Msun"].values)
        for i, (col, ylabel) in enumerate([
            ("slope_FeH_star", r"$\nabla$[Fe/H]$_*$"),
            ("slope_OFe_star", r"$\nabla$[O/Fe]$_*$"),
            ("slope_MgFe_star", r"$\nabla$[Mg/Fe]$_*$"),
        ]):
            axes[i].scatter(logm, sub[col].values, s=3, alpha=0.1, color=colors_z[z])
            bins_m = np.linspace(max(logm.min(), 9.5), min(logm.max(), 12.5), 8)
            c, med, lo, hi = binned_median(logm, sub[col].values, bins_m)
            axes[i].plot(c, med, "o-", color=colors_z[z], ms=4,
                         label=f"z={z:.1f}" if i == 0 else f"z={z:.1f}")
            axes[i].set_xlabel(r"$\log_{10}(M_*/M_\odot)$")
            axes[i].set_ylabel(ylabel + r" [dex/$R_e$]")

    for ax in axes:
        ax.axhline(0, color="grey", ls="--", lw=0.8)
    axes[0].legend(fontsize=6)
    axes[0].set_title("[Fe/H] gradient vs mass")
    axes[1].set_title("[O/Fe] gradient vs mass")
    axes[2].set_title("[Mg/Fe] gradient vs mass")

    fig.suptitle("Mass\u2013Gradient Relation at Different Redshifts", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "evolution_mass_gradient.png"),
                dpi=300, bbox_inches="tight")
    print("Saved evolution_mass_gradient.png")
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 3: Binned [O/Fe] vs [Fe/H] across redshift (with median tracks)
    # ══════════════════════════════════════════════════════════════════════
    redshifts_sorted = sorted(df["redshift"].unique())
    n_z = len(redshifts_sorted)
    ncols = 4
    nrows = (n_z + ncols - 1) // ncols

    fig = plt.figure(figsize=(7.014, 3.2 * nrows))
    gs = fig.add_gridspec(nrows, ncols + 1, width_ratios=[1]*ncols + [0.05],
                          hspace=0.25, wspace=0.08)

    sc = None
    axes_ofe = []
    for i, z in enumerate(redshifts_sorted):
        row, col = divmod(i, ncols)
        share_x = axes_ofe[0] if axes_ofe else None
        share_y = axes_ofe[0] if axes_ofe else None
        ax = fig.add_subplot(gs[row, col], sharex=share_x, sharey=share_y)
        axes_ofe.append(ax)
        sub = df[df["redshift"] == z]
        feh_g = sub["FeH_star_global"].values
        ofe_g = sub["OFe_star_global"].values
        logm = np.log10(sub["M_star_Msun"].clip(lower=1).values)

        # Hexbin background
        valid = np.isfinite(feh_g) & np.isfinite(ofe_g)
        if valid.sum() > 10:
            ax.hexbin(feh_g[valid], ofe_g[valid], gridsize=20, cmap="Greys",
                      mincnt=1, linewidths=0.2, alpha=0.6)
            # Median tracks by mass bin
            for j, (mlo, mhi, label, color) in enumerate(MASS_BINS):
                mm = valid & (sub["M_star_Msun"].values >= mlo) & (sub["M_star_Msun"].values < mhi)
                if mm.sum() < 10:
                    continue
                feh_bins = np.linspace(np.percentile(feh_g[mm], 5),
                                       np.percentile(feh_g[mm], 95), 7)
                c, med, lo, hi = binned_median(feh_g[mm], ofe_g[mm], feh_bins)
                ax.plot(c, med, "o-", ms=2.5, lw=1, color=color,
                        label=label if i == 0 else None)

        ax.set_title(f"z={z:.1f} (N={len(sub)})", fontsize=8)
        if col == 0:
            ax.set_ylabel(r"[O/Fe]$_*$")
        else:
            ax.tick_params(labelleft=False)
        if row == nrows - 1:
            ax.set_xlabel(r"[Fe/H]$_*$")
        else:
            ax.tick_params(labelbottom=False)

    axes_ofe[0].set_xlim(-1.4, 0.35)
    axes_ofe[0].set_ylim(0.3, 0.75)

    for i in range(n_z, nrows * ncols):
        row, col = divmod(i, ncols)
        ax = fig.add_subplot(gs[row, col])
        ax.set_visible(False)

    # Legend in first panel
    axes_ofe[0].legend(fontsize=4.5, loc="lower left")

    fig.suptitle("[O/Fe] vs [Fe/H] Across Redshift (binned medians by mass)",
                 fontsize=10)
    fig.savefig(os.path.join(outdir, "evolution_ofe_feh.png"),
                dpi=300, bbox_inches="tight")
    print("Saved evolution_ofe_feh.png")
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 4: Inner vs outer gradient evolution (2x3)
    # ══════════════════════════════════════════════════════════════════════
    has_inner = "slope_inner_FeH_star" in df.columns
    if has_inner:
        fig, axes = plt.subplots(2, 3, figsize=(7.014, 5.0))

        for col, (key, label) in enumerate([
            ("FeH_star", r"[Fe/H]$_*$"),
            ("OFe_star", r"[O/Fe]$_*$"),
            ("MgFe_star", r"[Mg/Fe]$_*$"),
        ]):
            # Row 1: inner gradient evolution
            plot_gradient_evolution(df, f"slope_inner_{key}",
                                   f"Inner $\\nabla${label}",
                                   f"Inner {label} ($<R_e$)",
                                   axes[0, col])
            if col == 0:
                axes[0, col].legend(loc="lower left", fontsize=5.5)

            # Row 2: outer gradient evolution
            plot_gradient_evolution(df, f"slope_outer_{key}",
                                   f"Outer $\\nabla${label}",
                                   f"Outer {label} ($>R_e$)",
                                   axes[1, col])

        fig.suptitle("Inner vs Outer Gradient Evolution", fontsize=12, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(os.path.join(outdir, "evolution_inner_outer.png"),
                    dpi=300, bbox_inches="tight")
        print("Saved evolution_inner_outer.png")
        plt.close(fig)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 5: CNO gradient evolution (1x3)
    # ══════════════════════════════════════════════════════════════════════
    has_cno = "slope_NO_star" in df.columns
    if has_cno:
        fig, axes = plt.subplots(1, 3, figsize=(7.014, 2.6))
        for col, (key, label) in enumerate([
            ("NO_star", r"[N/O]$_*$"),
            ("CO_star", r"[C/O]$_*$"),
            ("CN_star", r"[C/N]$_*$"),
        ]):
            plot_gradient_evolution(df, f"slope_{key}",
                                   f"$\\nabla${label} [dex/$R_e$]",
                                   f"{label} gradient evolution",
                                   axes[col])
            if col == 0:
                axes[col].legend(loc="best", fontsize=5.5)

        fig.suptitle("CNO Gradient Evolution", fontsize=11, y=1.02)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "evolution_cno.png"),
                    dpi=300, bbox_inches="tight")
        print("Saved evolution_cno.png")
        plt.close(fig)

    print("\nDone! All evolution plots saved.")


if __name__ == "__main__":
    main()
