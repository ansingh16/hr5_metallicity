"""
Diagnostic plots for z=0 (snap 99) metallicity profiles and gradients.

Figure 1 (diagnostic_z0.png): 2x3 panel
  1. Stacked [Fe/H] stellar profiles in mass bins
  2. Stacked [O/Fe] stellar profiles in mass bins
  3. Stacked [Mg/Fe] stellar profiles in mass bins
  4. [Fe/H] gradient slope vs stellar mass (+ sigma_v)
  5. Binned [O/Fe] vs [Fe/H] diagram
  6. Inner vs outer gradient comparison

Figure 2 (diagnostic_cno_z0.png): 1x3 panel
  1-3. N/O, C/O, C/N radial profiles by mass bin

Figure 3 (diagnostic_conditional_z0.png): 2x3 panel
  grad(Fe/H) vs mass binned by [O/Fe], age, environment
"""

import os
import sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

style_path = os.path.join(os.path.dirname(__file__), "..", "paper.mplstyle")
plt.style.use(style_path)

MASS_BINS = [
    (2e9,   1e10,  r"$2{\times}10^9 < M_* < 10^{10}$"),
    (1e10,  5e10,  r"$10^{10} < M_* < 5{\times}10^{10}$"),
    (5e10,  1e11,  r"$5{\times}10^{10} < M_* < 10^{11}$"),
    (1e11,  1e13,  r"$M_* > 10^{11}$"),
]


def load_data(snap=99):
    basedir = os.path.join(os.path.dirname(__file__), "..")
    catdir = os.path.join(basedir, cfg.CATALOG_DIR)
    grad = pd.read_csv(os.path.join(catdir, f"gradients_snap{snap:03d}.csv"))
    prof_file = os.path.join(catdir, f"profiles_snap{snap:03d}.hdf5")
    profiles = h5py.File(prof_file, "r")
    return grad, profiles


def median_profile_in_mass_bin(profiles, grad, mass_lo, mass_hi, field):
    mask = (grad["M_star_Msun"] >= mass_lo) & (grad["M_star_Msun"] < mass_hi)
    sub_ids = grad.loc[mask, "subhalo_id"].values
    all_prof = []
    for sid in sub_ids:
        key = str(sid)
        if key not in profiles or field not in profiles[key]:
            continue
        prof = profiles[key][field][:]
        if not np.all(np.isnan(prof)):
            all_prof.append(prof)
    if len(all_prof) < 3:
        return None, None, None, 0
    arr = np.array(all_prof)
    return np.nanmedian(arr, axis=0), np.nanpercentile(arr, 16, axis=0), \
           np.nanpercentile(arr, 84, axis=0), len(all_prof)


def binned_median(x, y, bins, min_count=5):
    """Compute median and 16/84 percentiles in bins."""
    centers = 0.5 * (bins[:-1] + bins[1:])
    meds, lo, hi, ns = [], [], [], []
    for i in range(len(bins) - 1):
        m = (x >= bins[i]) & (x < bins[i + 1]) & np.isfinite(y)
        n = m.sum()
        ns.append(n)
        if n >= min_count:
            meds.append(np.median(y[m]))
            lo.append(np.percentile(y[m], 16))
            hi.append(np.percentile(y[m], 84))
        else:
            meds.append(np.nan)
            lo.append(np.nan)
            hi.append(np.nan)
    return centers, np.array(meds), np.array(lo), np.array(hi), np.array(ns)


def main():
    grad, profiles = load_data()

    first_key = list(profiles.keys())[0]
    r_bins = profiles[first_key]["r_bins"][:]

    outdir = os.path.join(os.path.dirname(__file__), "..", cfg.PLOT_DIR)
    os.makedirs(outdir, exist_ok=True)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 1: Main diagnostic (2x3)
    # ══════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(2, 3, figsize=(7.014, 5.0))

    # Row 1: Profiles by mass bin
    for col, (field, ylabel, title) in enumerate([
        ("FeH_star", r"[Fe/H]$_*$", "Stellar [Fe/H]"),
        ("OFe_star", r"[O/Fe]$_*$", "Stellar [O/Fe]"),
        ("MgFe_star", r"[Mg/Fe]$_*$", "Stellar [Mg/Fe]"),
    ]):
        ax = axes[0, col]
        for (mlo, mhi, label) in MASS_BINS:
            med, p16, p84, n = median_profile_in_mass_bin(
                profiles, grad, mlo, mhi, field)
            if med is None:
                continue
            ax.plot(r_bins, med, lw=1.5, label=f"{label} (N={n})")
            ax.fill_between(r_bins, p16, p84, alpha=0.15)
        ax.set_xlabel(r"$R / R_{\rm half}$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(0, 3)
        ax.axhline(0, color="gray", ls="--", lw=0.5)
        if col == 0:
            ax.legend(fontsize=5, loc="lower left")

    # Panel 4: gradient vs mass AND vs sigma_v
    ax = axes[1, 0]
    valid = grad["slope_FeH_star"].notna()
    logm = np.log10(grad.loc[valid, "M_star_Msun"])
    slopes = grad.loc[valid, "slope_FeH_star"].values

    ax.scatter(logm, slopes, s=6, alpha=0.3, color="C0", edgecolors="none")
    bins_m = np.arange(9.3, 12.5, 0.3)
    c, med, lo, hi, _ = binned_median(logm.values, slopes, bins_m)
    ax.plot(c, med, "ko-", lw=1.5, ms=4, zorder=5, label=r"Median vs $M_*$")
    ax.fill_between(c, lo, hi, alpha=0.15, color="k")

    # Overlay sigma_v on twin axis if available
    if "vel_disp_kms" in grad.columns:
        ax2 = ax.twiny()
        logsig = np.log10(grad.loc[valid, "vel_disp_kms"].clip(lower=1))
        bins_sig = np.arange(1.4, 2.8, 0.15)
        c_s, med_s, lo_s, hi_s, _ = binned_median(logsig.values, slopes, bins_sig)
        ax2.plot(c_s, med_s, "s--", color="C3", lw=1.2, ms=3, zorder=4,
                 label=r"Median vs $\sigma_v$")
        ax2.set_xlabel(r"$\log_{10}(\sigma_v\,/\,{\rm km\,s^{-1}})$",
                       fontsize=8, color="C3")
        ax2.tick_params(axis="x", labelcolor="C3", labelsize=7)

    ax.axhline(0, color="gray", ls="--", lw=0.5)
    ax.set_xlabel(r"$\log_{10}(M_*/M_\odot)$")
    ax.set_ylabel(r"$\nabla$[Fe/H]$_*$ [dex/$R_e$]")
    ax.set_title(r"$\nabla$[Fe/H] vs $M_*$ and $\sigma_v$")
    ax.legend(fontsize=6, loc="lower right")

    # Panel 5: Binned [O/Fe] vs [Fe/H]
    ax = axes[1, 1]
    has_both = grad["FeH_star_global"].notna() & grad["OFe_star_global"].notna()
    feh_g = grad.loc[has_both, "FeH_star_global"].values
    ofe_g = grad.loc[has_both, "OFe_star_global"].values
    logm_g = np.log10(grad.loc[has_both, "M_star_Msun"].values)

    # Hexbin background
    hb = ax.hexbin(feh_g, ofe_g, gridsize=25, cmap="Greys", mincnt=1,
                   linewidths=0.2)

    # Median tracks by mass bin
    for i, (mlo, mhi, label) in enumerate(MASS_BINS):
        m = (grad.loc[has_both, "M_star_Msun"].values >= mlo) & \
            (grad.loc[has_both, "M_star_Msun"].values < mhi)
        if m.sum() < 10:
            continue
        feh_bins = np.linspace(np.percentile(feh_g[m], 5),
                               np.percentile(feh_g[m], 95), 8)
        c, med, lo, hi, _ = binned_median(feh_g[m], ofe_g[m], feh_bins)
        ax.plot(c, med, "o-", ms=3, lw=1.2, label=label)

    ax.set_xlabel(r"[Fe/H]$_*$ (global)")
    ax.set_ylabel(r"[O/Fe]$_*$ (global)")
    ax.set_title("[O/Fe] vs [Fe/H] (binned)")
    ax.legend(fontsize=5, loc="upper right")

    # Panel 6: Inner vs outer gradient
    ax = axes[1, 2]
    has_io = grad["slope_inner_FeH_star"].notna() & grad["slope_outer_FeH_star"].notna()
    if has_io.sum() > 10:
        inner = grad.loc[has_io, "slope_inner_FeH_star"].values
        outer = grad.loc[has_io, "slope_outer_FeH_star"].values
        ax.scatter(inner, outer, s=6, alpha=0.3, c="C0", edgecolors="none")
        lim = max(abs(np.nanpercentile(inner, [2, 98])).max(),
                  abs(np.nanpercentile(outer, [2, 98])).max()) * 1.1
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, alpha=0.5)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel(r"Inner $\nabla$[Fe/H]$_*$ ($<R_e$)")
        ax.set_ylabel(r"Outer $\nabla$[Fe/H]$_*$ ($>R_e$)")
        ax.set_title("Inner vs outer gradient")
        ax.axhline(0, color="gray", ls="--", lw=0.5)
        ax.axvline(0, color="gray", ls="--", lw=0.5)
        # Annotate medians
        ax.text(0.05, 0.95, f"med inner={np.nanmedian(inner):.3f}\n"
                f"med outer={np.nanmedian(outer):.3f}",
                transform=ax.transAxes, fontsize=6, va="top")

    fig.suptitle(f"TNG-Cluster Ellipticals at z=0 (N={len(grad)})",
                 fontsize=12, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    path1 = os.path.join(outdir, "diagnostic_z0.png")
    fig.savefig(path1, dpi=300, bbox_inches="tight")
    print(f"Saved {path1}")
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 2: CNO profiles (1x3)
    # ══════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(7.014, 2.6))

    for col, (field, ylabel, title) in enumerate([
        ("NO_star", r"[N/O]$_*$", "Stellar [N/O]"),
        ("CO_star", r"[C/O]$_*$", "Stellar [C/O]"),
        ("CN_star", r"[C/N]$_*$", "Stellar [C/N]"),
    ]):
        ax = axes[col]
        for (mlo, mhi, label) in MASS_BINS:
            med, p16, p84, n = median_profile_in_mass_bin(
                profiles, grad, mlo, mhi, field)
            if med is None:
                continue
            ax.plot(r_bins, med, lw=1.5, label=f"{label} (N={n})")
            ax.fill_between(r_bins, p16, p84, alpha=0.15)
        ax.set_xlabel(r"$R / R_{\rm half}$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(0, 3)
        ax.axhline(0, color="gray", ls="--", lw=0.5)
        if col == 0:
            ax.legend(fontsize=5, loc="best")

    fig.suptitle("TNG-Cluster: CNO Abundance Profiles at z=0", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    path2 = os.path.join(outdir, "diagnostic_cno_z0.png")
    fig.savefig(path2, dpi=300, bbox_inches="tight")
    print(f"Saved {path2}")
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════════════
    # Figure 3: Conditional gradient-mass relation (2x3)
    # ══════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(2, 3, figsize=(7.014, 5.0))

    logm_all = np.log10(grad["M_star_Msun"].values)
    slope_feh = grad["slope_FeH_star"].values
    bins_m = np.arange(9.5, 12.5, 0.4)

    def plot_conditional(ax, cond_col, cond_label, n_bins=3):
        """Split sample by percentiles of cond_col, plot grad vs mass."""
        valid_mask = np.isfinite(slope_feh) & np.isfinite(grad[cond_col].values)
        vals = grad.loc[valid_mask, cond_col].values
        lm = logm_all[valid_mask]
        sl = slope_feh[valid_mask]

        percs = np.linspace(0, 100, n_bins + 1)
        boundaries = np.percentile(vals, percs)

        for j in range(n_bins):
            m = (vals >= boundaries[j]) & (vals < boundaries[j + 1])
            if m.sum() < 10:
                continue
            lo_val = boundaries[j]
            hi_val = boundaries[j + 1]
            c, med, plo, phi, _ = binned_median(lm[m], sl[m], bins_m)
            ax.plot(c, med, "o-", ms=3, lw=1.2,
                    label=f"{lo_val:.2f}–{hi_val:.2f}")
            ax.fill_between(c, plo, phi, alpha=0.1)

        ax.axhline(0, color="gray", ls="--", lw=0.5)
        ax.set_xlabel(r"$\log_{10}(M_*/M_\odot)$")
        ax.set_ylabel(r"$\nabla$[Fe/H]$_*$")
        ax.set_title(f"Binned by {cond_label}")
        ax.legend(fontsize=5.5, title=cond_label, title_fontsize=6)

    # Row 1: [O/Fe], [Mg/Fe], mean age
    if "OFe_star_global" in grad.columns:
        plot_conditional(axes[0, 0], "OFe_star_global", r"[O/Fe]$_*$")
    if "MgFe_star_global" in grad.columns:
        plot_conditional(axes[0, 1], "MgFe_star_global", r"[Mg/Fe]$_*$")
    if "age_star_global" in grad.columns:
        plot_conditional(axes[0, 2], "age_star_global", "Age [Gyr]")
    else:
        axes[0, 2].text(0.5, 0.5, "Age not\navailable",
                        transform=axes[0, 2].transAxes, ha="center")

    # Row 2: R/R200, cluster M200, vel_disp
    if "R_over_R200" in grad.columns and grad["R_over_R200"].notna().sum() > 50:
        plot_conditional(axes[1, 0], "R_over_R200", r"$R/R_{200}$")
    else:
        axes[1, 0].text(0.5, 0.5, "R/R200 not\navailable",
                        transform=axes[1, 0].transAxes, ha="center")
        axes[1, 0].set_title(r"Binned by $R/R_{200}$")

    if "cluster_M200_Msun" in grad.columns and grad["cluster_M200_Msun"].notna().sum() > 50:
        # Use log M200 for binning
        grad["log_cluster_M200"] = np.log10(grad["cluster_M200_Msun"].clip(lower=1))
        plot_conditional(axes[1, 1], "log_cluster_M200",
                         r"$\log M_{200}$")
    else:
        axes[1, 1].text(0.5, 0.5, "M200 not\navailable",
                        transform=axes[1, 1].transAxes, ha="center")
        axes[1, 1].set_title(r"Binned by $M_{200}$")

    if "vel_disp_kms" in grad.columns:
        grad["log_sigma_v"] = np.log10(grad["vel_disp_kms"].clip(lower=1))
        plot_conditional(axes[1, 2], "log_sigma_v",
                         r"$\log\sigma_v$")
    else:
        axes[1, 2].text(0.5, 0.5, "sigma_v not\navailable",
                        transform=axes[1, 2].transAxes, ha="center")

    fig.suptitle(r"$\nabla$[Fe/H] vs $M_*$: conditional dependences (z=0)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path3 = os.path.join(outdir, "diagnostic_conditional_z0.png")
    fig.savefig(path3, dpi=300, bbox_inches="tight")
    print(f"Saved {path3}")
    plt.close(fig)

    profiles.close()

    # Summary stats
    print(f"\n--- Summary Statistics (z=0) ---")
    print(f"Total galaxies: {len(grad)}")
    print(f"Mass range: {grad['M_star_Msun'].min():.2e} - "
          f"{grad['M_star_Msun'].max():.2e} Msun")
    for key in ["slope_FeH_star", "slope_OFe_star", "slope_MgFe_star"]:
        if key in grad.columns:
            vals = grad[key].dropna()
            print(f"{key}: median={vals.median():.4f}, "
                  f"mean={vals.mean():.4f}, std={vals.std():.4f} (N={len(vals)})")
    for key in ["slope_inner_FeH_star", "slope_outer_FeH_star"]:
        if key in grad.columns:
            vals = grad[key].dropna()
            print(f"{key}: median={vals.median():.4f}, std={vals.std():.4f}")


if __name__ == "__main__":
    main()
