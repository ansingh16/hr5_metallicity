#!/usr/bin/env python3
"""
Generate TNG-Cluster figures matching the HR5 paper figures.

Figure 1 (metal_distribution.png):
    3-panel vertical: individual [Fe/H], Z_star, Z_gas profiles + median at z=0

Figure 2 (metal_evolution_redshift.png):
    3-panel vertical: median [Fe/H], Z_star, Z_gas profiles at different redshifts

Figure 3 (metal_distribution_mass.png):
    3-panel vertical: median [Fe/H], Z_star, Z_gas profiles by mass bin at z=0

Figure 4 (slope_evolution.png):
    3-panel vertical: slope of [Fe/H], Z_star, Z_gas vs redshift

Figure 5 (OFe_vs_Fe_evolution.png):
    Median [O/Fe] vs [Fe/H] at different redshifts

Figure 6 (SFR_redshift.png):
    SFR median profile at different redshifts

Figure 7 (slope_evolution_SFR.png):
    SFR slope vs redshift
"""

import os
import sys
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

style_path = os.path.join(os.path.dirname(__file__), "..", "paper.mplstyle")
plt.style.use(style_path)

BASEDIR = os.path.join(os.path.dirname(__file__), "..")
CATDIR = os.path.join(BASEDIR, cfg.CATALOG_DIR)
OUTDIR = os.path.join(BASEDIR, cfg.PLOT_DIR)

MASS_BINS = [
    (2e9, 1e10, r"$\log\,M_\star(M_\odot)$: 9--10"),
    (1e10, 1e11, r"$\log\,M_\star(M_\odot)$: 10--11"),
    (1e11, 1e13, r"$\log\,M_\star(M_\odot)$: 11--12"),
]

REDSHIFT_COLORS = {
    0.0: "#0072B2",
    0.5: "#D55E00",
    1.0: "#009E73",
    2.0: "#E69F00",
    3.0: "#CC79A7",
}


def load_profiles(snap=99):
    path = os.path.join(CATDIR, f"profiles_snap{snap:03d}.hdf5")
    return h5py.File(path, "r")


def load_gradients():
    return pd.read_csv(os.path.join(CATDIR, "gradients_evolution.csv"))


def load_gradients_snap(snap):
    return pd.read_csv(os.path.join(CATDIR, f"gradients_snap{snap:03d}.csv"))


def collect_profiles(profiles, grad, field):
    """Collect all radial profiles for a given field from the HDF5."""
    all_prof = []
    sub_ids = grad["subhalo_id"].values
    for sid in sub_ids:
        key = str(sid)
        if key not in profiles or field not in profiles[key]:
            continue
        prof = profiles[key][field][:]
        if not np.all(np.isnan(prof)):
            all_prof.append(prof)
    return np.array(all_prof) if all_prof else None


def median_profile(profiles, grad, field, mass_lo=None, mass_hi=None):
    """Compute median + 16/84th percentile of radial profiles."""
    if mass_lo is not None and mass_hi is not None:
        mask = (grad["M_star_Msun"] >= mass_lo) & (grad["M_star_Msun"] < mass_hi)
        sub_ids = grad.loc[mask, "subhalo_id"].values
    else:
        sub_ids = grad["subhalo_id"].values

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
    return (np.nanmedian(arr, axis=0),
            np.nanpercentile(arr, 16, axis=0),
            np.nanpercentile(arr, 84, axis=0),
            len(all_prof))


def reconstructed_profile(grad_snap, slope_col, intercept_col, r_bins,
                          mass_lo=None, mass_hi=None):
    """Reconstruct median linear profile from slope + intercept."""
    mask = grad_snap[slope_col].notna()
    if mass_lo is not None:
        mask &= (grad_snap["M_star_Msun"] >= mass_lo)
    if mass_hi is not None:
        mask &= (grad_snap["M_star_Msun"] < mass_hi)
    sub = grad_snap[mask]
    if len(sub) < 3:
        return None
    slopes = sub[slope_col].values
    intercepts = sub[intercept_col].values
    lines = slopes[:, None] * r_bins[None, :] + intercepts[:, None]
    return np.nanmedian(lines, axis=0)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 1: metal_distribution.png
# ═══════════════════════════════════════════════════════════════════════════
def figure1():
    profiles = load_profiles()
    grad = load_gradients_snap(99)
    first_key = list(profiles.keys())[0]
    r_bins = profiles[first_key]["r_bins"][:]

    fig, axes = plt.subplots(3, 1, figsize=(3.321, 7.0), sharex=True)

    panels = [
        ("FeH_star", r"[Fe/H]", "Median Fe/H"),
        ("Z_star", r"$\log\,Z_{\rm star}/Z_\odot$", r"Median $Z_{\rm star}$"),
        ("Z_gas", r"$\log\,Z_{\rm gas}/Z_\odot$", r"Median $Z_{\rm gas}$"),
    ]

    for ax, (field, ylabel, med_label) in zip(axes, panels):
        all_prof = collect_profiles(profiles, grad, field)
        if all_prof is not None:
            for i in range(len(all_prof)):
                ax.plot(r_bins, all_prof[i], color="0.8", lw=0.3, alpha=0.5)
            med = np.nanmedian(all_prof, axis=0)
            ax.plot(r_bins, med, "o-", color="#E63946", ms=3, lw=1.5,
                    zorder=5, label=med_label)
        ax.set_ylabel(ylabel)
        ax.legend(loc="upper right")
        ax.set_xlim(0, 3)

    axes[-1].set_xlabel(r"$R/R_{\rm half}$")
    fig.tight_layout()
    path = os.path.join(OUTDIR, "metal_distribution.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)
    profiles.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2: metal_evolution_redshift.png
# ═══════════════════════════════════════════════════════════════════════════
def figure2():
    snap_z = {99: 0.0, 67: 0.5, 50: 1.0, 33: 2.0, 25: 3.0}

    fig, axes = plt.subplots(3, 1, figsize=(3.321, 7.0), sharex=True)

    panels = [
        ("FeH_star", r"[Fe/H]"),
        ("Z_star", r"$\log\,Z_{\rm star}/Z_\odot$"),
        ("Z_gas", r"$\log\,Z_{\rm gas}/Z_\odot$"),
    ]

    for ax, (field, ylabel) in zip(axes, panels):
        for snap, z in sorted(snap_z.items(), key=lambda x: x[1]):
            color = REDSHIFT_COLORS[z]
            prof_path = os.path.join(CATDIR, f"profiles_snap{snap:03d}.hdf5")
            if not os.path.exists(prof_path):
                continue
            profiles = h5py.File(prof_path, "r")
            grad = load_gradients_snap(snap)
            first_key = list(profiles.keys())[0]
            r_bins = profiles[first_key]["r_bins"][:]

            med, _, _, n = median_profile(profiles, grad, field)
            profiles.close()
            if med is not None:
                ax.plot(r_bins, med, "o-", color=color, ms=3, lw=1.2,
                        label=f"z={z:.1f} (N={n})")

        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 3)

    axes[0].legend(loc="lower left", title="Redshift")
    axes[-1].set_xlabel(r"$R/R_{\rm half}$")
    fig.tight_layout()
    path = os.path.join(OUTDIR, "metal_evolution_redshift.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3: metal_distribution_mass.png
# ═══════════════════════════════════════════════════════════════════════════
def figure3():
    profiles = load_profiles()
    grad = load_gradients_snap(99)
    first_key = list(profiles.keys())[0]
    r_bins = profiles[first_key]["r_bins"][:]

    fig, axes = plt.subplots(3, 1, figsize=(3.321, 7.0), sharex=True)

    panels = [
        ("FeH_star", r"[Fe/H]"),
        ("Z_star", r"$\log\,Z_{\rm star}/Z_\odot$"),
        ("Z_gas", r"$\log\,Z_{\rm gas}/Z_\odot$"),
    ]

    colors = ["#E63946", "#0072B2", "#009E73"]

    for ax, (field, ylabel) in zip(axes, panels):
        for (mlo, mhi, label), color in zip(MASS_BINS, colors):
            med, p16, p84, n = median_profile(profiles, grad, field,
                                               mlo, mhi)
            if med is None:
                continue
            ax.plot(r_bins, med, "o-", color=color, ms=3, lw=1.2,
                    label=f"{label} (N={n})")
            ax.fill_between(r_bins, p16, p84, alpha=0.15, color=color)
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 3)

    axes[0].legend(loc="lower left", fontsize=6)
    axes[-1].set_xlabel(r"$R/R_{\rm half}$")
    fig.tight_layout()
    path = os.path.join(OUTDIR, "metal_distribution_mass.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)
    profiles.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure 4: slope_evolution.png
# ═══════════════════════════════════════════════════════════════════════════
def figure4():
    df = load_gradients()

    fig, axes = plt.subplots(3, 1, figsize=(3.321, 7.0), sharex=True)

    panels = [
        ("slope_FeH_star", r"[Fe/H] slope"),
        ("slope_Z_star", r"$Z_{\rm star}$ slope"),
        ("slope_Z_gas", r"$Z_{\rm gas}$ slope"),
    ]

    redshifts = sorted(df["redshift"].unique(), reverse=True)

    for ax, (col, ylabel) in zip(axes, panels):
        # Individual slopes as scatter
        valid = df[col].notna()
        ax.scatter(df.loc[valid, "redshift"], df.loc[valid, col],
                   s=6, color="#CC79A7", alpha=0.3, marker="x",
                   linewidths=0.5, zorder=1)

        # Median + 25-75th percentile at each redshift
        meds, q25s, q75s, zs = [], [], [], []
        for z in redshifts:
            vals = df.loc[(df["redshift"] == z) & valid, col].values
            if len(vals) < 5:
                continue
            meds.append(np.median(vals))
            q25s.append(np.percentile(vals, 25))
            q75s.append(np.percentile(vals, 75))
            zs.append(z)

        ax.fill_between(zs, q25s, q75s, alpha=0.3, color="#7EB2DD",
                         label="25th--75th Quantiles")
        ax.plot(zs, meds, "-", color="#E63946", lw=1.8,
                label="running median", zorder=3)

        ax.set_ylabel(ylabel)
        ax.invert_xaxis()

    axes[0].legend(loc="lower left")
    axes[-1].set_xlabel("Redshift")
    fig.tight_layout()
    path = os.path.join(OUTDIR, "slope_evolution.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 5: OFe_vs_Fe_evolution.png
# ═══════════════════════════════════════════════════════════════════════════
def figure5():
    df = load_gradients()

    fig, ax = plt.subplots(figsize=(3.321, 3.321))

    redshifts = sorted(df["redshift"].unique())
    cmap = plt.cm.tab10
    colors = {z: cmap(i) for i, z in enumerate(redshifts)}

    for z in redshifts:
        sub = df[df["redshift"] == z]
        feh = sub["FeH_star_global"].dropna()
        ofe = sub["OFe_star_global"].dropna()
        common = sub.loc[feh.index.intersection(ofe.index)]
        if len(common) < 5:
            continue
        med_feh = common["FeH_star_global"].median()
        med_ofe = common["OFe_star_global"].median()
        ax.scatter(med_feh, med_ofe, s=60, color=colors[z], marker="*",
                   edgecolors="k", linewidths=0.3, zorder=5,
                   label=f"z={z:.1f}")

    ax.set_xlabel("[Fe/H]")
    ax.set_ylabel("[O/Fe]")
    ax.legend(loc="lower left", fontsize=6, title="Redshift",
              title_fontsize=7)
    fig.tight_layout()
    path = os.path.join(OUTDIR, "OFe_vs_Fe_evolution.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 6: SFR_redshift.png
# ═══════════════════════════════════════════════════════════════════════════
def figure6():
    snap_z = {99: 0.0, 67: 0.5, 50: 1.0, 33: 2.0, 25: 3.0}

    fig, ax = plt.subplots(figsize=(3.321, 3.321))

    for snap, z in sorted(snap_z.items(), key=lambda x: x[1]):
        color = REDSHIFT_COLORS[z]
        prof_path = os.path.join(CATDIR, f"profiles_snap{snap:03d}.hdf5")
        if not os.path.exists(prof_path):
            continue
        profiles = h5py.File(prof_path, "r")
        grad = load_gradients_snap(snap)
        first_key = list(profiles.keys())[0]
        r_bins = profiles[first_key]["r_bins"][:]

        med, _, _, n = median_profile(profiles, grad, "SFR_profile")
        profiles.close()
        if med is not None:
            ax.plot(r_bins, med, "o-", color=color, ms=3, lw=1.2,
                    label=f"z={z:.1f} (N={n})")

    ax.set_xlabel(r"$R/R_{\rm half}$")
    ax.set_ylabel(r"SFR ($M_\odot\,{\rm yr}^{-1}$)")
    ax.set_xlim(0, 3)
    ax.legend(loc="upper right", title="Redshift")
    fig.tight_layout()
    path = os.path.join(OUTDIR, "SFR_redshift.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 7: slope_evolution_SFR.png
# ═══════════════════════════════════════════════════════════════════════════
def figure7():
    df = load_gradients()

    fig, ax = plt.subplots(figsize=(3.321, 3.321))

    col = "slope_SFR_profile"
    valid = df[col].notna()
    ax.scatter(df.loc[valid, "redshift"], df.loc[valid, col],
               s=6, color="#CC79A7", alpha=0.3, marker="x",
               linewidths=0.5, zorder=1)

    redshifts = sorted(df["redshift"].unique(), reverse=True)
    meds, q25s, q75s, zs = [], [], [], []
    for z in redshifts:
        vals = df.loc[(df["redshift"] == z) & valid, col].values
        if len(vals) < 5:
            continue
        meds.append(np.median(vals))
        q25s.append(np.percentile(vals, 25))
        q75s.append(np.percentile(vals, 75))
        zs.append(z)

    ax.fill_between(zs, q25s, q75s, alpha=0.3, color="#7EB2DD",
                     label="25th--75th Quantiles")
    ax.plot(zs, meds, "-", color="#E63946", lw=1.8,
            label="running median", zorder=3)

    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"SFR ($M_\odot\,{\rm yr}^{-1}$) slope")
    ax.legend(loc="lower left")
    ax.invert_xaxis()
    fig.tight_layout()
    path = os.path.join(OUTDIR, "slope_evolution_SFR.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("Generating HR5-equivalent figures for TNG-Cluster\n")

    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()
    figure7()

    print("\nAll figures saved to", OUTDIR)


if __name__ == "__main__":
    main()
