"""
Step 4: Compute radial metallicity profiles for each elliptical galaxy.

Reads particle cutout HDF5 files and computes radial profiles of:
  - [Fe/H]_star   (stellar iron-to-hydrogen ratio)
  - [O/Fe]_star   (stellar alpha enhancement)
  - [Mg/Fe]_star  (stellar magnesium-to-iron ratio)
  - Z_star        (stellar total metallicity)
  - Z_gas         (gas-phase total metallicity)
  - [Fe/H]_gas    (gas-phase iron-to-hydrogen)
  - SFR surface density
  - Stellar age (mass-weighted)

All profiles are in bins of R / R_half.
Uses multiprocessing (--workers N) for CPU-bound profile computation.

Output: catalogs/profiles_snap{NN}.hdf5
        catalogs/gradients_snap{NN}.csv
"""

import os
import sys
import argparse
import h5py
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.compute_profiles_lib import process_galaxy, fit_gradient


def process_one(filepath, sub_id, r_half, z, halo_id, mstar,
                vel_disp=np.nan, r_over_r200=np.nan, cluster_m200=np.nan):
    """
    Process a single galaxy: compute profiles and fit gradients.
    Runs in a worker process. Returns (sub_id, result_dict, grad_row) or
    (sub_id, None, None) on failure.
    """
    result = process_galaxy(filepath, sub_id, r_half, z)
    if result is None:
        return sub_id, None, None

    bins = result["r_bins"]
    grad_row = {
        "subhalo_id": sub_id,
        "halo_id": halo_id,
        "redshift": z,
        "M_star_Msun": mstar,
        "R_half_pkpc": r_half,
        "vel_disp_kms": vel_disp,
        "R_over_R200": r_over_r200,
        "cluster_M200_Msun": cluster_m200,
    }

    # Stellar abundance keys: original + new CNO ratios
    star_keys = ["FeH_star", "OFe_star", "MgFe_star", "Z_star",
                 "NO_star", "CO_star", "CN_star"]

    for key in star_keys:
        if key in result:
            # Full gradient (0.5 - 2.5 R_half)
            s, ic, se = fit_gradient(bins, result[key],
                                     cfg.PROFILE_FIT_RMIN,
                                     cfg.PROFILE_FIT_RMAX)
            grad_row[f"slope_{key}"] = s
            grad_row[f"intercept_{key}"] = ic
            grad_row[f"slope_err_{key}"] = se

            # Inner gradient (0.1 - 1.0 R_half)
            si, ici, sei = fit_gradient(bins, result[key],
                                        cfg.PROFILE_FIT_INNER_MIN,
                                        cfg.PROFILE_FIT_RBREAK)
            grad_row[f"slope_inner_{key}"] = si
            grad_row[f"slope_err_inner_{key}"] = sei

            # Outer gradient (1.0 - 2.5 R_half)
            so, ico, seo = fit_gradient(bins, result[key],
                                        cfg.PROFILE_FIT_RBREAK,
                                        cfg.PROFILE_FIT_OUTER_MAX)
            grad_row[f"slope_outer_{key}"] = so
            grad_row[f"slope_err_outer_{key}"] = seo

    if result.get("has_gas", False):
        for key in ["FeH_gas", "Z_gas", "SFR_profile"]:
            if key in result:
                s, ic, se = fit_gradient(bins, result[key],
                                         cfg.PROFILE_FIT_RMIN,
                                         cfg.PROFILE_FIT_RMAX)
                grad_row[f"slope_{key}"] = s
                grad_row[f"intercept_{key}"] = ic
                grad_row[f"slope_err_{key}"] = se

    # Global properties
    for gkey in ["FeH_star_global", "OFe_star_global", "MgFe_star_global",
                 "NO_star_global", "CO_star_global", "CN_star_global",
                 "age_star_global"]:
        if gkey in result:
            grad_row[gkey] = result[gkey]

    return sub_id, result, grad_row


def main():
    parser = argparse.ArgumentParser(
        description="Compute radial metallicity profiles"
    )
    parser.add_argument("--snap", type=int, default=99)
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel processes (default: CPU count)")
    args = parser.parse_args()

    if args.workers is None:
        args.workers = min(os.cpu_count() or 4, 16)

    z = cfg.SNAPSHOTS[args.snap]
    basedir = os.path.join(os.path.dirname(__file__), "..")
    catdir = os.path.join(basedir, cfg.CATALOG_DIR)
    datadir = os.path.join(basedir, cfg.DATA_DIR, f"snap{args.snap:03d}")

    ellip_file = os.path.join(catdir, f"elliptical_catalog_snap{args.snap:03d}.csv")
    df = pd.read_csv(ellip_file)
    print(f"Processing {len(df)} ellipticals at snap {args.snap} (z={z:.2f})")
    print(f"Using {args.workers} worker processes")

    # Load cluster catalog for R200 and cluster positions
    cluster_file = os.path.join(catdir, f"cluster_catalog_snap{args.snap:03d}.csv")
    cluster_info = {}
    if os.path.exists(cluster_file):
        cluster_df = pd.read_csv(cluster_file)
        a = 1.0 / (1.0 + z)
        for _, cr in cluster_df.iterrows():
            cluster_info[int(cr["halo_id"])] = {
                "R200": cr["R200_pkpc"],
                "M200": cr["M200_Msun"],
                "pos": np.array([cr["pos_x"], cr["pos_y"], cr["pos_z"]]) * a / cfg.h,
            }

    # Build job list (skip missing files upfront)
    jobs = []
    for _, row in df.iterrows():
        sub_id = int(row["subhalo_id"])
        filepath = os.path.join(datadir, f"subhalo_{sub_id}.hdf5")
        if not os.path.exists(filepath):
            continue
        halo_id = int(row["halo_id"])
        vel_disp = row.get("vel_disp_kms", np.nan)
        r_over_r200 = np.nan
        m200 = np.nan
        ci = cluster_info.get(halo_id)
        if ci is not None:
            m200 = ci["M200"]
            gal_pos = np.array([row["pos_x"], row["pos_y"], row["pos_z"]])
            gal_pos_pkpc = gal_pos * a / cfg.h
            dist = np.sqrt(np.sum((gal_pos_pkpc - ci["pos"]) ** 2))
            if ci["R200"] > 0:
                r_over_r200 = dist / ci["R200"]
        jobs.append((filepath, sub_id, row["R_half_pkpc"], z,
                      halo_id, row["M_star_Msun"], vel_disp, r_over_r200, m200))

    print(f"Found {len(jobs)} HDF5 files (skipped {len(df) - len(jobs)} missing)\n")

    all_profiles = {}
    gradient_rows = []
    done = 0
    skipped = 0

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_one, *job): job[1]  # job[1] = sub_id
            for job in jobs
        }

        for future in as_completed(futures):
            sub_id = futures[future]
            done += 1

            try:
                sid, result, grad_row = future.result()
            except Exception as e:
                print(f"  [{done}/{len(jobs)}] Subhalo {sub_id}: ERROR - {e}")
                skipped += 1
                continue

            if result is None:
                skipped += 1
                if done % 100 == 0:
                    print(f"  [{done}/{len(jobs)}] ... "
                          f"{len(gradient_rows)} processed, {skipped} skipped")
                continue

            all_profiles[sid] = result
            gradient_rows.append(grad_row)

            if done % 50 == 0:
                print(f"  [{done}/{len(jobs)}] "
                      f"{len(gradient_rows)} processed, {skipped} skipped")

    print(f"\nComputed profiles for {len(gradient_rows)} galaxies "
          f"({skipped} skipped)")

    # Save profiles to HDF5
    profile_file = os.path.join(catdir, f"profiles_snap{args.snap:03d}.hdf5")
    with h5py.File(profile_file, "w") as hf:
        for sub_id, prof in all_profiles.items():
            grp = hf.create_group(str(sub_id))
            for key, val in prof.items():
                if isinstance(val, np.ndarray):
                    grp.create_dataset(key, data=val)
                elif isinstance(val, (int, float, bool)):
                    grp.attrs[key] = val

    print(f"Saved profiles to {profile_file}")

    grad_df = pd.DataFrame(gradient_rows)
    grad_file = os.path.join(catdir, f"gradients_snap{args.snap:03d}.csv")
    grad_df.to_csv(grad_file, index=False)
    print(f"Saved gradients to {grad_file}")


if __name__ == "__main__":
    main()
