"""
Step 5: Redshift evolution of metallicity profiles and gradients.

At each snapshot (z=0 to z~3), independently:
  1. Downloads group catalog + offsets (cached after first run).
  2. Downloads ~352 snapshot chunks per snapshot (bulk approach).
  3. Computes Sersic indices and selects ellipticals.
  4. Saves full star+gas particle data for ellipticals.
  5. Computes radial metallicity profiles and fits gradients.

No individual galaxy progenitor tracking -- each redshift is an
independent elliptical population.

API-efficient: ~352 file downloads per snapshot instead of ~64,000
per-subhalo cutout API calls.

Output: catalogs/elliptical_catalog_snap{NN}.csv  (per-snapshot)
        catalogs/gradients_snap{NN}.csv           (per-snapshot)
        catalogs/gradients_evolution.csv           (combined)
        data/snap{NN}/subhalo_{ID}.hdf5           (particle data)
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.api_utils import select_and_extract_bulk
from scripts.compute_profiles_lib import (
    process_galaxy, fit_gradient, compute_sersic_index,
)


# ── Profile computation ───────────────────────────────────────────────────

def compute_one_profile(filepath, sub_id, r_half_pkpc, z,
                        halo_id, mstar, snap_num, subhalo_pos=None):
    """Compute profiles and gradients for one galaxy."""
    result = process_galaxy(filepath, sub_id, r_half_pkpc, z,
                            subhalo_pos=subhalo_pos)
    if result is None:
        return None

    bins = result["r_bins"]
    grad_row = {
        "subhalo_id": sub_id,
        "halo_id": halo_id,
        "snap_num": snap_num,
        "redshift": z,
        "M_star_Msun": mstar,
        "R_half_pkpc": r_half_pkpc,
    }

    star_keys = ["FeH_star", "OFe_star", "MgFe_star", "Z_star",
                 "NO_star", "CO_star", "CN_star"]
    for key in star_keys:
        if key in result:
            s, ic, se = fit_gradient(
                bins, result[key], cfg.PROFILE_FIT_RMIN, cfg.PROFILE_FIT_RMAX)
            grad_row[f"slope_{key}"] = s
            grad_row[f"intercept_{key}"] = ic
            grad_row[f"slope_err_{key}"] = se

            si, ici, sei = fit_gradient(
                bins, result[key],
                cfg.PROFILE_FIT_INNER_MIN, cfg.PROFILE_FIT_RBREAK)
            grad_row[f"slope_inner_{key}"] = si
            grad_row[f"slope_err_inner_{key}"] = sei

            so, ico, seo = fit_gradient(
                bins, result[key],
                cfg.PROFILE_FIT_RBREAK, cfg.PROFILE_FIT_OUTER_MAX)
            grad_row[f"slope_outer_{key}"] = so
            grad_row[f"slope_err_outer_{key}"] = seo

    if result.get("has_gas", False):
        for key in ["FeH_gas", "Z_gas", "SFR_profile"]:
            if key in result:
                s, ic, se = fit_gradient(
                    bins, result[key],
                    cfg.PROFILE_FIT_RMIN, cfg.PROFILE_FIT_RMAX)
                grad_row[f"slope_{key}"] = s
                grad_row[f"intercept_{key}"] = ic
                grad_row[f"slope_err_{key}"] = se

    for gkey in ["FeH_star_global", "OFe_star_global", "MgFe_star_global",
                 "NO_star_global", "CO_star_global", "CN_star_global",
                 "age_star_global"]:
        if gkey in result:
            grad_row[gkey] = result[gkey]

    return grad_row


# ── Per-snapshot pipeline ─────────────────────────────────────────────────

def process_snapshot(snap_num, basedir, catdir, cpu_workers=None):
    """Full pipeline for one snapshot: select + extract + profiles."""
    z = cfg.SNAPSHOTS[snap_num]
    datadir = os.path.join(basedir, cfg.DATA_DIR, f"snap{snap_num:03d}")
    gcdir = os.path.join(basedir, cfg.DATA_DIR, "groupcats")
    os.makedirs(datadir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Snapshot {snap_num} (z={z:.2f})")
    print(f"{'='*60}")

    grad_file = os.path.join(catdir, f"gradients_snap{snap_num:03d}.csv")
    if os.path.exists(grad_file):
        print(f"  Gradients already exist at {grad_file}, loading...")
        grad_df = pd.read_csv(grad_file)
        print(f"  Loaded {len(grad_df)} existing gradient measurements")
        return grad_df.to_dict("records")

    # Select ellipticals + extract particles (bulk snapshot)
    ellip_file = os.path.join(
        catdir, f"elliptical_catalog_snap{snap_num:03d}.csv")
    if os.path.exists(ellip_file):
        ellip_df = pd.read_csv(ellip_file)
        print(f"  Loaded {len(ellip_df)} cached ellipticals from {ellip_file}")
    else:
        ellip_df = select_and_extract_bulk(
            snap_num, gcdir, datadir, compute_sersic_index)
        if len(ellip_df) == 0:
            print(f"  No ellipticals found, skipping")
            return []
        ellip_df.to_csv(ellip_file, index=False)
        print(f"  Saved {len(ellip_df)} ellipticals to {ellip_file}")

    if len(ellip_df) == 0:
        return []

    # Compute profiles (CPU-bound, parallel)
    if cpu_workers is None:
        cpu_workers = min(os.cpu_count() or 4, 12)

    print(f"\n  Computing profiles ({len(ellip_df)} galaxies, "
          f"{cpu_workers} workers)...")

    a = 1.0 / (1.0 + z)
    jobs = []
    for _, row in ellip_df.iterrows():
        sid = int(row["subhalo_id"])
        filepath = os.path.join(datadir, f"subhalo_{sid}.hdf5")
        if not os.path.exists(filepath):
            continue
        sub_pos = [row["pos_x"], row["pos_y"], row["pos_z"]]
        jobs.append((filepath, sid, row["R_half_pkpc"], z,
                      int(row["halo_id"]), row["M_star_Msun"],
                      snap_num, sub_pos))

    print(f"  {len(jobs)} particle files found "
          f"({len(ellip_df) - len(jobs)} missing)")

    gradient_rows = []
    done = 0

    with ProcessPoolExecutor(max_workers=cpu_workers) as executor:
        futures = {executor.submit(compute_one_profile, *job): job[1]
                   for job in jobs}
        for future in as_completed(futures):
            done += 1
            try:
                grad_row = future.result()
                if grad_row is not None:
                    gradient_rows.append(grad_row)
            except Exception as e:
                print(f"    Profile error: {e}")
            if done % 50 == 0:
                print(f"    {done}/{len(jobs)} profiles, "
                      f"{len(gradient_rows)} successful")

    print(f"  {len(gradient_rows)} profiles computed at z={z:.2f}")

    if gradient_rows:
        grad_df = pd.DataFrame(gradient_rows)
        grad_df.to_csv(grad_file, index=False)
        print(f"  Saved to {grad_file}")

    return gradient_rows


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Redshift evolution of metallicity profiles "
                    "(independent populations, bulk snapshot approach)"
    )
    parser.add_argument("--snaps", type=int, nargs="+",
                        default=[91, 84, 78, 67, 50, 33, 25],
                        help="Snapshot numbers to process (default: z=0.1..3)")
    parser.add_argument("--cpu-workers", type=int, default=None)
    args = parser.parse_args()

    snaps_to_process = []
    for s in args.snaps:
        if s not in cfg.SNAPSHOTS:
            print(f"Unknown snapshot {s}. "
                  f"Available: {list(cfg.SNAPSHOTS.keys())}")
            sys.exit(1)
        if cfg.SNAPSHOTS[s] > 3.0:
            print(f"Snapshot {s} (z={cfg.SNAPSHOTS[s]:.1f}) beyond z=3, "
                  f"skipping")
        else:
            snaps_to_process.append(s)

    basedir = os.path.join(os.path.dirname(__file__), "..")
    catdir = os.path.join(basedir, cfg.CATALOG_DIR)

    all_gradients = []

    grad_z0 = os.path.join(catdir, "gradients_snap099.csv")
    if os.path.exists(grad_z0):
        g0 = pd.read_csv(grad_z0)
        all_gradients.append(g0)
        print(f"Loaded {len(g0)} z=0 gradient measurements")

    for target_snap in snaps_to_process:
        rows = process_snapshot(
            target_snap, basedir, catdir, args.cpu_workers)
        if rows:
            all_gradients.append(pd.DataFrame(rows))

    if all_gradients:
        combined = pd.concat(all_gradients, ignore_index=True)
        combined_file = os.path.join(catdir, "gradients_evolution.csv")
        combined.to_csv(combined_file, index=False)
        print(f"\nCombined evolution: {len(combined)} rows -> {combined_file}")


if __name__ == "__main__":
    main()
