"""
Step 6: Recompute and save full radial profile HDF5 files for all snapshots.

The evolution pipeline (step 5) computed gradients but discarded the full
profile arrays.  This script re-processes the already-downloaded cutouts
using the gradient CSVs as the galaxy list and saves profiles_snap{NN}.hdf5
for every snapshot that has data.

Usage:
    python 06_save_profiles_all_snaps.py [--snaps 67 50 33 25] [--workers 8]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import h5py
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.compute_profiles_lib import process_galaxy

BASEDIR = os.path.join(os.path.dirname(__file__), "..")
CATDIR = os.path.join(BASEDIR, cfg.CATALOG_DIR)


def process_one(filepath, sub_id, r_half, z):
    result = process_galaxy(filepath, sub_id, r_half, z)
    return sub_id, result


def run_snapshot(snap, n_workers):
    z = cfg.SNAPSHOTS[snap]
    datadir = os.path.join(BASEDIR, cfg.DATA_DIR, f"snap{snap:03d}")
    profile_file = os.path.join(CATDIR, f"profiles_snap{snap:03d}.hdf5")

    if os.path.exists(profile_file):
        with h5py.File(profile_file, "r") as hf:
            n_existing = len(hf.keys())
        print(f"  snap{snap:03d} (z={z:.1f}): profiles already exist "
              f"({n_existing} galaxies), skipping")
        return

    grad_file = os.path.join(CATDIR, f"gradients_snap{snap:03d}.csv")
    if not os.path.exists(grad_file):
        print(f"  snap{snap:03d}: no gradient CSV found, skipping")
        return

    df = pd.read_csv(grad_file)
    print(f"  snap{snap:03d} (z={z:.1f}): {len(df)} galaxies to process")

    jobs = []
    for _, row in df.iterrows():
        sub_id = int(row["subhalo_id"])
        filepath = os.path.join(datadir, f"subhalo_{sub_id}.hdf5")
        if os.path.exists(filepath):
            jobs.append((filepath, sub_id, row["R_half_pkpc"], z))

    if not jobs:
        print(f"    No HDF5 cutout files found in {datadir}")
        return

    print(f"    Found {len(jobs)} cutout files, computing profiles...")

    all_profiles = {}
    done = 0

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(process_one, *job): job[1]
            for job in jobs
        }
        for future in as_completed(futures):
            done += 1
            try:
                sub_id, result = future.result()
            except Exception as e:
                print(f"    ERROR: {e}")
                continue
            if result is not None:
                all_profiles[sub_id] = result
            if done % 100 == 0:
                print(f"    [{done}/{len(jobs)}] {len(all_profiles)} successful")

    with h5py.File(profile_file, "w") as hf:
        for sub_id, prof in all_profiles.items():
            grp = hf.create_group(str(sub_id))
            for key, val in prof.items():
                if isinstance(val, np.ndarray):
                    grp.create_dataset(key, data=val)
                elif isinstance(val, (int, float, bool)):
                    grp.attrs[key] = val

    print(f"    Saved {len(all_profiles)} profiles to {profile_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Save full profile HDF5 files for all snapshots"
    )
    all_snaps = sorted(cfg.SNAPSHOTS.keys(), reverse=True)
    parser.add_argument("--snaps", type=int, nargs="+", default=all_snaps,
                        help=f"Snapshots to process (default: all)")
    parser.add_argument("--workers", type=int,
                        default=min(os.cpu_count() or 4, 16))
    parser.add_argument("--force", action="store_true",
                        help="Re-process even if profile HDF5 already exists")
    args = parser.parse_args()

    print(f"Processing {len(args.snaps)} snapshots with {args.workers} workers\n")

    for snap in args.snaps:
        if snap not in cfg.SNAPSHOTS:
            print(f"  Unknown snapshot {snap}, skipping")
            continue
        if args.force:
            profile_file = os.path.join(CATDIR, f"profiles_snap{snap:03d}.hdf5")
            if os.path.exists(profile_file):
                os.remove(profile_file)
        run_snapshot(snap, args.workers)


if __name__ == "__main__":
    main()
