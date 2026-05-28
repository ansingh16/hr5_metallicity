"""
Step 2+3: Select elliptical galaxies and extract particle data.

Uses bulk snapshot file downloads (~352 per snapshot) instead of
per-subhalo cutout API calls (~64,000 per snapshot).  For each
snapshot chunk, extracts stellar particles, computes Sersic index,
and saves full star+gas particle data for ellipticals.

This replaces both the old Step 2 (Sersic selection via cutouts)
and Step 3 (particle extraction via cutouts) with a single pass.

Output: catalogs/elliptical_catalog_snap{NN}.csv
        data/snap{NN}/subhalo_{ID}.hdf5
"""

import os
import sys
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.api_utils import select_and_extract_bulk
from scripts.compute_profiles_lib import compute_sersic_index


def main():
    parser = argparse.ArgumentParser(
        description="Select ellipticals + extract particles (bulk snapshot)")
    parser.add_argument("--snap", type=int, default=99)
    args = parser.parse_args()

    if args.snap not in cfg.SNAPSHOTS:
        print(f"Snapshot {args.snap} not configured. "
              f"Available: {list(cfg.SNAPSHOTS.keys())}")
        sys.exit(1)

    basedir = os.path.join(os.path.dirname(__file__), "..")
    catdir = os.path.join(basedir, cfg.CATALOG_DIR)
    gcdir = os.path.join(basedir, cfg.DATA_DIR, "groupcats")
    datadir = os.path.join(basedir, cfg.DATA_DIR, f"snap{args.snap:03d}")
    os.makedirs(catdir, exist_ok=True)

    z = cfg.SNAPSHOTS[args.snap]
    print(f"Snapshot {args.snap} (z={z:.2f})")
    print(f"Using bulk snapshot approach (352 file downloads)\n")

    df = select_and_extract_bulk(
        args.snap, gcdir, datadir, compute_sersic_index)

    outfile = os.path.join(
        catdir, f"elliptical_catalog_snap{args.snap:03d}.csv")
    df.to_csv(outfile, index=False)
    print(f"\nSaved {len(df)} ellipticals to {outfile}")
    print(f"Particle data in {datadir}/")


if __name__ == "__main__":
    main()
