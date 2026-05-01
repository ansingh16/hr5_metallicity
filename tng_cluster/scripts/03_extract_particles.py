"""
Step 3: Download particle data for selected elliptical galaxies.

For each elliptical galaxy, downloads stellar and gas particle cutouts
from the TNG API and saves them as individual HDF5 files.

Uses parallel downloads (--workers N, default 4) for speed.
Resumable: skips already-downloaded files.

Output: data/snap{NN}/subhalo_{ID}.hdf5
"""

import os
import sys
import argparse
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.api_utils import download_cutout


# Fields to download for each particle type
STAR_FIELDS = ",".join([
    "Coordinates",
    "Masses",
    "GFM_InitialMass",
    "GFM_Metallicity",
    "GFM_Metals",
    "GFM_StellarFormationTime",
    "Velocities",
    "ParticleIDs",
])

GAS_FIELDS = ",".join([
    "Coordinates",
    "Masses",
    "GFM_Metallicity",
    "GFM_Metals",
    "StarFormationRate",
    "Velocities",
    "Density",
    "InternalEnergy",
    "ElectronAbundance",
    "ParticleIDs",
])


def download_one(snap_num, subhalo_id, save_dir):
    """
    Download stellar and gas particle cutouts for one galaxy.
    Returns (subhalo_id, success_bool, file_size_MB_or_error_msg).
    """
    outpath = os.path.join(save_dir, f"subhalo_{subhalo_id}.hdf5")

    # Skip if already downloaded
    if os.path.exists(outpath):
        size_mb = os.path.getsize(outpath) / 1e6
        return subhalo_id, "skipped", size_mb

    particle_fields = {
        "stars": STAR_FIELDS,
        "gas": GAS_FIELDS,
    }

    try:
        filename = download_cutout(snap_num, subhalo_id, particle_fields,
                                   save_dir=save_dir)
        # Rename to consistent name
        if filename != outpath:
            os.rename(filename, outpath)
        size_mb = os.path.getsize(outpath) / 1e6
        return subhalo_id, "ok", size_mb
    except Exception as e:
        return subhalo_id, "error", str(e)


def main():
    parser = argparse.ArgumentParser(
        description="Download particle data for elliptical galaxies"
    )
    parser.add_argument("--snap", type=int, default=99)
    parser.add_argument("--start", type=int, default=0,
                        help="Start index in catalog (for resuming)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Max number of galaxies to download")
    parser.add_argument("--workers", type=int, default=4,
                        help="Parallel downloads (default: 4)")
    args = parser.parse_args()

    catdir = os.path.join(os.path.dirname(__file__), "..", cfg.CATALOG_DIR)
    ellip_file = os.path.join(catdir, f"elliptical_catalog_snap{args.snap:03d}.csv")

    if not os.path.exists(ellip_file):
        print(f"Elliptical catalog not found: {ellip_file}")
        print("Run 02_select_ellipticals.py first.")
        sys.exit(1)

    df = pd.read_csv(ellip_file)
    print(f"Loaded {len(df)} ellipticals from {ellip_file}")

    save_dir = os.path.join(os.path.dirname(__file__), "..",
                            cfg.DATA_DIR, f"snap{args.snap:03d}")
    os.makedirs(save_dir, exist_ok=True)

    end = args.start + args.limit if args.limit else len(df)
    subset = df.iloc[args.start:end]
    total = len(subset)

    print(f"Downloading {total} galaxies with {args.workers} workers")
    print(f"Save directory: {save_dir}\n")

    downloaded = 0
    skipped = 0
    errors = 0
    total_mb = 0.0
    done = 0
    lock = threading.Lock()

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_one, args.snap, int(row["subhalo_id"]), save_dir): (
                int(row["subhalo_id"]), row["M_star_Msun"]
            )
            for _, row in subset.iterrows()
        }

        for future in as_completed(futures):
            sub_id, mstar = futures[future]
            sid, status, info = future.result()
            done += 1

            with lock:
                if status == "ok":
                    downloaded += 1
                    total_mb += info
                    print(f"  [{done}/{total}] Subhalo {sid}: "
                          f"{info:.1f} MB (M*={mstar:.2e})")
                elif status == "skipped":
                    skipped += 1
                    total_mb += info
                    if done % 50 == 0:
                        print(f"  [{done}/{total}] ... {skipped} skipped so far")
                else:
                    errors += 1
                    print(f"  [{done}/{total}] Subhalo {sid}: ERROR - {info}")

    print(f"\nDone: {downloaded} downloaded, {skipped} skipped, {errors} errors")
    print(f"Total data: {total_mb:.0f} MB")


if __name__ == "__main__":
    main()
