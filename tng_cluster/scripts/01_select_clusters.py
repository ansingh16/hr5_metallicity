"""
Step 1: Select galaxy clusters from TNG-Cluster.

Iterates through FoF halos (sorted by mass, descending) and selects those
with M200 above the mass threshold. Stops when masses drop below the cut.

Filters out halos with high low-res DM contamination.

Output: catalogs/cluster_catalog_snap{NN}.csv
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.api_utils import get_halo_info


def mass_api_to_msun(mass_api):
    """Convert API units (10^10 Msun/h) to solar masses."""
    return mass_api * 1e10 / cfg.h


def length_api_to_pkpc(length_ckpc_h, z):
    """Convert comoving kpc/h to physical kpc."""
    a = 1.0 / (1.0 + z)
    return length_ckpc_h * a / cfg.h


def select_clusters(snap_num):
    """
    Select clusters by iterating FoF halos from ID 0 upward.
    Halos in TNG are sorted by total mass (descending), so we stop
    once M200 drops below our threshold.
    """
    z = cfg.SNAPSHOTS[snap_num]
    print(f"Selecting clusters at snapshot {snap_num} (z={z:.2f})")
    print(f"Mass threshold: M200 > {cfg.MIN_GROUP_M200_MSUN:.2e} Msun")
    print(f"Contamination threshold: < {cfg.MAX_CONTAMINATION_FRAC}")

    clusters = []
    consecutive_below = 0
    halo_id = 0

    while True:
        try:
            info = get_halo_info(snap_num, halo_id)
        except Exception as e:
            print(f"  Halo {halo_id}: API error ({e}), stopping")
            break

        m200_api = info.get("Group_M_Crit200", 0)
        m200_msun = mass_api_to_msun(m200_api)

        # Stop if we've seen several consecutive halos below threshold
        # (handles occasional non-monotonic ordering from substructure)
        if m200_msun < cfg.MIN_GROUP_M200_MSUN:
            consecutive_below += 1
            if consecutive_below > 10:
                print(f"  Stopped at halo {halo_id} "
                      f"(M200={m200_msun:.2e} below threshold)")
                break
            halo_id += 1
            continue
        consecutive_below = 0

        # Contamination check
        contam = info.get("GroupContaminationFracByMass", 0.0)
        if contam > cfg.MAX_CONTAMINATION_FRAC:
            print(f"  Halo {halo_id}: M200={m200_msun:.2e}, "
                  f"contamination={contam:.3f} -> SKIP")
            halo_id += 1
            continue

        r200 = info.get("Group_R_Crit200", 0)
        r200_pkpc = length_api_to_pkpc(r200, z)
        pos = info.get("GroupPos", [0, 0, 0])
        nsubs = info.get("GroupNsubs", 0)
        first_sub = info.get("GroupFirstSub", -1)
        sfr = info.get("GroupSFR", 0)
        primary = info.get("GroupPrimaryZoomTarget", 0)

        clusters.append({
            "halo_id": halo_id,
            "snap_num": snap_num,
            "redshift": z,
            "M200_Msun": m200_msun,
            "R200_pkpc": r200_pkpc,
            "pos_x": pos[0],
            "pos_y": pos[1],
            "pos_z": pos[2],
            "n_subhalos": nsubs,
            "first_subhalo_id": first_sub,
            "SFR": sfr,
            "contamination_frac": contam,
            "is_primary_zoom": primary,
        })
        print(f"  Halo {halo_id}: M200={m200_msun:.2e} Msun, "
              f"R200={r200_pkpc:.1f} pkpc, Nsubs={nsubs}, "
              f"primary={primary}")

        halo_id += 1

    df = pd.DataFrame(clusters)
    print(f"\nSelected {len(df)} clusters at z={z:.2f}")
    if len(df) > 0:
        print(f"  Mass range: {df['M200_Msun'].min():.2e} -- "
              f"{df['M200_Msun'].max():.2e} Msun")
        print(f"  Primary zoom targets: {df['is_primary_zoom'].sum()}")
    return df


def main():
    parser = argparse.ArgumentParser(description="Select TNG-Cluster clusters")
    parser.add_argument("--snap", type=int, default=99,
                        help="Snapshot number (default: 99, z=0)")
    args = parser.parse_args()

    if args.snap not in cfg.SNAPSHOTS:
        print(f"Snapshot {args.snap} not in configured list. "
              f"Available: {list(cfg.SNAPSHOTS.keys())}")
        sys.exit(1)

    df = select_clusters(args.snap)

    outdir = os.path.join(os.path.dirname(__file__), "..", cfg.CATALOG_DIR)
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"cluster_catalog_snap{args.snap:03d}.csv")
    df.to_csv(outfile, index=False)
    print(f"Saved to {outfile}")


if __name__ == "__main__":
    main()
