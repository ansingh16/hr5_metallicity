"""
Step 1: Select galaxy clusters from TNG-Cluster.

Uses the group catalog HDF5 files (bulk download) to identify primary
zoom-target clusters via the GroupPrimaryZoomTarget field.  No per-halo
API calls needed.

Output: catalogs/cluster_catalog_snap{NN}.csv
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.api_utils import (
    download_group_catalog, load_group_catalog, find_clusters_from_catalog,
)


def select_clusters(snap_num, gcdir):
    z = cfg.SNAPSHOTS[snap_num]
    print(f"Selecting clusters at snapshot {snap_num} (z={z:.2f})")

    compact_path = download_group_catalog(snap_num, gcdir)
    sub, grp = load_group_catalog(compact_path)

    cluster_ids = find_clusters_from_catalog(grp)
    print(f"Found {len(cluster_ids)} primary zoom-target clusters (M200 >= "
          f"{cfg.MIN_GROUP_M200_MSUN:.1e} Msun)")

    a = 1.0 / (1.0 + z)
    clusters = []
    for cid in cluster_ids:
        m200_msun = float(grp["Group_M_Crit200"][cid] * 1e10 / cfg.h)
        r200_pkpc = float(grp["Group_R_Crit200"][cid] * a / cfg.h)
        pos = grp["GroupPos"][cid]
        nsubs = int(grp["GroupNsubs"][cid])
        first_sub = int(grp["GroupFirstSub"][cid])

        clusters.append({
            "halo_id": int(cid),
            "snap_num": snap_num,
            "redshift": z,
            "M200_Msun": m200_msun,
            "R200_pkpc": r200_pkpc,
            "pos_x": float(pos[0]),
            "pos_y": float(pos[1]),
            "pos_z": float(pos[2]),
            "n_subhalos": nsubs,
            "first_subhalo_id": first_sub,
            "is_primary_zoom": 1,
        })

    df = pd.DataFrame(clusters)
    print(f"\nSelected {len(df)} clusters")
    if len(df) > 0:
        print(f"  Mass range: {df['M200_Msun'].min():.2e} -- "
              f"{df['M200_Msun'].max():.2e} Msun")
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

    basedir = os.path.join(os.path.dirname(__file__), "..")
    catdir = os.path.join(basedir, cfg.CATALOG_DIR)
    gcdir = os.path.join(basedir, cfg.DATA_DIR, "groupcats")
    os.makedirs(catdir, exist_ok=True)

    df = select_clusters(args.snap, gcdir)

    outfile = os.path.join(catdir, f"cluster_catalog_snap{args.snap:03d}.csv")
    df.to_csv(outfile, index=False)
    print(f"Saved to {outfile}")


if __name__ == "__main__":
    main()
