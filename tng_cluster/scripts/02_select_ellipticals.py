"""
Step 2: Identify elliptical galaxies within selected clusters.

Strategy (optimized for speed):
  1. Use the subhalo meta endpoint (fast, ~0.3s) for initial filtering
     on mass, flag, and sSFR — this avoids the slow info.json call.
  2. Only fetch full info.json (with element fractions) for galaxies
     that pass all selection cuts.
  3. Parallelize with ThreadPoolExecutor.

Selection criteria:
  - Stellar mass > MIN_STELLAR_MASS_MSUN
  - SubhaloFlag == 1 (cosmological origin)
  - sSFR < 10^-11 yr^-1 (quiescent / elliptical proxy)
  - Subhalo belongs to a selected cluster (grnr matches)

Output: catalogs/elliptical_catalog_snap{NN}.csv
"""

import os
import sys
import argparse
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.api_utils import get_subhalo_meta, get_subhalo_info


def mass_api_to_msun(val):
    """10^10 Msun/h -> Msun."""
    return val * 1e10 / cfg.h


def length_api_to_pkpc(val, z):
    """Comoving kpc/h -> physical kpc."""
    a = 1.0 / (1.0 + z)
    return val * a / cfg.h


def prefilter_subhalo(snap_num, sub_id, cluster_halo_ids, z):
    """
    Quick check using the meta endpoint (~0.3s) to see if this subhalo
    passes basic cuts. Returns (sub_id, True/False, meta_dict).
    """
    try:
        meta = get_subhalo_meta(snap_num, sub_id)
    except Exception:
        return sub_id, False, None

    # Must be cosmological
    if not meta.get("subhaloflag", 0):
        return sub_id, False, None

    # Must belong to a selected cluster
    grnr = meta.get("grnr", -1)
    if grnr not in cluster_halo_ids:
        return sub_id, False, None

    # Stellar mass cut (mass_stars is in 10^10 Msun/h)
    mstar_api = meta.get("mass_stars", 0)
    mstar_msun = mass_api_to_msun(mstar_api)
    if mstar_msun < cfg.MIN_STELLAR_MASS_MSUN:
        return sub_id, False, None

    # Star particle count
    nstars = meta.get("len_stars", 0)
    if nstars < cfg.MIN_STAR_PARTICLES:
        return sub_id, False, None

    # sSFR cut
    sfr = meta.get("sfr", 0.0)
    ssfr = sfr / mstar_msun if mstar_msun > 0 else 1.0
    if ssfr >= 1e-11:
        return sub_id, False, None

    return sub_id, True, meta


def fetch_full_info(snap_num, sub_id, meta, z, cluster_lookup):
    """
    For a subhalo that passed prefiltering, fetch the full info.json
    to get element fractions.
    """
    try:
        info = get_subhalo_info(snap_num, sub_id)
    except Exception:
        return None

    grnr = meta["grnr"]
    mstar_msun = mass_api_to_msun(meta["mass_stars"])
    sfr = meta.get("sfr", 0.0)
    ssfr = sfr / mstar_msun

    rhalf_pkpc = length_api_to_pkpc(meta.get("halfmassrad_stars", 0), z)
    pos_x = meta.get("pos_x", 0)
    pos_y = meta.get("pos_y", 0)
    pos_z = meta.get("pos_z", 0)
    vel_disp = meta.get("veldisp", 0)

    z_star = info.get("SubhaloStarMetallicity", 0)
    z_gas = info.get("SubhaloGasMetallicity", 0)
    star_metals = info.get("SubhaloStarMetalFractions", [0] * 10)
    gas_metals = info.get("SubhaloGasMetalFractions", [0] * 10)

    cl = cluster_lookup[grnr]

    return {
        "subhalo_id": sub_id,
        "halo_id": grnr,
        "snap_num": snap_num,
        "redshift": z,
        "M_star_Msun": mstar_msun,
        "R_half_pkpc": rhalf_pkpc,
        "pos_x": pos_x,
        "pos_y": pos_y,
        "pos_z": pos_z,
        "vel_disp_kms": vel_disp,
        "SFR": sfr,
        "sSFR": ssfr,
        "Z_star": z_star,
        "Z_gas": z_gas,
        "Fe_star_frac": star_metals[cfg.ELEM_FE],
        "O_star_frac": star_metals[cfg.ELEM_O],
        "Mg_star_frac": star_metals[cfg.ELEM_MG],
        "Fe_gas_frac": gas_metals[cfg.ELEM_FE],
        "O_gas_frac": gas_metals[cfg.ELEM_O],
        "n_star_particles": meta.get("len_stars", 0),
        "cluster_M200_Msun": cl["M200_Msun"],
        "cluster_R200_pkpc": cl["R200_pkpc"],
    }


def select_ellipticals(snap_num, cluster_catalog, n_workers):
    """
    Two-pass approach:
      Pass 1: Parallel prefilter using meta endpoint (fast).
      Pass 2: Fetch full info.json only for candidates (few hundred).
    """
    z = cfg.SNAPSHOTS[snap_num]

    # Build lookup: halo_id -> cluster row
    cluster_halo_ids = set(cluster_catalog["halo_id"].astype(int))
    cluster_lookup = {
        int(row["halo_id"]): row
        for _, row in cluster_catalog.iterrows()
    }

    # Compute total subhalo range across all selected clusters
    first_sub = int(cluster_catalog["first_subhalo_id"].min())
    last_sub = max(
        int(row["first_subhalo_id"]) + int(row["n_subhalos"])
        for _, row in cluster_catalog.iterrows()
    )
    total = last_sub - first_sub
    print(f"Scanning subhalo IDs {first_sub} to {last_sub - 1} ({total} total)")
    print(f"Using {n_workers} parallel workers")

    # Pass 1: Prefilter
    candidates = []
    done = 0
    lock = threading.Lock()

    print("\n--- Pass 1: Prefiltering (meta endpoint) ---")
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(prefilter_subhalo, snap_num, sid, cluster_halo_ids, z): sid
            for sid in range(first_sub, last_sub)
        }

        for future in as_completed(futures):
            sub_id, passed, meta = future.result()
            done += 1

            if passed:
                candidates.append((sub_id, meta))

            if done % 1000 == 0:
                with lock:
                    print(f"  {done}/{total} checked, "
                          f"{len(candidates)} candidates so far")

    print(f"\nPass 1 complete: {len(candidates)} candidates from {total} subhalos")

    # Pass 2: Fetch full info for candidates
    print(f"\n--- Pass 2: Fetching element fractions ({len(candidates)} galaxies) ---")
    ellipticals = []

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(fetch_full_info, snap_num, sid, meta, z, cluster_lookup): sid
            for sid, meta in candidates
        }

        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                ellipticals.append(result)

    df = pd.DataFrame(ellipticals)
    if len(df) > 0:
        df = df.sort_values("M_star_Msun", ascending=False).reset_index(drop=True)
    print(f"\nTotal ellipticals selected: {len(df)}")
    return df


def main():
    parser = argparse.ArgumentParser(description="Select ellipticals in TNG-Cluster")
    parser.add_argument("--snap", type=int, default=99)
    parser.add_argument("--max-clusters", type=int, default=None,
                        help="Process only the first N clusters (for testing)")
    parser.add_argument("--workers", type=int, default=8,
                        help="Parallel API requests (default: 8)")
    args = parser.parse_args()

    catdir = os.path.join(os.path.dirname(__file__), "..", cfg.CATALOG_DIR)
    cluster_file = os.path.join(catdir, f"cluster_catalog_snap{args.snap:03d}.csv")

    if not os.path.exists(cluster_file):
        print(f"Cluster catalog not found: {cluster_file}")
        print("Run 01_select_clusters.py first.")
        sys.exit(1)

    cluster_cat = pd.read_csv(cluster_file)
    print(f"Loaded {len(cluster_cat)} clusters from {cluster_file}")

    if args.max_clusters:
        cluster_cat = cluster_cat.head(args.max_clusters)
        print(f"Processing first {len(cluster_cat)} clusters only")

    df = select_ellipticals(args.snap, cluster_cat, args.workers)

    outfile = os.path.join(catdir, f"elliptical_catalog_snap{args.snap:03d}.csv")
    df.to_csv(outfile, index=False)
    print(f"Saved to {outfile}")


if __name__ == "__main__":
    main()
