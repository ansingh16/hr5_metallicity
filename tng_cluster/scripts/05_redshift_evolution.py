"""
Step 5: Run the full pipeline across multiple snapshots to study
the redshift evolution of metallicity profiles and gradients.

For each snapshot:
  1. Loads the z=0 (snap 99) elliptical catalog as the reference sample.
  2. Traces progenitors via SubLink merger tree (parallel API calls).
  3. Downloads particle data at earlier snapshots (parallel downloads).
  4. Computes profiles (parallel multiprocessing).

Output: catalogs/gradients_snap{NN}.csv    (per-snapshot)
        catalogs/gradients_evolution.csv   (combined)
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import h5py
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg
from scripts.api_utils import get, get_subhalo_info, download_cutout
from scripts.compute_profiles_lib import process_galaxy, fit_gradient


# ── Phase 1: Trace progenitors (I/O bound, threaded) ────────────────────

def trace_one_progenitor(sub_id_z0, target_snap):
    """
    Query the merger tree (main progenitor branch) for one galaxy.
    Uses the simple.json endpoint which returns flat [snap, subfindID] pairs.
    Returns (sub_id_z0, prog_id_or_None).
    """
    url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/99/"
           f"subhalos/{sub_id_z0}/sublink/simple.json")
    try:
        tree = get(url)
    except Exception:
        return sub_id_z0, None

    # simple.json returns {"Main": [[snap, subfindID], ...], "Mergers": [...]}
    main_branch = tree.get("Main", [])
    for snap_num, subfind_id in main_branch:
        if snap_num == target_snap:
            return sub_id_z0, subfind_id

    return sub_id_z0, None


def trace_all_progenitors(ref_df, target_snap, n_workers=8):
    """Trace progenitors for all reference galaxies in parallel."""
    sub_ids_z0 = ref_df["subhalo_id"].astype(int).tolist()
    results = {}
    done = 0

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(trace_one_progenitor, sid, target_snap): sid
            for sid in sub_ids_z0
        }
        for future in as_completed(futures):
            sub_id_z0, prog_id = future.result()
            if prog_id is not None:
                results[sub_id_z0] = prog_id
            done += 1
            if done % 100 == 0:
                print(f"    Traced {done}/{len(sub_ids_z0)}, "
                      f"{len(results)} found")

    return results


# ── Phase 2: Get progenitor properties + filter (I/O, threaded) ─────────

def get_progenitor_props(prog_id, target_snap, z):
    """Fetch meta for a progenitor. Returns (prog_id, props_dict_or_None)."""
    try:
        meta = get_subhalo_info(target_snap, prog_id)
        if not isinstance(meta, dict):
            # API returned non-JSON (HTML page, etc.) — skip this progenitor
            return prog_id, None
    except Exception:
        return prog_id, None

    # info.json returns full catalog field names; stellar = index 4
    mass_type = meta.get("SubhaloMassType", [0]*6)
    len_type = meta.get("SubhaloLenType", [0]*6)
    rad_type = meta.get("SubhaloHalfmassRadType", [0]*6)

    mstar_msun = mass_type[4] * 1e10 / cfg.h
    if mstar_msun < cfg.MIN_STELLAR_MASS_MSUN:
        return prog_id, None

    nstars = len_type[4]
    if nstars < cfg.MIN_STAR_PARTICLES:
        return prog_id, None

    a = 1.0 / (1.0 + z)
    rhalf_pkpc = rad_type[4] * a / cfg.h
    if rhalf_pkpc < 0.5:
        return prog_id, None

    return prog_id, {
        "M_star_Msun": mstar_msun,
        "R_half_pkpc": rhalf_pkpc,
        "n_stars": nstars,
    }


def filter_progenitors(progenitor_map, target_snap, z, n_workers=8):
    """Fetch properties and filter progenitors in parallel."""
    valid = {}
    done = 0

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(get_progenitor_props, pid, target_snap, z): (sid_z0, pid)
            for sid_z0, pid in progenitor_map.items()
        }
        for future in as_completed(futures):
            sid_z0, pid = futures[future]
            prog_id, props = future.result()
            if props is not None:
                valid[sid_z0] = (prog_id, props)
            done += 1
            if done % 100 == 0:
                print(f"    Checked {done}/{len(progenitor_map)}, "
                      f"{len(valid)} pass cuts")

    return valid


# ── Phase 3: Download cutouts (I/O, threaded) ───────────────────────────

def download_one_progenitor(prog_id, target_snap, datadir):
    """Download particle cutout for one progenitor. Returns (prog_id, filepath_or_None)."""
    filepath = os.path.join(datadir, f"subhalo_{prog_id}.hdf5")
    if os.path.exists(filepath):
        return prog_id, filepath

    star_fields = ("Coordinates,Masses,GFM_InitialMass,"
                   "GFM_Metallicity,GFM_Metals,"
                   "GFM_StellarFormationTime,Velocities,ParticleIDs")
    gas_fields = ("Coordinates,Masses,GFM_Metallicity,"
                  "GFM_Metals,StarFormationRate,Velocities,"
                  "Density,InternalEnergy,ElectronAbundance,ParticleIDs")
    try:
        fn = download_cutout(
            target_snap, prog_id,
            {"stars": star_fields, "gas": gas_fields},
            save_dir=datadir
        )
        if fn != filepath:
            os.rename(fn, filepath)
        return prog_id, filepath
    except Exception as e:
        print(f"    Download failed for {prog_id}: {e}")
        return prog_id, None


def download_all(valid_progenitors, target_snap, datadir, n_workers=4):
    """Download cutouts for all valid progenitors in parallel."""
    paths = {}
    done = 0
    total = len(valid_progenitors)

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(download_one_progenitor, prog_id, target_snap, datadir): sid_z0
            for sid_z0, (prog_id, _) in valid_progenitors.items()
        }
        for future in as_completed(futures):
            sid_z0 = futures[future]
            prog_id, filepath = future.result()
            if filepath is not None:
                paths[sid_z0] = filepath
            done += 1
            if done % 50 == 0:
                print(f"    Downloaded {done}/{total}, {len(paths)} success")

    return paths


# ── Phase 4: Compute profiles (CPU bound, multiprocessing) ──────────────

def compute_one_profile(filepath, prog_id, r_half_pkpc, z,
                        sub_id_z0, halo_id, mstar, snap_num):
    """Compute profiles and gradients for one progenitor galaxy."""
    result = process_galaxy(filepath, prog_id, r_half_pkpc, z)
    if result is None:
        return None

    bins = result["r_bins"]
    grad_row = {
        "subhalo_id": prog_id,
        "subhalo_id_z0": sub_id_z0,
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
                bins, result[key], cfg.PROFILE_FIT_INNER_MIN, cfg.PROFILE_FIT_RBREAK)
            grad_row[f"slope_inner_{key}"] = si
            grad_row[f"slope_err_inner_{key}"] = sei

            so, ico, seo = fit_gradient(
                bins, result[key], cfg.PROFILE_FIT_RBREAK, cfg.PROFILE_FIT_OUTER_MAX)
            grad_row[f"slope_outer_{key}"] = so
            grad_row[f"slope_err_outer_{key}"] = seo

    if result.get("has_gas", False):
        for key in ["FeH_gas", "Z_gas", "SFR_profile"]:
            if key in result:
                s, ic, se = fit_gradient(
                    bins, result[key], cfg.PROFILE_FIT_RMIN, cfg.PROFILE_FIT_RMAX)
                grad_row[f"slope_{key}"] = s
                grad_row[f"intercept_{key}"] = ic
                grad_row[f"slope_err_{key}"] = se

    for gkey in ["FeH_star_global", "OFe_star_global", "MgFe_star_global",
                 "NO_star_global", "CO_star_global", "CN_star_global",
                 "age_star_global"]:
        if gkey in result:
            grad_row[gkey] = result[gkey]

    return grad_row


def compute_all_profiles(valid_progenitors, paths, ref_df, target_snap, z,
                         n_workers=None):
    """Compute profiles for all downloaded progenitors in parallel."""
    if n_workers is None:
        n_workers = min(os.cpu_count() or 4, 12)

    # Build ref lookup
    ref_lookup = ref_df.set_index("subhalo_id").to_dict("index")

    jobs = []
    for sid_z0, filepath in paths.items():
        prog_id, props = valid_progenitors[sid_z0]
        ref_row = ref_lookup.get(sid_z0, {})
        halo_id = int(ref_row.get("halo_id", -1))
        jobs.append((filepath, prog_id, props["R_half_pkpc"], z,
                      sid_z0, halo_id, props["M_star_Msun"], target_snap))

    gradient_rows = []
    done = 0

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(compute_one_profile, *job): job[4] for job in jobs}
        for future in as_completed(futures):
            done += 1
            try:
                grad_row = future.result()
                if grad_row is not None:
                    gradient_rows.append(grad_row)
            except Exception as e:
                print(f"    Profile error: {e}")

            if done % 50 == 0:
                print(f"    Profiles: {done}/{len(jobs)}, "
                      f"{len(gradient_rows)} successful")

    return gradient_rows


# ── Main pipeline ────────────────────────────────────────────────────────

def process_snapshot(ref_df, target_snap, basedir, catdir,
                     api_workers=8, dl_workers=4, cpu_workers=None):
    """Full pipeline for one snapshot."""
    z = cfg.SNAPSHOTS[target_snap]
    datadir = os.path.join(basedir, cfg.DATA_DIR, f"snap{target_snap:03d}")
    os.makedirs(datadir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Snapshot {target_snap} (z={z:.2f})")
    print(f"{'='*60}")

    # Phase 1: Trace progenitors
    print(f"\n  Phase 1: Tracing progenitors ({len(ref_df)} galaxies)...")
    progenitor_map = trace_all_progenitors(ref_df, target_snap, api_workers)
    print(f"  Found {len(progenitor_map)} progenitors")

    if not progenitor_map:
        print("  No progenitors found, skipping")
        return []

    # Phase 2: Filter by mass/size
    print(f"\n  Phase 2: Filtering progenitors...")
    valid = filter_progenitors(progenitor_map, target_snap, z, api_workers)
    print(f"  {len(valid)} pass mass/size cuts")

    if not valid:
        print("  No valid progenitors, skipping")
        return []

    # Phase 3: Download cutouts
    print(f"\n  Phase 3: Downloading particle data...")
    paths = download_all(valid, target_snap, datadir, dl_workers)
    print(f"  {len(paths)} cutouts ready")

    if not paths:
        print("  No cutouts available, skipping")
        return []

    # Phase 4: Compute profiles
    print(f"\n  Phase 4: Computing profiles...")
    gradient_rows = compute_all_profiles(
        valid, paths, ref_df, target_snap, z, cpu_workers)
    print(f"  {len(gradient_rows)} profiles computed")

    # Save per-snapshot results
    if gradient_rows:
        snap_grad_df = pd.DataFrame(gradient_rows)
        snap_grad_file = os.path.join(catdir, f"gradients_snap{target_snap:03d}.csv")
        snap_grad_df.to_csv(snap_grad_file, index=False)
        print(f"  Saved to {snap_grad_file}")

    return gradient_rows


def main():
    parser = argparse.ArgumentParser(
        description="Compute redshift evolution of metallicity profiles"
    )
    parser.add_argument("--snaps", type=int, nargs="+",
                        default=[91, 84, 78, 67, 50, 33, 25],
                        help="Snapshot numbers to process")
    parser.add_argument("--api-workers", type=int, default=8,
                        help="Parallel API requests (default: 8)")
    parser.add_argument("--dl-workers", type=int, default=4,
                        help="Parallel downloads (default: 4)")
    parser.add_argument("--cpu-workers", type=int, default=None,
                        help="Parallel profile computations (default: CPU count)")
    args = parser.parse_args()

    for s in args.snaps:
        if s not in cfg.SNAPSHOTS:
            print(f"Unknown snapshot {s}. Available: {list(cfg.SNAPSHOTS.keys())}")
            sys.exit(1)

    basedir = os.path.join(os.path.dirname(__file__), "..")
    catdir = os.path.join(basedir, cfg.CATALOG_DIR)

    ref_file = os.path.join(catdir, "elliptical_catalog_snap099.csv")
    if not os.path.exists(ref_file):
        print("Reference catalog (snap099) not found. Run steps 01-04 first.")
        sys.exit(1)

    ref_df = pd.read_csv(ref_file)
    print(f"Reference sample: {len(ref_df)} ellipticals at z=0")

    all_gradients = []

    # Include z=0 gradients
    grad_z0 = os.path.join(catdir, "gradients_snap099.csv")
    if os.path.exists(grad_z0):
        g0 = pd.read_csv(grad_z0)
        all_gradients.append(g0)
        print(f"Loaded {len(g0)} gradient measurements at z=0")

    for target_snap in args.snaps:
        rows = process_snapshot(
            ref_df, target_snap, basedir, catdir,
            args.api_workers, args.dl_workers, args.cpu_workers
        )
        if rows:
            all_gradients.append(pd.DataFrame(rows))

    # Combined evolution table
    if all_gradients:
        combined = pd.concat(all_gradients, ignore_index=True)
        combined_file = os.path.join(catdir, "gradients_evolution.csv")
        combined.to_csv(combined_file, index=False)
        print(f"\nCombined evolution: {len(combined)} rows -> {combined_file}")


if __name__ == "__main__":
    main()
