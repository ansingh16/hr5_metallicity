"""
Utility functions for interacting with the TNG API.
"""

import os
import sys
import time
import requests
import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

HEADERS = {"api-key": cfg.API_KEY}


def get(path, params=None, retries=5, backoff=10):
    """
    Generic GET request to the TNG API with retry logic.
    Returns JSON dict, or saves and returns filename for HDF5 cutouts.
    """
    for attempt in range(retries):
        try:
            r = requests.get(path, params=params, headers=HEADERS, timeout=300)
            r.raise_for_status()

            ct = r.headers.get("content-type", "")
            if "application/json" in ct:
                return r.json()

            if "content-disposition" in r.headers:
                filename = r.headers["content-disposition"].split("filename=")[1]
                with open(filename, "wb") as f:
                    f.write(r.content)
                return filename

            return r

        except (requests.exceptions.RequestException, requests.exceptions.Timeout) as e:
            if attempt < retries - 1:
                wait = backoff * (2 ** attempt)
                print(f"  Request failed ({e}), retrying in {wait}s...")
                time.sleep(wait)
            else:
                raise


def get_halo_info(snap_num, halo_id):
    """
    Return full FoF group catalog properties for one halo.
    Uses the info.json endpoint which returns all GroupCat fields.
    """
    url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/"
           f"halos/{halo_id}/info.json")
    return get(url)


def get_halo_meta(snap_num, halo_id):
    """
    Return halo metadata (child_subhalos list, vis links, etc).
    This is the default halo endpoint (not info.json).
    """
    url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/"
           f"halos/{halo_id}/")
    return get(url)


def get_subhalo_info(snap_num, subhalo_id):
    """Return full SubhaloFields for one subhalo via info.json."""
    url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/"
           f"subhalos/{subhalo_id}/info.json")
    return get(url)


def get_subhalo_meta(snap_num, subhalo_id):
    """Return subhalo metadata (links, basic info)."""
    url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/"
           f"subhalos/{subhalo_id}/")
    return get(url)


def search_subhalos(snap_num, params):
    """
    Search subhalos with filter parameters. Returns paginated JSON.
    Example: search_subhalos(99, {'mass_stars__gt': 0.1, 'limit': 100,
                                   'order_by': '-mass_stars'})
    """
    url = f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/subhalos/"
    return get(url, params=params)


def get_snapshot_info(snap_num):
    """Return metadata dict for a specific snapshot."""
    return get(f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/")


def download_cutout(snap_num, subhalo_id, particle_fields, save_dir=None):
    """
    Download particle cutout HDF5 for a subhalo.
    Thread-safe: writes directly to save_dir without os.chdir.

    Parameters
    ----------
    snap_num : int
    subhalo_id : int
    particle_fields : dict
        e.g. {'stars': 'Coordinates,Masses,GFM_Metals,...',
              'gas': 'Coordinates,Masses,GFM_Metals,...'}
    save_dir : str or None

    Returns
    -------
    str : path to the downloaded HDF5 file.
    """
    url = (
        f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/"
        f"subhalos/{subhalo_id}/cutout.hdf5"
    )

    r = requests.get(url, params=particle_fields, headers=HEADERS, timeout=300)
    r.raise_for_status()

    # Extract filename from Content-Disposition header
    filename = r.headers["content-disposition"].split("filename=")[1]

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        filepath = os.path.join(save_dir, filename)
    else:
        filepath = filename

    with open(filepath, "wb") as f:
        f.write(r.content)

    return filepath


def paginate_subhalos(snap_num, params, max_items=None):
    """
    Iterate through paginated subhalo search results.
    Yields individual result dicts with keys: id, sfr, mass_log_msun, url.
    """
    url = f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/subhalos/"
    count = 0
    while url is not None:
        page = get(url, params=params)
        params = None  # only needed for first request
        for item in page["results"]:
            yield item
            count += 1
            if max_items and count >= max_items:
                return
        url = page.get("next")


def get_halo_infos_parallel(snap_num, halo_ids, n_workers=16):
    """Fetch halo info.json for multiple halo IDs in parallel."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    results = {}
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(get_halo_info, snap_num, hid): hid
            for hid in halo_ids
        }
        for future in as_completed(futures):
            hid = futures[future]
            try:
                results[hid] = future.result()
            except Exception as e:
                print(f"  Halo {hid}: failed ({e})")
    return results


def search_subhalos_in_halo(snap_num, halo_id, min_len_stars=20):
    """Search for subhalos belonging to a specific FoF halo via grnr filter."""
    params = {
        "grnr": halo_id,
        "len_stars__gt": min_len_stars - 1,
        "limit": 50000,
        "order_by": "-mass_stars",
    }
    return list(paginate_subhalos(snap_num, params))


def iterate_child_subhalos(snap_num, halo_id, max_items=None):
    """
    Iterate through child subhalos of a FoF halo.
    The halo endpoint returns paginated child_subhalos.
    Yields dicts with keys: id, url.
    """
    url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/snapshots/{snap_num}/"
           f"halos/{halo_id}/")
    count = 0
    while url is not None:
        page = get(url)
        children = page.get("child_subhalos", {})
        for item in children.get("results", []):
            yield item
            count += 1
            if max_items and count >= max_items:
                return
        url = children.get("next")


# ── Group catalog bulk download ──────────────────────────────────────────

GROUPCAT_SUBHALO_FIELDS = [
    "SubhaloGrNr", "SubhaloMassType", "SubhaloLenType", "SubhaloFlag",
    "SubhaloHalfmassRadType", "SubhaloPos", "SubhaloSFR", "SubhaloVelDisp",
    "SubhaloStarMetallicity", "SubhaloGasMetallicity",
    "SubhaloStarMetalFractions", "SubhaloGasMetalFractions",
]

GROUPCAT_GROUP_FIELDS = [
    "GroupPrimaryZoomTarget", "Group_M_Crit200", "Group_R_Crit200",
    "GroupFirstSub", "GroupMassType", "GroupNsubs", "GroupPos",
]


def download_group_catalog(snap_num, save_dir, fields_sub=None, fields_grp=None):
    """
    Download TNG-Cluster group catalog chunks and extract needed fields
    into a single compact HDF5 file. Processes chunks one at a time to
    minimise disk usage. Skips download if compact catalog already exists.

    Returns path to the compact catalog HDF5.
    """
    if fields_sub is None:
        fields_sub = GROUPCAT_SUBHALO_FIELDS
    if fields_grp is None:
        fields_grp = GROUPCAT_GROUP_FIELDS

    os.makedirs(save_dir, exist_ok=True)
    compact_path = os.path.join(save_dir, f"groupcat_compact_{snap_num:03d}.hdf5")

    if os.path.exists(compact_path):
        print(f"  Group catalog already cached: {compact_path}")
        return compact_path

    listing_url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/files/groupcat-{snap_num}/")
    listing = get(listing_url)
    chunk_urls = listing["files"]
    n_chunks = len(chunk_urls)
    print(f"  Downloading group catalog: {n_chunks} chunks for snap {snap_num}")

    sub_arrays = {f: [] for f in fields_sub}
    grp_arrays = {f: [] for f in fields_grp}
    failed_chunks = []

    for i, url in enumerate(chunk_urls):
        tmp_path = os.path.join(save_dir, f"_chunk_{snap_num}_{i}.hdf5")
        max_retries = 8
        backoff = 60
        success = False

        for attempt in range(max_retries):
            try:
                r = requests.get(url, headers=HEADERS, timeout=600, stream=True)
                r.raise_for_status()
                with open(tmp_path, "wb") as f:
                    for block in r.iter_content(chunk_size=1024 * 1024):
                        f.write(block)

                with h5py.File(tmp_path, "r") as hf:
                    if "Subhalo" in hf:
                        for field in fields_sub:
                            if field in hf["Subhalo"]:
                                sub_arrays[field].append(hf["Subhalo"][field][:])
                    if "Group" in hf:
                        for field in fields_grp:
                            if field in hf["Group"]:
                                grp_arrays[field].append(hf["Group"][field][:])
                success = True
                break
            except Exception as e:
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass
                if attempt < max_retries - 1:
                    wait = backoff * (2 ** attempt)
                    print(f"    Chunk {i}/{n_chunks} failed ({e}), retry in {wait}s...")
                    time.sleep(wait)
                else:
                    print(f"    Chunk {i}/{n_chunks} FAILED after {max_retries} attempts: {e}")

        if success:
            try:
                os.remove(tmp_path)
            except OSError:
                pass
            time.sleep(2)
        else:
            failed_chunks.append(i)

        if (i + 1) % 50 == 0 or i == n_chunks - 1:
            print(f"    {i+1}/{n_chunks} chunks processed")

    if failed_chunks:
        print(f"  FATAL: {len(failed_chunks)} group catalog chunks failed: {failed_chunks}")
        print(f"  Cannot build compact catalog — SubhaloGrNr indexing would be broken.")
        raise RuntimeError(
            f"Group catalog download incomplete: {len(failed_chunks)}/{n_chunks} "
            f"chunks failed for snap {snap_num}. Retry later."
        )

    with h5py.File(compact_path, "w") as hf:
        sg = hf.create_group("Subhalo")
        for field in fields_sub:
            if sub_arrays[field]:
                sg.create_dataset(field, data=np.concatenate(sub_arrays[field]))
        gg = hf.create_group("Group")
        for field in fields_grp:
            if grp_arrays[field]:
                gg.create_dataset(field, data=np.concatenate(grp_arrays[field]))

    n_sub = len(sub_arrays[fields_sub[0]]) if sub_arrays[fields_sub[0]] else 0
    total_sub = sum(len(a) for a in sub_arrays[fields_sub[0]]) if sub_arrays[fields_sub[0]] else 0
    print(f"  Saved compact catalog: {total_sub} subhalos -> {compact_path}")
    return compact_path


def load_group_catalog(compact_path):
    """Load compact group catalog HDF5 into dicts of numpy arrays."""
    sub = {}
    grp = {}
    with h5py.File(compact_path, "r") as hf:
        for field in hf["Subhalo"]:
            sub[field] = hf["Subhalo"][field][:]
        for field in hf["Group"]:
            grp[field] = hf["Group"][field][:]
    return sub, grp


def find_clusters_from_catalog(grp, min_m200_msun=None):
    """
    Identify primary zoom-target clusters from group catalog arrays.
    Returns array of group indices (= halo IDs in the unified catalog).
    """
    if min_m200_msun is None:
        min_m200_msun = cfg.MIN_GROUP_M200_MSUN

    is_zoom = grp["GroupPrimaryZoomTarget"].astype(bool)
    m200_msun = grp["Group_M_Crit200"] * 1e10 / cfg.h
    massive = m200_msun >= min_m200_msun

    cluster_mask = is_zoom & massive
    cluster_ids = np.where(cluster_mask)[0]
    return cluster_ids


def filter_candidates_from_catalog(sub, cluster_ids, z):
    """
    Apply all non-Sérsic selection cuts using local group catalog arrays.
    Returns boolean mask over the subhalo arrays.
    """
    n = len(sub["SubhaloGrNr"])

    in_cluster = np.isin(sub["SubhaloGrNr"], cluster_ids)
    has_flag = sub["SubhaloFlag"].astype(bool)
    nstars = sub["SubhaloLenType"][:, 4]
    enough_stars = nstars >= cfg.MIN_STAR_PARTICLES
    mstar = sub["SubhaloMassType"][:, 4] * 1e10 / cfg.h
    massive_enough = mstar >= cfg.MIN_STELLAR_MASS_MSUN

    total_mass = sub["SubhaloMassType"].sum(axis=1)
    dm_lowres = sub["SubhaloMassType"][:, 2]
    safe = np.ones(n, dtype=bool)
    pos_mass = total_mass > 0
    safe[pos_mass] = (dm_lowres[pos_mass] / total_mass[pos_mass]) <= cfg.MAX_CONTAMINATION_FRAC

    a = 1.0 / (1.0 + z)
    rhalf = sub["SubhaloHalfmassRadType"][:, 4] * a / cfg.h
    big_enough = rhalf >= 0.5

    mask = in_cluster & has_flag & enough_stars & massive_enough & safe & big_enough
    return mask


def download_cutout_safe(snap_num, subhalo_id, particle_fields, save_dir=None,
                         retries=3, backoff=10):
    """download_cutout with retry logic and rate-friendly pacing."""
    for attempt in range(retries):
        try:
            return download_cutout(snap_num, subhalo_id, particle_fields,
                                   save_dir=save_dir)
        except Exception as e:
            if attempt < retries - 1:
                wait = backoff * (2 ** attempt)
                print(f"    Cutout {subhalo_id} failed ({e}), retry in {wait}s...")
                time.sleep(wait)
            else:
                raise


# ── Snapshot bulk download + extraction ────────────────────────────────────

def download_offsets(snap_num, save_dir):
    """Download and cache the particle offsets file for a snapshot."""
    os.makedirs(save_dir, exist_ok=True)
    outpath = os.path.join(save_dir, f"offsets_{snap_num:03d}.hdf5")
    if os.path.exists(outpath):
        print(f"  Offsets already cached: {outpath}")
        return outpath
    url = f"{cfg.BASE_URL}{cfg.SIM_NAME}/files/offsets.{snap_num}.hdf5"
    print(f"  Downloading offsets for snap {snap_num}...")
    r = requests.get(url, headers=HEADERS, timeout=1800, stream=True)
    r.raise_for_status()
    with open(outpath, "wb") as f:
        for block in r.iter_content(1024 * 1024):
            f.write(block)
    size_mb = os.path.getsize(outpath) / 1e6
    print(f"  Saved offsets ({size_mb:.0f} MB)")
    return outpath


def build_chunk_map(offsets_path, candidate_indices, sub_lentype):
    """
    Map candidate subhalos to snapshot chunk files using particle offsets.

    Returns (chunk_map, gas_extra) where:
      chunk_map[file_idx] = list of candidate dicts with star (+ gas if same file)
      gas_extra[file_idx] = list of candidate dicts needing gas from a different file
    """
    with h5py.File(offsets_path, "r") as hf:
        sub_snap = hf["Subhalo/SnapByType"][:]
        file_off = hf["FileOffsets/SnapByType"][:]

    star_globals = sub_snap[candidate_indices, 4]
    gas_globals = sub_snap[candidate_indices, 0]

    star_files = np.searchsorted(file_off[:, 4], star_globals, side='right') - 1
    gas_files = np.searchsorted(file_off[:, 0], gas_globals, side='right') - 1

    star_locals = star_globals - file_off[star_files, 4]
    gas_locals = gas_globals - file_off[gas_files, 0]

    n_stars = sub_lentype[candidate_indices, 4]
    n_gas = sub_lentype[candidate_indices, 0]

    chunk_map = {}
    gas_extra = {}

    for i, idx in enumerate(candidate_indices):
        sf = int(star_files[i])
        gf = int(gas_files[i])

        entry = {
            'idx': int(idx),
            'star_offset': int(star_locals[i]),
            'n_stars': int(n_stars[i]),
        }

        if gf == sf:
            entry['gas_offset'] = int(gas_locals[i])
            entry['n_gas'] = int(n_gas[i])
        else:
            entry['gas_offset'] = None
            entry['n_gas'] = 0
            if int(n_gas[i]) > 0:
                gas_extra.setdefault(gf, []).append({
                    'idx': int(idx),
                    'gas_offset': int(gas_locals[i]),
                    'n_gas': int(n_gas[i]),
                })

        chunk_map.setdefault(sf, []).append(entry)

    return chunk_map, gas_extra


def download_snapshot_chunk(snap_num, file_idx, save_dir,
                            retries=3, backoff=30):
    """Download a single snapshot chunk file with retry logic."""
    os.makedirs(save_dir, exist_ok=True)
    outpath = os.path.join(save_dir, f"snap_{snap_num:03d}.{file_idx}.hdf5")
    url = (f"{cfg.BASE_URL}{cfg.SIM_NAME}/files/"
           f"snapshot-{snap_num}.{file_idx}.hdf5")
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=HEADERS, timeout=1800, stream=True)
            r.raise_for_status()
            with open(outpath, "wb") as f:
                for block in r.iter_content(1024 * 1024):
                    f.write(block)
            return outpath
        except Exception as e:
            if attempt < retries - 1:
                wait = backoff * (2 ** attempt)
                print(f"    Chunk {file_idx} failed ({e}), retry in {wait}s...")
                time.sleep(wait)
            else:
                raise


def select_and_extract_bulk(snap_num, gcdir, particle_save_dir, sersic_func,
                            dl_workers=1):
    """
    Select ellipticals and extract full particle data via bulk snapshot files.

    Downloads ~352 snapshot chunks (one per zoom region) instead of ~64,000
    individual cutout API calls. For each chunk, computes Sersic indices and
    saves full star+gas particle data for ellipticals.

    Parameters
    ----------
    snap_num : int
    gcdir : str
        Directory for group catalog and offsets cache.
    particle_save_dir : str
        Directory to save subhalo_*.hdf5 particle files for ellipticals.
    sersic_func : callable
        compute_sersic_index(coords_3d, masses, r_half) -> float or NaN.

    Returns
    -------
    pd.DataFrame with elliptical catalog columns.
    """
    import tempfile
    import pandas as pd
    from concurrent.futures import ThreadPoolExecutor, as_completed

    z = cfg.SNAPSHOTS[snap_num]
    a = 1.0 / (1.0 + z)

    compact_path = download_group_catalog(snap_num, gcdir)
    sub, grp = load_group_catalog(compact_path)

    cluster_ids = find_clusters_from_catalog(grp)
    print(f"  {len(cluster_ids)} primary zoom-target clusters")

    if len(cluster_ids) == 0:
        return pd.DataFrame()

    cluster_lookup = {}
    for cid in cluster_ids:
        cluster_lookup[int(cid)] = {
            "M200_Msun": float(grp["Group_M_Crit200"][cid] * 1e10 / cfg.h),
            "R200_pkpc": float(grp["Group_R_Crit200"][cid] * a / cfg.h),
        }

    candidate_mask = filter_candidates_from_catalog(sub, cluster_ids, z)
    candidate_indices = np.where(candidate_mask)[0]
    n_candidates = len(candidate_indices)
    print(f"  {n_candidates} candidates after local filtering")

    if n_candidates == 0:
        return pd.DataFrame()

    offsets_path = download_offsets(snap_num, gcdir)
    chunk_map, gas_extra = build_chunk_map(
        offsets_path, candidate_indices, sub["SubhaloLenType"])
    n_chunks = len(chunk_map)
    print(f"  Mapped to {n_chunks} snapshot chunks"
          + (f" (+{len(gas_extra)} gas-only)" if gas_extra else ""))

    os.makedirs(particle_save_dir, exist_ok=True)
    tmpdir = tempfile.mkdtemp(prefix="tng_snap_")

    ellipticals = []
    n_processed = 0
    n_sersic_fail = 0
    elliptical_ids = set()

    sorted_chunks = sorted(chunk_map.keys())

    with ThreadPoolExecutor(max_workers=dl_workers) as dl_pool:
        dl_futures = {
            dl_pool.submit(
                download_snapshot_chunk, snap_num, fi, tmpdir): fi
            for fi in sorted_chunks
        }

        for ci, future in enumerate(as_completed(dl_futures)):
            file_idx = dl_futures[future]
            cands = chunk_map[file_idx]
            chunk_path = None

            try:
                chunk_path = future.result()
            except Exception as e:
                print(f"    Chunk {file_idx} download failed: {e}")
                n_processed += len(cands)
                if (ci + 1) % 10 == 0:
                    print(f"    [{ci+1}/{n_chunks}] "
                          f"{n_processed}/{n_candidates} processed, "
                          f"{len(ellipticals)} ellipticals")
                continue

            chunk_mb = os.path.getsize(chunk_path) / 1e6

            try:
                with h5py.File(chunk_path, "r") as hf:
                    has_pt4 = "PartType4" in hf
                    has_pt0 = "PartType0" in hf

                    for cand in cands:
                        idx = cand['idx']
                        n_processed += 1

                        if not has_pt4 or cand['n_stars'] == 0:
                            continue

                        s0 = cand['star_offset']
                        s1 = s0 + cand['n_stars']

                        coords = hf["PartType4"]["Coordinates"][s0:s1]
                        masses = hf["PartType4"]["Masses"][s0:s1]
                        a_form = hf["PartType4"]["GFM_StellarFormationTime"][s0:s1]

                        star_mask = a_form > 0
                        if star_mask.sum() < cfg.MIN_STAR_PARTICLES:
                            continue

                        rhalf = float(
                            sub["SubhaloHalfmassRadType"][idx, 4] * a / cfg.h)

                        try:
                            sersic_n = sersic_func(
                                coords[star_mask] * a / cfg.h,
                                masses[star_mask] * 1e10 / cfg.h,
                                rhalf)
                        except Exception:
                            sersic_n = np.nan

                        if not np.isfinite(sersic_n):
                            n_sersic_fail += 1
                            continue
                        if sersic_n <= cfg.SERSIC_ELLIPTICAL_CUT:
                            continue

                        outpath = os.path.join(
                            particle_save_dir, f"subhalo_{idx}.hdf5")
                        with h5py.File(outpath, "w") as out:
                            sg = out.create_group("PartType4")
                            for field in hf["PartType4"]:
                                sg.create_dataset(
                                    field,
                                    data=hf["PartType4"][field][s0:s1])

                            if (has_pt0
                                    and cand.get('gas_offset') is not None
                                    and cand['n_gas'] > 0):
                                g0 = cand['gas_offset']
                                g1 = g0 + cand['n_gas']
                                gg = out.create_group("PartType0")
                                for field in hf["PartType0"]:
                                    gg.create_dataset(
                                        field,
                                        data=hf["PartType0"][field][g0:g1])

                        mstar = float(
                            sub["SubhaloMassType"][idx, 4] * 1e10 / cfg.h)
                        pos = sub["SubhaloPos"][idx]
                        halo_id = int(sub["SubhaloGrNr"][idx])
                        sfr = float(sub["SubhaloSFR"][idx])
                        cl = cluster_lookup.get(
                            halo_id, {"M200_Msun": 0, "R200_pkpc": 0})

                        ellipticals.append({
                            "subhalo_id": idx,
                            "halo_id": halo_id,
                            "snap_num": snap_num,
                            "redshift": z,
                            "M_star_Msun": mstar,
                            "R_half_pkpc": rhalf,
                            "sersic_n": sersic_n,
                            "pos_x": float(pos[0]),
                            "pos_y": float(pos[1]),
                            "pos_z": float(pos[2]),
                            "vel_disp_kms": float(sub["SubhaloVelDisp"][idx]),
                            "SFR": sfr,
                            "sSFR": sfr / mstar if mstar > 0 else 0.0,
                            "Z_star": float(
                                sub["SubhaloStarMetallicity"][idx]),
                            "Z_gas": float(
                                sub["SubhaloGasMetallicity"][idx]),
                            "Fe_star_frac": float(
                                sub["SubhaloStarMetalFractions"][
                                    idx, cfg.ELEM_FE]),
                            "O_star_frac": float(
                                sub["SubhaloStarMetalFractions"][
                                    idx, cfg.ELEM_O]),
                            "Mg_star_frac": float(
                                sub["SubhaloStarMetalFractions"][
                                    idx, cfg.ELEM_MG]),
                            "Fe_gas_frac": float(
                                sub["SubhaloGasMetalFractions"][
                                    idx, cfg.ELEM_FE]),
                            "O_gas_frac": float(
                                sub["SubhaloGasMetalFractions"][
                                    idx, cfg.ELEM_O]),
                            "n_star_particles": int(
                                sub["SubhaloLenType"][idx, 4]),
                            "cluster_M200_Msun": cl["M200_Msun"],
                            "cluster_R200_pkpc": cl["R200_pkpc"],
                        })
                        elliptical_ids.add(idx)
            except Exception as e:
                print(f"    Error processing chunk {file_idx}: {e}")
            finally:
                if chunk_path:
                    try:
                        os.remove(chunk_path)
                    except OSError:
                        pass

            if (ci + 1) % 10 == 0 or ci == len(sorted_chunks) - 1:
                print(f"    [{ci+1}/{n_chunks}] {chunk_mb:.0f} MB | "
                      f"{n_processed}/{n_candidates} processed, "
                      f"{len(ellipticals)} ellipticals")

    # Second pass: gas data from different files (rare)
    if gas_extra and elliptical_ids:
        needed = {fid: [c for c in cs if c['idx'] in elliptical_ids]
                  for fid, cs in gas_extra.items()}
        needed = {fid: cs for fid, cs in needed.items() if cs}
        if needed:
            print(f"  Fetching gas from {len(needed)} additional chunks...")
            for file_idx, cands in sorted(needed.items()):
                try:
                    chunk_path = download_snapshot_chunk(
                        snap_num, file_idx, tmpdir)
                    with h5py.File(chunk_path, "r") as hf:
                        if "PartType0" in hf:
                            for cand in cands:
                                g0 = cand['gas_offset']
                                g1 = g0 + cand['n_gas']
                                fpath = os.path.join(
                                    particle_save_dir,
                                    f"subhalo_{cand['idx']}.hdf5")
                                with h5py.File(fpath, "a") as out:
                                    if "PartType0" not in out:
                                        gg = out.create_group("PartType0")
                                        for field in hf["PartType0"]:
                                            gg.create_dataset(
                                                field,
                                                data=hf["PartType0"][field][
                                                    g0:g1])
                except Exception as e:
                    print(f"    Gas chunk {file_idx} failed: {e}")
                finally:
                    try:
                        os.remove(chunk_path)
                    except OSError:
                        pass

    try:
        os.rmdir(tmpdir)
    except OSError:
        pass

    df = pd.DataFrame(ellipticals)
    if len(df) > 0:
        df = df.sort_values(
            "M_star_Msun", ascending=False).reset_index(drop=True)
    print(f"\n  Ellipticals: {len(df)} "
          f"(Sersic fail: {n_sersic_fail}, "
          f"download fail: {n_candidates - n_processed})")
    return df
