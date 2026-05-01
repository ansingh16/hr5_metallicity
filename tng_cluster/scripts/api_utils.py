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
