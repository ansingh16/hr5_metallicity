"""
Run missing snapshots one at a time, validating elliptical counts after each.
Retries snapshots that fail due to server errors (corrupt group catalog).
"""
import os
import sys
import subprocess
import time
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

BASEDIR = os.path.join(os.path.dirname(__file__), "..")
CATDIR = os.path.join(BASEDIR, "catalogs")
GCDIR = os.path.join(BASEDIR, "data", "groupcats")

KNOWN_COUNTS = {
    99: (0.0, 4186),
    91: (0.1, 4269),
    84: (0.2, 3992),
    78: (0.3, 2850),
    67: (0.5, 2482),
    59: (0.7, 640),
}

EXPECTED_RANGE = {
    72: (0.4, 1500, 3500),
    50: (1.0, 80, 500),
    40: (1.5, 20, 250),
    33: (2.0, 5, 60),
    25: (2.0, 1, 40),
}

SNAPS_TO_RUN = [50, 40, 33, 25]  # 72 skipped: TNG server chunks 150-160 broken
MAX_SNAP_RETRIES = 5
RETRY_WAIT = 600


def get_elliptical_count(snap):
    path = os.path.join(CATDIR, f"elliptical_catalog_snap{snap:03d}.csv")
    if not os.path.exists(path):
        return None
    return len(pd.read_csv(path))


def clean_snap_artifacts(snap):
    """Remove corrupt catalogs so the pipeline re-downloads them."""
    files_to_remove = [
        os.path.join(GCDIR, f"groupcat_compact_{snap:03d}.hdf5"),
        os.path.join(CATDIR, f"elliptical_catalog_snap{snap:03d}.csv"),
        os.path.join(CATDIR, f"gradients_snap{snap:03d}.csv"),
    ]
    for f in files_to_remove:
        if os.path.exists(f):
            os.remove(f)
            print(f"    Removed {os.path.basename(f)}")


def validate(snap):
    count = get_elliptical_count(snap)
    if count is None:
        print(f"  FAIL: no elliptical catalog produced for snap {snap}")
        return False

    z, lo, hi = EXPECTED_RANGE[snap]
    print(f"  snap {snap} (z={z}): {count} ellipticals")

    all_counts = dict(KNOWN_COUNTS)
    all_counts[snap] = (z, count)
    sorted_snaps = sorted(all_counts.items(), key=lambda x: x[1][0])
    print(f"  Full trend:")
    for s, (zz, n) in sorted_snaps:
        marker = " <--" if s == snap else ""
        print(f"    z={zz:.1f} (snap {s:3d}): {n:>5}{marker}")

    if lo <= count <= hi:
        print(f"  PASS: {count} is within expected range [{lo}, {hi}]")
        KNOWN_COUNTS[snap] = (z, count)
        return True
    elif count > 0:
        print(f"  WARNING: {count} outside expected [{lo}, {hi}] — "
              f"check if physically reasonable given neighbors")
        KNOWN_COUNTS[snap] = (z, count)
        return True
    else:
        print(f"  FAIL: 0 ellipticals — likely corrupt group catalog")
        return False


def run_snap(snap):
    z, lo, hi = EXPECTED_RANGE[snap]

    for attempt in range(MAX_SNAP_RETRIES):
        label = f" (attempt {attempt+1}/{MAX_SNAP_RETRIES})" if attempt > 0 else ""
        print(f"\n{'='*60}")
        print(f"Running snap {snap} (z={z}), expecting {lo}-{hi} ellipticals{label}")
        print(f"{'='*60}")

        if attempt > 0:
            print(f"  Cleaning up artifacts from previous failed attempt...")
            clean_snap_artifacts(snap)

        cmd = [
            sys.executable, "-u",
            os.path.join(os.path.dirname(__file__), "05_redshift_evolution.py"),
            "--snaps", str(snap),
        ]
        result = subprocess.run(cmd, cwd=BASEDIR)

        if result.returncode != 0:
            print(f"  Pipeline returned exit code {result.returncode}")
            if attempt < MAX_SNAP_RETRIES - 1:
                print(f"  Server likely unstable. Waiting {RETRY_WAIT}s before retry...")
                time.sleep(RETRY_WAIT)
                continue
            else:
                print(f"  All {MAX_SNAP_RETRIES} attempts exhausted for snap {snap}.")
                return False

        ok = validate(snap)
        if ok:
            return True

        if attempt < MAX_SNAP_RETRIES - 1:
            print(f"  Validation failed. Waiting {RETRY_WAIT}s before retry...")
            time.sleep(RETRY_WAIT)
        else:
            print(f"  All {MAX_SNAP_RETRIES} attempts exhausted for snap {snap}.")

    return False


if __name__ == "__main__":
    print("Current elliptical counts:")
    for s in sorted(KNOWN_COUNTS, key=lambda x: KNOWN_COUNTS[x][0]):
        z, n = KNOWN_COUNTS[s]
        print(f"  snap {s:3d} (z={z:.1f}): {n}")

    results = {}
    for snap in SNAPS_TO_RUN:
        ok = run_snap(snap)
        results[snap] = ok
        if not ok:
            print(f"\n  Snap {snap} failed after all retries. Stopping.")
            break

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for snap, ok in results.items():
        z = EXPECTED_RANGE[snap][0]
        count = get_elliptical_count(snap)
        status = "OK" if ok else "FAILED"
        print(f"  snap {snap} (z={z}): {count or 0} ellipticals — {status}")

    remaining = [s for s in SNAPS_TO_RUN if s not in results]
    if remaining:
        print(f"  Not run: {remaining}")
