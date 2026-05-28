"""
Check completeness of the redshift evolution pipeline.
Flags snapshots with suspiciously low counts that may need re-running.
"""

import os
import sys
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

BASEDIR = os.path.join(os.path.dirname(__file__), "..")
CATDIR = os.path.join(BASEDIR, cfg.CATALOG_DIR)
DATADIR = os.path.join(BASEDIR, cfg.DATA_DIR)

SNAPS = [99, 91, 84, 78, 72, 67, 59, 50, 40, 33, 25]

# Minimum expected ellipticals per snapshot (rough lower bound)
# Based on snaps 84/78 having ~3000-4000; higher z should have fewer
MIN_EXPECTED = {
    99: 2000, 91: 2000, 84: 2000, 78: 2000,
    72: 1500, 67: 1000, 59: 800, 50: 500,
    40: 300, 33: 150, 25: 50,
}


def check():
    issues = []
    print(f"{'Snap':>5} {'z':>5} {'Ellipticals':>12} {'Gradients':>10} "
          f"{'Particles':>10} {'Status':>10}")
    print("-" * 65)

    for snap in SNAPS:
        z = cfg.SNAPSHOTS.get(snap, "?")

        ellip_file = os.path.join(
            CATDIR, f"elliptical_catalog_snap{snap:03d}.csv")
        grad_file = os.path.join(CATDIR, f"gradients_snap{snap:03d}.csv")
        snap_dir = os.path.join(DATADIR, f"snap{snap:03d}")

        n_ellip = 0
        n_grad = 0
        n_particles = 0

        if os.path.exists(ellip_file):
            df = pd.read_csv(ellip_file)
            n_ellip = len(df)

        if os.path.exists(grad_file):
            df = pd.read_csv(grad_file)
            n_grad = len(df)

        if os.path.exists(snap_dir):
            n_particles = len([
                f for f in os.listdir(snap_dir) if f.startswith("subhalo_")])

        status = "OK"
        min_exp = MIN_EXPECTED.get(snap, 50)

        if not os.path.exists(ellip_file):
            status = "MISSING"
            issues.append(f"snap {snap}: elliptical catalog missing")
        elif not os.path.exists(grad_file):
            status = "MISSING"
            issues.append(f"snap {snap}: gradient file missing")
        elif n_ellip < min_exp:
            status = "LOW"
            issues.append(
                f"snap {snap}: only {n_ellip} ellipticals "
                f"(expected >= {min_exp})")
        elif n_grad < n_ellip * 0.8:
            status = "LOW"
            issues.append(
                f"snap {snap}: only {n_grad} gradients for "
                f"{n_ellip} ellipticals ({n_grad/n_ellip:.0%})")

        print(f"{snap:>5} {z:>5} {n_ellip:>12} {n_grad:>10} "
              f"{n_particles:>10} {status:>10}")

    print()
    if issues:
        print(f"ISSUES FOUND ({len(issues)}):")
        for iss in issues:
            print(f"  - {iss}")
        print()
        print("To re-run a snapshot, delete its catalogs and run:")
        print("  rm catalogs/elliptical_catalog_snap{NN}.csv "
              "catalogs/gradients_snap{NN}.csv")
        print("  python scripts/05_redshift_evolution.py --snaps {NN}")
    else:
        print("ALL SNAPSHOTS COMPLETE")

    return len(issues) == 0


if __name__ == "__main__":
    ok = check()
    sys.exit(0 if ok else 1)
