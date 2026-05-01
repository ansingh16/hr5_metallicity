#!/bin/bash
#
# Full TNG-Cluster pipeline: cluster selection, elliptical identification,
# particle download, profile computation — at ALL redshifts.
#
# Usage:
#   bash run_evolution.sh                 # all snapshots
#   bash run_evolution.sh 91 84 67        # specific snapshots only
#   bash run_evolution.sh --skip-z0       # skip z=0 steps (already done)
#
# Snapshots and redshifts:
#   99 -> z=0.0    91 -> z=0.1    84 -> z=0.2    78 -> z=0.3
#   72 -> z=0.4    67 -> z=0.5    59 -> z=0.7    50 -> z=1.0
#   40 -> z=1.5    33 -> z=2.0    25 -> z=3.0    21 -> z=4.0
#   17 -> z=5.0

set -e

cd "$(dirname "$0")"

API_WORKERS=8
DL_WORKERS=4
LOGDIR="catalogs"
mkdir -p "$LOGDIR"

SKIP_Z0=false

# Parse arguments
SNAPS=()
for arg in "$@"; do
    if [ "$arg" = "--skip-z0" ]; then
        SKIP_Z0=true
    else
        SNAPS+=("$arg")
    fi
done

# Default snapshots if none given
if [ ${#SNAPS[@]} -eq 0 ]; then
    SNAPS=(91 84 78 67 50 33 25)
fi

echo "=============================================" | tee "$LOGDIR/pipeline_run.log"
echo "TNG-Cluster Full Pipeline"                     | tee -a "$LOGDIR/pipeline_run.log"
echo "=============================================" | tee -a "$LOGDIR/pipeline_run.log"
echo "Skip z=0: $SKIP_Z0"                            | tee -a "$LOGDIR/pipeline_run.log"
echo "Evolution snapshots: ${SNAPS[*]}"               | tee -a "$LOGDIR/pipeline_run.log"
echo "API workers: $API_WORKERS"                      | tee -a "$LOGDIR/pipeline_run.log"
echo "Download workers: $DL_WORKERS"                  | tee -a "$LOGDIR/pipeline_run.log"
echo "Started: $(date)"                               | tee -a "$LOGDIR/pipeline_run.log"
echo ""                                               | tee -a "$LOGDIR/pipeline_run.log"

# ── z=0 baseline (snap 99) ─────────────────────────────────────────────

if [ "$SKIP_Z0" = false ]; then

    echo "===== Step 1: Select clusters (snap 99, z=0) =====" | tee -a "$LOGDIR/pipeline_run.log"
    if [ -f "catalogs/cluster_catalog_snap099.csv" ]; then
        echo "  Already exists, skipping." | tee -a "$LOGDIR/pipeline_run.log"
    else
        python -u scripts/01_select_clusters.py --snap 99 \
            2>&1 | tee -a "$LOGDIR/pipeline_run.log"
    fi
    echo "" | tee -a "$LOGDIR/pipeline_run.log"

    echo "===== Step 2: Select ellipticals (snap 99, z=0) =====" | tee -a "$LOGDIR/pipeline_run.log"
    if [ -f "catalogs/elliptical_catalog_snap099.csv" ]; then
        echo "  Already exists, skipping." | tee -a "$LOGDIR/pipeline_run.log"
    else
        python -u scripts/02_select_ellipticals.py --snap 99 --workers $API_WORKERS \
            2>&1 | tee -a "$LOGDIR/pipeline_run.log"
    fi
    echo "" | tee -a "$LOGDIR/pipeline_run.log"

    echo "===== Step 3: Download particles (snap 99, z=0) =====" | tee -a "$LOGDIR/pipeline_run.log"
    python -u scripts/03_extract_particles.py --snap 99 --workers $DL_WORKERS \
        2>&1 | tee -a "$LOGDIR/pipeline_run.log"
    echo "" | tee -a "$LOGDIR/pipeline_run.log"

    echo "===== Step 4: Compute profiles (snap 99, z=0) =====" | tee -a "$LOGDIR/pipeline_run.log"
    if [ -f "catalogs/gradients_snap099.csv" ]; then
        echo "  Already exists, skipping." | tee -a "$LOGDIR/pipeline_run.log"
    else
        python -u scripts/04_compute_profiles.py --snap 99 \
            2>&1 | tee -a "$LOGDIR/pipeline_run.log"
    fi
    echo "" | tee -a "$LOGDIR/pipeline_run.log"

fi

# ── Redshift evolution (other snapshots) ────────────────────────────────

echo "===== Step 5: Redshift evolution =====" | tee -a "$LOGDIR/pipeline_run.log"
echo "  Snapshots: ${SNAPS[*]}" | tee -a "$LOGDIR/pipeline_run.log"
echo "" | tee -a "$LOGDIR/pipeline_run.log"

python -u scripts/05_redshift_evolution.py \
    --snaps ${SNAPS[*]} \
    --api-workers $API_WORKERS \
    --dl-workers $DL_WORKERS \
    2>&1 | tee -a "$LOGDIR/pipeline_run.log"

# ── Summary ─────────────────────────────────────────────────────────────

echo "" | tee -a "$LOGDIR/pipeline_run.log"
echo "=============================================" | tee -a "$LOGDIR/pipeline_run.log"
echo "Pipeline complete: $(date)"                     | tee -a "$LOGDIR/pipeline_run.log"
echo "=============================================" | tee -a "$LOGDIR/pipeline_run.log"
echo "" | tee -a "$LOGDIR/pipeline_run.log"
echo "Output files:" | tee -a "$LOGDIR/pipeline_run.log"
ls -lh catalogs/gradients_*.csv catalogs/profiles_*.hdf5 2>/dev/null | tee -a "$LOGDIR/pipeline_run.log"
echo "" | tee -a "$LOGDIR/pipeline_run.log"
echo "Combined evolution table:" | tee -a "$LOGDIR/pipeline_run.log"
wc -l catalogs/gradients_evolution.csv 2>/dev/null | tee -a "$LOGDIR/pipeline_run.log"
echo "Log: $LOGDIR/pipeline_run.log" | tee -a "$LOGDIR/pipeline_run.log"
