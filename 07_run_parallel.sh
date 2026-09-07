#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${SCRIPT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.env}"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi

R_BIN="${R_BIN:-Rscript}"
RESULT_ROOT="${RESULT_ROOT:-$SCRIPT_DIR/results}"
LOG_DIR="${LOG_DIR:-$RESULT_ROOT/logs}"
MAX_JOBS="${MAX_JOBS:-80}"
TOTAL_CORES="${TOTAL_CORES:-$(nproc)}"
CORE_OFFSET="${CORE_OFFSET:-0}"
REPLICATIONS="${REPLICATIONS:-100}"
DRY_RUN="${DRY_RUN:-0}"

DGP_LIST="${DGP_LIST:-M11,M12,M13,M2}"
GRAPH_LIST="${GRAPH_LIST:-ER,SBM,SWN}"
TRAIN_GRID="${TRAIN_GRID:-20,40,80,160}"

export BGNAR_SCRIPT_DIR="$SCRIPT_DIR"
export RESULT_ROOT LOG_DIR REPLICATIONS DGP_LIST GRAPH_LIST TRAIN_GRID
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "$LOG_DIR" "$RESULT_ROOT/checkpoints" "$RESULT_ROOT/summary"

IFS=',' read -r -a dgp_array <<< "$DGP_LIST"
IFS=',' read -r -a graph_array <<< "$GRAPH_LIST"
IFS=',' read -r -a train_array <<< "$TRAIN_GRID"

echo "[$(date '+%F %T')] SCRIPT_DIR=$SCRIPT_DIR"
echo "[$(date '+%F %T')] MAX_JOBS=$MAX_JOBS TOTAL_CORES=$TOTAL_CORES"
echo "[$(date '+%F %T')] DGP=$DGP_LIST GRAPH=$GRAPH_LIST T=$TRAIN_GRID REP=$REPLICATIONS"
echo "[$(date '+%F %T')] METHODS=${METHODS:-all configured methods}; GNAR_MODE=${GNAR_MODE:-centered}"

idx=0
pids=()
job_ids=()

for dgp in "${dgp_array[@]}"; do
    for graph in "${graph_array[@]}"; do
        for train_t in "${train_array[@]}"; do
            for rep in $(seq 1 "$REPLICATIONS"); do
                while [[ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]]; do
                    sleep 2
                done

                core=$(( CORE_OFFSET + idx % TOTAL_CORES ))
                job_id="${dgp}_${graph}_T$(printf '%03d' "$train_t")_r$(printf '%03d' "$rep")"
                out="$LOG_DIR/${job_id}.out"
                err="$LOG_DIR/${job_id}.err"
                cmd=("$R_BIN" "$SCRIPT_DIR/04_run_one.R" "$dgp" "$graph" "$train_t" "$rep")

                echo "[$(date '+%F %T')] Start: $job_id on core $core"
                if [[ "$DRY_RUN" == "1" ]]; then
                    printf 'taskset -c %q ' "$core"
                    printf '%q ' "${cmd[@]}"
                    printf '\n'
                elif command -v taskset >/dev/null 2>&1; then
                    nohup taskset -c "$core" "${cmd[@]}" >"$out" 2>"$err" &
                    pids+=("$!")
                    job_ids+=("$job_id")
                else
                    nohup "${cmd[@]}" >"$out" 2>"$err" &
                    pids+=("$!")
                    job_ids+=("$job_id")
                fi
                idx=$((idx + 1))
            done
        done
    done
done

if [[ "$DRY_RUN" == "1" ]]; then
    echo "[$(date '+%F %T')] Dry run complete: $idx jobs"
    exit 0
fi

status=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "[$(date '+%F %T')] FAILED: ${job_ids[$i]} (see logs)" >&2
        status=1
    fi
done

echo "[$(date '+%F %T')] Simulation jobs finished; collecting checkpoints."
if ! "$R_BIN" "$SCRIPT_DIR/05_collect_results.R"; then
    status=1
fi
if ! "$R_BIN" "$SCRIPT_DIR/06_validate_results.R"; then
    status=1
fi

if [[ "$status" -ne 0 ]]; then
    echo "[$(date '+%F %T')] Completed with failures. Check $LOG_DIR and validation_issues.csv." >&2
    exit "$status"
fi
echo "[$(date '+%F %T')] All simulations complete and validated."
