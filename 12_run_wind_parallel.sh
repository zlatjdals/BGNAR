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
LOG_DIR="${WIND_LOG_DIR:-$RESULT_ROOT/wind/logs}"
MAX_JOBS="${WIND_MAX_JOBS:-15}"
TOTAL_CORES="${TOTAL_CORES:-$(nproc)}"
CORE_OFFSET="${CORE_OFFSET:-0}"
WIND_N_ORIGINS="${WIND_N_ORIGINS:-15}"
DRY_RUN="${DRY_RUN:-0}"

export BGNAR_SCRIPT_DIR="$SCRIPT_DIR"
export RESULT_ROOT WIND_N_ORIGINS
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p "$LOG_DIR" "$RESULT_ROOT/wind/checkpoints" "$RESULT_ROOT/wind/summary"

echo "[$(date '+%F %T')] Wind origins=$WIND_N_ORIGINS MAX_JOBS=$MAX_JOBS"
echo "[$(date '+%F %T')] WIND_METHODS=${WIND_METHODS:-${METHODS:-all configured methods}}; GNAR_MODE=${GNAR_MODE:-centered}"

pids=()
job_ids=()
for origin_id in $(seq 1 "$WIND_N_ORIGINS"); do
    while [[ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]]; do sleep 2; done
    core=$(( CORE_OFFSET + (origin_id - 1) % TOTAL_CORES ))
    job_id="wind_origin_$(printf '%03d' "$origin_id")"
    out="$LOG_DIR/${job_id}.out"
    err="$LOG_DIR/${job_id}.err"
    cmd=("$R_BIN" "$SCRIPT_DIR/10_run_wind_one.R" "$origin_id")
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
done

if [[ "$DRY_RUN" == "1" ]]; then
    echo "[$(date '+%F %T')] Wind dry run complete."
    exit 0
fi

status=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "[$(date '+%F %T')] FAILED: ${job_ids[$i]}" >&2
        status=1
    fi
done
if ! "$R_BIN" "$SCRIPT_DIR/11_collect_wind_results.R"; then status=1; fi
if [[ "$status" -ne 0 ]]; then
    echo "[$(date '+%F %T')] Wind analysis completed with failures." >&2
    exit "$status"
fi
echo "[$(date '+%F %T')] Wind analysis complete and validated."
