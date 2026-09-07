#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${SCRIPT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config.env}"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi
export BGNAR_SCRIPT_DIR="$SCRIPT_DIR"
export RESULT_ROOT="${RESULT_ROOT:-$SCRIPT_DIR/results}"
"${R_BIN:-Rscript}" "$SCRIPT_DIR/11_collect_wind_results.R"
