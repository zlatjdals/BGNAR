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
R_BIN="${R_BIN:-Rscript}"

"$R_BIN" "$SCRIPT_DIR/05_collect_results.R"
"$R_BIN" "$SCRIPT_DIR/06_validate_results.R"
