#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)

function usage {
    cat << EOF
Usage: $0 [--fast]

Run lint, formatting, type checks and the test scripts in tests/.
The test scripts run in parallel, at most GENESIS_CHECK_JOBS at once (default 4), and each
script runs up to that many independent Nextflow runs at once. Set GENESIS_CHECK_JOBS=1 on
shared machines.

ARGS:
    -h|--help: Show this message and exit.
    --fast: Skip the Nextflow pipeline regressions, for quick checks while editing.
EOF
}

# Worker mode, used by the parallel run below: one test script, its output in LOGS/NAME.log.
function run_test {
    local start=$SECONDS
    if uv run --frozen --project "$root/genesis_tools" python "$root/tests/$1.py" \
        > "$logs/$1.log" 2>&1; then
        printf 'PASS %s (%ss)\n' "$1" "$((SECONDS - start))"
    else
        touch -- "$logs/$1.failed"
        printf 'FAIL %s (%ss)\n' "$1" "$((SECONDS - start))"
    fi
}

fast=false
case "${1:-}" in
    --run-test)
        logs=$3
        run_test "$2"
        exit 0
        ;;
    -h | --help)
        usage
        exit 0
        ;;
    --fast) fast=true ;;
    "") ;;
    *)
        usage >&2
        exit 2
        ;;
esac
jobs=${GENESIS_CHECK_JOBS:-4}
if [[ ! $jobs =~ ^[1-9][0-9]*$ ]]; then
    printf 'GENESIS_CHECK_JOBS must be a positive integer: %s\n' "$jobs" >&2
    exit 2
fi
export GENESIS_CHECK_JOBS=$jobs

shellcheck "$root"/scripts/*.sh
uv run --frozen --project "$root/genesis_tools" ruff check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests"
uv run --frozen --project "$root/genesis_tools" ruff format --check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools" \
    --extra-search-path "$root/tests" "$root/tests"

# Slowest first, so the pipeline regressions start before the quick checks.
tests=()
if ! $fast; then
    tests+=(verify_pipeline verify_rna verify_atac verify_scatac verify_tss)
fi
tests+=(verify_deployment verify_nextflow_style verify_build_tags
    verify_multiqc verify_metro verify_readme)

# Each test writes to its own log, printed in full if it fails.
logs=$(mktemp -d)
trap 'rm -rf -- "$logs"' EXIT
printf '%s\n' "${tests[@]}" | xargs -P "$jobs" -I '{}' "$0" --run-test '{}' "$logs"

failed=0
for test in "${tests[@]}"; do
    if [[ -e $logs/$test.failed ]]; then
        printf '\n===== %s output =====\n' "$test"
        cat -- "$logs/$test.log"
        failed=1
    fi
done
exit "$failed"
