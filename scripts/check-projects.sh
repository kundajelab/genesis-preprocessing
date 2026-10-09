#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
shellcheck "$root"/scripts/*.sh
uv run --frozen --project "$root/genesis_tools" ruff check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests"
uv run --frozen --project "$root/genesis_tools" ruff format --check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools" "$root/tests"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_nextflow_style.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_build_tags.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_pipeline.py"
