#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
shellcheck "$root"/scripts/*.sh
uv run --frozen --project "$root/genesis_tools" ruff check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests"
uv run --frozen --project "$root/genesis_tools" ruff format --check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools" \
    --extra-search-path "$root/tests" "$root/tests"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_nextflow_style.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_publication.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_multiqc.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_build_tags.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_pipeline.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_atac.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_metro.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_supported_atac.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_deployment.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_readme.py"
