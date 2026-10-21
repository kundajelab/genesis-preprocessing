#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
shellcheck "$root"/scripts/*.sh "$root"/benchmarks/*.sh
uv run --frozen --project "$root/genesis_tools" ruff check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests" "$root/examples/curation"
uv run --frozen --project "$root/genesis_tools" ruff format --check \
    --config "$root/genesis_tools/pyproject.toml" "$root/genesis_tools" "$root/tests" "$root/examples/curation"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools"
uv run --frozen --project "$root/genesis_tools" ty check --project "$root/genesis_tools" \
    --extra-search-path "$root/tests" "$root/tests"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_nextflow_style.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_publication.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_multiqc.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_build_tags.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_pipeline.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_benchmark.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_benchmark_spec.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_benchmark_methods.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_benchmark_references.py"
uv run --frozen --project "$root/genesis_tools" python "$root/benchmarks/experiments/reference-comparison-20261014/check_reference_comparison.py"
uv run --frozen --project "$root/genesis_tools" python "$root/benchmarks/experiments/reference-comparison-20261014/check_sra_sampling.py"
uv run --frozen --project "$root/genesis_tools" python "$root/benchmarks/experiments/reference-comparison-20261014/check_queue_stop.py"
uv run --frozen --project "$root/genesis_tools" python "$root/benchmarks/experiments/dap-finalization/check_pilot.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_atac.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_metro.py"

uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_execution.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_cli.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_supported_atac.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_deployment.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_curation.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_curation_artifacts.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_curation_adapters.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_ai.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_llm_deployment.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/evaluate_metadata.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_readme.py"

uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_study.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_scatac_inventory.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_scatac_fragments.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_scatac_pseudobulk.py"
uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_scatac_exports.py"

uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_scatac_inspection.py"

uv run --frozen --project "$root/genesis_tools" python "$root/tests/verify_macs3_launcher.py"
