#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
case "${1:-}" in
    "") ;;
    -h|--help)
        printf '%s\n' 'Usage: pixi run install-all' 'Install locked dependencies.'
        exit 0
        ;;
    *) printf 'Unknown installation option: %s\n' "$1" >&2; exit 2 ;;
esac
pixi install --manifest-path "$root/pixi.toml" --locked
# genesis_tools pins a "+gil" Python, which uv selects but cannot download; install the plain
# version so a managed free-threaded build of the same version is never chosen instead.
uv python install "$(sed 's/+gil$//' "$root/genesis_tools/.python-version")"
uv sync --project "$root/genesis_tools" --locked
# YAML environments are created by Nextflow's conda profile, or by Docker builds.
