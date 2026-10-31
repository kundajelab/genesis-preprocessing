#!/usr/bin/env bash
# Perlmutter is detected from NERSC_HOST; any other cluster must be selected explicitly.
set -euo pipefail
if [[ ${NERSC_HOST:-} == perlmutter ]]; then
    echo "nersc"
else
    echo "local,docker"
fi
