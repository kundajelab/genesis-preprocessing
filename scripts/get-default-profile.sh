#!/usr/bin/env bash
# A cluster must always be selected explicitly.
set -euo pipefail

if command -v sbatch &> /dev/null; then
    # on slurm cluster
    if command -v podman-hpc &> /dev/null; then
        # slurm and podman, use nersc (perlmutter) profile
        profile=nersc
    elif command -v apptainer &> /dev/null; then
        # slurm and apptainer, use sherlock profile
        profile=sherlock
    else
        1>&2 echo "Unknown profile, defaulting to slurm and docker"
        profile="slurm,docker"
    fi
else
    # run local docker, check if we're on an ARM mac
    if [[ $(/usr/bin/arch) == "arm64" ]]; then
        profile="local,docker,arm"
    else
        profile="local,docker"
    fi
fi
echo "$profile"
