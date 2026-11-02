#!/usr/bin/env bash
# Print this host's default profile: nersc on Perlmutter (NERSC_HOST, or Slurm with podman-hpc),
# sherlock for Slurm with Apptainer, otherwise local Docker.
set -euo pipefail

if [[ ${NERSC_HOST:-} == perlmutter ]]; then
    profile=nersc
elif command -v sbatch &> /dev/null; then
    # on slurm cluster
    if command -v podman-hpc &> /dev/null; then
        # slurm and podman, use nersc (perlmutter) profile
        profile=nersc
    elif command -v apptainer &> /dev/null; then
        # slurm and apptainer, use sherlock profile
        profile=sherlock
    else
        1>&2 echo "Slurm found without podman-hpc or Apptainer; defaulting to local,docker." \
            "Choose a profile with -p to use the cluster."
        profile="local,docker"
    fi
else
    profile="local,docker"
fi
echo "$profile"
