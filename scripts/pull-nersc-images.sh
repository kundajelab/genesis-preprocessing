#!/usr/bin/env bash
# Pull every image pinned in conf/genesis-images.lock with podman-hpc, on a NERSC login node.
# Compute nodes only run images that podman-hpc has pulled and squashed into $SCRATCH; tasks
# that pull for themselves run out of memory. An image already pulled costs little, and an
# interrupted earlier pull is completed.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
while IFS='=' read -r _key image; do
    podman-hpc pull "docker.io/$image"
done < <(grep '_IMAGE=' "$repo_dir/conf/genesis-images.lock")
