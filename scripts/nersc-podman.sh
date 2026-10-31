#!/usr/bin/env bash
# Run a command with a private Podman-HPC adapter; never modify global executables.
set -euo pipefail
command -v podman-hpc >/dev/null || { echo 'podman-hpc is unavailable on this host' >&2; exit 2; }
[[ $# -gt 0 ]] || { echo 'Usage: nersc-podman.sh COMMAND [ARGUMENTS...]' >&2; exit 2; }
# Scheduler tasks need this path on their nodes; controller-local /tmp is unsuitable.
runtime_root=${GENESIS_RUNTIME_DIR:-$PWD/workspace/runtime}
[[ "$runtime_root" == /* ]] || { echo 'GENESIS_RUNTIME_DIR must be an absolute shared path' >&2; exit 2; }
mkdir -p -- "$runtime_root"
adapter=$(mktemp -d "$runtime_root/podman.XXXXXXXX")
finish() {
    adapter_status=$?
    if (( adapter_status == 0 )); then
        rm -rf -- "$adapter"
    else
        printf 'Runtime adapter retained for recovery: %s\n' "$adapter" >&2
    fi
}
trap finish EXIT
# podman-hpc runs podman helpers by name, so the adapter calls it with the original PATH.
printf '%s\n' '#!/usr/bin/env bash' "PATH=${PATH@Q} exec podman-hpc \"\$@\"" > "$adapter/podman"
chmod 700 "$adapter/podman"
export PATH="$adapter:$PATH"
"$@"
