#!/usr/bin/env bash
# Print the default references folder for a profile (e.g. nersc,test): project storage on the
# clusters, where home folders are too small for genomes and their indexes, otherwise the
# checkout's references/.
#   nersc:    $CFS/PROJECT/$USER/references, PROJECT being your default Slurm account
#   sherlock: $OAK/$USER/genesis/references
set -euo pipefail
profile=${1:?Usage: get-default-references.sh PROFILE}
if [[ ,$profile, == *,nersc,* ]]; then
    project=$(sacctmgr -n -P show user "$USER" format=defaultaccount 2>/dev/null || true)
    # Perlmutter GPU allocations are named PROJECT_g; their storage is the project's.
    project=${project%_g}
    if [[ -z $project ]]; then
        echo "No default NERSC project found with sacctmgr; pass -R/--references" >&2
        exit 2
    fi
    echo "${CFS:-/global/cfs/cdirs}/$project/$USER/references"
elif [[ ,$profile, == *,sherlock,* ]]; then
    if [[ -z ${OAK:-} ]]; then
        echo "OAK is not set, so there is no default references folder; pass -R/--references" >&2
        exit 2
    fi
    echo "$OAK/$USER/genesis/references"
else
    root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
    echo "$root/references"
fi
