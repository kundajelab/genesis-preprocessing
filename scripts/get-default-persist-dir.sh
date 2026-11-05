#!/usr/bin/env bash
# Print the default folder for files kept long term for a profile (e.g. nersc,test): project
# storage on the clusters, whose scratch is purged; nothing elsewhere, where runs keep their
# output in the workspace. References default to references/ beside it (get-default-references.sh).
#   nersc:    $CFS/PROJECT/$USER/runs, PROJECT being your default Slurm account
#   sherlock: $OAK/$USER/genesis/runs
# Exits 2, with a message, on a cluster where there is no default.
set -euo pipefail
profile=${1:?Usage: get-default-persist-dir.sh PROFILE}
if [[ ,$profile, == *,nersc,* ]]; then
    project=$(sacctmgr -n -P show user "$USER" format=defaultaccount 2>/dev/null || true)
    # Perlmutter GPU allocations are named PROJECT_g; their storage is the project's.
    project=${project%_g}
    if [[ -z $project ]]; then
        echo "No default NERSC project found with sacctmgr, so there is no default folder" >&2
        exit 2
    fi
    echo "${CFS:-/global/cfs/cdirs}/$project/$USER/runs"
elif [[ ,$profile, == *,sherlock,* ]]; then
    if [[ -z ${OAK:-} ]]; then
        echo "OAK is not set, so there is no default folder" >&2
        exit 2
    fi
    echo "$OAK/$USER/genesis/runs"
fi
