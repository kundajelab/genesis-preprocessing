#!/usr/bin/env bash
# Print the default references folder for a profile (e.g. nersc,test): references/ beside the
# default persist folder (scripts/get-default-persist-dir.sh, ending in runs/) on the clusters,
# where home folders are too small for genomes and their indexes, else the checkout's references/.
#   nersc:    $CFS/PROJECT/$USER/references, PROJECT being your default Slurm account
#   sherlock: $OAK/$USER/genesis/references
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
profile=${1:?Usage: get-default-references.sh PROFILE}
if ! persist_dir=$("$script_dir/get-default-persist-dir.sh" "$profile"); then
    echo "No default references folder; pass -R/--references" >&2
    exit 2
fi
if [[ -n $persist_dir ]]; then
    echo "$(dirname -- "$persist_dir")/references"
else
    echo "$(dirname -- "$script_dir")/references"
fi
