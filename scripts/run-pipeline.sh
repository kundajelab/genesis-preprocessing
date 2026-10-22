#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(dirname -- "$script_dir")

profile=$("$script_dir/get-default-profile.sh")
workspace=$("$script_dir/get-default-workspace.sh")
run_name=""

function usage {
    cat << EOF
Usage: $0 [ARGS] [--] INPUT.tsv [NEXTFLOW_ARGS]

Run the pipeline on exactly one sample sheet.
* Work, published output, trace files and the Nextflow log go in WORKSPACE/RUN_NAME/.
* Nextflow is launched from WORKSPACE/RUN_NAME, so -resume picks up that run's history.
  Relative paths in NEXTFLOW_ARGS resolve against that folder; prefer absolute paths.
* NEXTFLOW_ARGS are passed directly to nextflow unchecked (e.g. -resume, --slurm_queue owners).
* If launching via pixi, be sure to use full paths.

ARGS:
    -h|--help: Show this message and exit.
    -p|--profile: Use comma-separated nextflow profiles. Default inferred from environment ($profile)
    -w|--workspace: Where to put run folders. Default inferred from environment ($workspace)
    -n|--run-name: Name of the run folder in the workspace. Default: INPUT.tsv basename without .tsv
EOF
}

while (($# >= 1)); do
    case "$1" in
        -h | --help)
            usage
            exit 0
            ;;
        -p | --profile)
            profile=$2
            shift 2
            ;;
        -w | --workspace)
            workspace=$2
            shift 2
            ;;
        -n | --run-name)
            run_name=$2
            shift 2
            ;;
        --)
            shift 1
            break
            ;;
        -*)
            printf 'Unknown argument: %s\n\n' "$1" >&2
            usage >&2
            exit 2
            ;;
        *)
            break
            ;;
    esac
done

if (($# < 1)); then
    usage >&2
    exit 2
fi
sheet=$1
shift
if [[ ! -f $sheet ]]; then printf 'Missing input TSV: %s\n' "$sheet" >&2; exit 2; fi
if [[ $sheet != /* ]]; then sheet="$PWD/$sheet"; fi
# Make the workspace absolute: nextflow is launched from inside the run folder.
if [[ $workspace != /* ]]; then workspace="$PWD/$workspace"; fi
if [[ -z $run_name ]]; then
    run_name=$(basename -- "$sheet" .tsv)
fi

run_folder="$workspace/$run_name"
mkdir -p -- "$run_folder"
# Launch from the run folder so .nextflow/ history and .nextflow.log are kept with the run.
cd -- "$run_folder"
# The workflow lives in workflows/, so pass the repository's config explicitly: Nextflow only
# looks for nextflow.config beside the script and in the launch folder.
exec pixi run --manifest-path "$repo_dir" \
    nextflow run "$repo_dir/workflows/dap-seq.nf" \
    -c "$repo_dir/nextflow.config" \
    --input "$sheet" \
    -profile "$profile" \
    --workspace "$workspace" \
    --run_name "$run_name" \
    "$@"
