#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(dirname -- "$script_dir")

profile=$("$script_dir/get-default-profile.sh")
workspace=$("$script_dir/get-default-workspace.sh")
references=""
run_name=""
workflow="dap-seq"

function usage {
    cat << EOF
Usage: $0 [ARGS] [--] INPUT.tsv [NEXTFLOW_ARGS]

Run a pipeline on exactly one sample sheet.
* Work, published output, trace files and the Nextflow log go in WORKSPACE/RUN_NAME/.
* Nextflow is launched from WORKSPACE/RUN_NAME, so -resume picks up that run's history.
  Relative paths in NEXTFLOW_ARGS resolve against that folder; prefer absolute paths.
* NEXTFLOW_ARGS are passed directly to nextflow unchecked (e.g. -resume, --queue owners).
* If launching via pixi, be sure to use full paths.

ARGS:
    -h|--help: Show this message and exit.
    -p|--profile: Use comma-separated nextflow profiles. Default inferred from environment ($profile)
    -w|--workspace: Where to put run folders. Default inferred from environment ($workspace)
    -R|--references: Folder holding the FASTAs and annotations, where runs also publish reference
        indexes; created if missing. Default from scripts/get-default-references.sh, by profile:
        \$CFS/PROJECT/\$USER/references with nersc (PROJECT: your default Slurm account),
        \$OAK/\$USER/genesis/references with sherlock, otherwise the checkout's references/
    -n|--run-name: Name of the run folder in the workspace. Default: INPUT.tsv basename without .tsv
    --workflow: Workflow in workflows/ to run: dap-seq, atac-seq, scatac-seq, scatac-reads, rna-seq or tss-seq. Default: $workflow
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
        -R | --references)
            references=$2
            shift 2
            ;;
        -n | --run-name)
            run_name=$2
            shift 2
            ;;
        --workflow)
            workflow=$2
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
if [[ ! -f $repo_dir/workflows/$workflow.nf || ! -f $repo_dir/conf/$workflow.config ]]; then
    printf 'Unknown workflow: %s\n' "$workflow" >&2
    exit 2
fi
if [[ $sheet != /* ]]; then sheet="$PWD/$sheet"; fi
# Make the workspace absolute: nextflow is launched from inside the run folder.
if [[ $workspace != /* ]]; then workspace="$PWD/$workspace"; fi
# The default is looked up only when needed, so -R works where there is none (e.g. no $OAK).
if [[ -z $references ]]; then
    references=$("$script_dir/get-default-references.sh" "$profile")
fi
if [[ $references != /* ]]; then references="$PWD/$references"; fi
# Runs publish reference products such as indexes into this folder, so it must exist.
mkdir -p -- "$references"
if [[ -z $run_name ]]; then
    run_name=$(basename -- "$sheet" .tsv)
fi

if [[ ,$profile, == *,nersc,* ]]; then
    # Compute nodes only run images that podman-hpc has pulled and squashed into $SCRATCH; tasks
    # that pull for themselves run out of memory. Pull every pinned image here, on the login
    # node: an image already pulled costs little, and an interrupted earlier pull is completed.
    while IFS='=' read -r _key image; do
        podman-hpc pull "docker.io/$image"
    done < <(grep '_IMAGE=' "$repo_dir/conf/genesis-images.lock")
    # The `podman` Nextflow runs must be podman-hpc (see scripts/podman-hpc-bin/podman).
    export PATH="$script_dir/podman-hpc-bin:$PATH"
fi

run_folder="$workspace/$run_name"
mkdir -p -- "$run_folder"
# Launch from the run folder so .nextflow/ history and .nextflow.log are kept with the run.
cd -- "$run_folder"
# The workflow lives in workflows/, so pass its config explicitly: Nextflow only looks for
# nextflow.config beside the script and in the launch folder. conf/WORKFLOW.config includes the
# shared nextflow.config.
exec pixi run --manifest-path "$repo_dir" \
    nextflow run "$repo_dir/workflows/$workflow.nf" \
    -c "$repo_dir/conf/$workflow.config" \
    --input "$sheet" \
    -profile "$profile" \
    --workspace "$workspace" \
    --run_name "$run_name" \
    --references "$references" \
    "$@"
