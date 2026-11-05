#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

# The Nextflow manager mostly waits on Slurm, so it needs a long walltime rather than many CPUs.
manager_queue=""
time_limit="2-00:00:00"
cpus=2
memory="8G"
log_dir="$HOME/logs/genesis"
job_name=""

function usage {
    cat << EOF
Usage: $0 [ARGS] [--] [RUN_PIPELINE_ARGS] INPUT.tsv [NEXTFLOW_ARGS]

Submit scripts/run-pipeline.sh as a Slurm job, so the Nextflow manager runs in the queue instead
of on the login node; it then submits each task as its own job. Returns once the job is queued.
* ARGS come first. Everything else is passed to run-pipeline.sh unchanged (see its --help), e.g.
  -p, -w, -R, -P, -n, --workflow, the sample sheet, then Nextflow arguments such as -resume.
* The job runs from the current folder, so relative paths work as with run-pipeline.sh.
* The walltime caps the whole run. A resubmitted run resumes where the last one stopped.

ARGS:
    -h|--help: Show this message and exit.
    -Q|--manager-queue: Where the manager job runs (not the tasks; for those pass --queue to
        Nextflow): the QOS with nersc (default shared, at most 48 hours; workflow, for runs of up
        to 90 days, needs NERSC's approval and may require scrontab instead of sbatch), the
        partition with sherlock (default normal) or slurm (default regular).
    -T|--time: Walltime of the manager job. Default: $time_limit
    -C|--cpus: CPUs for the manager job. Default: $cpus
    -M|--mem: Memory for the manager job; Nextflow gets 3/4 of it. Default: $memory
    -L|--log-dir: Folder for the manager's log. Default: $log_dir
    -J|--job-name: Slurm job name. Default: genesis-INPUT, INPUT the sample sheet's basename
EOF
}

while (($# >= 1)); do
    case "$1" in
        -h | --help)
            usage
            exit 0
            ;;
        -Q | --manager-queue)
            manager_queue=$2
            shift 2
            ;;
        -T | --time)
            time_limit=$2
            shift 2
            ;;
        -C | --cpus)
            cpus=$2
            shift 2
            ;;
        -M | --mem)
            memory=$2
            shift 2
            ;;
        -L | --log-dir)
            log_dir=$2
            shift 2
            ;;
        -J | --job-name)
            job_name=$2
            shift 2
            ;;
        --)
            shift 1
            break
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

# Read the profile and sample sheet from run-pipeline.sh's arguments, without changing them.
profile=""
sheet=""
arguments=("$@")
i=0
while ((i < ${#arguments[@]})); do
    case "${arguments[i]}" in
        -p | --profile)
            profile=${arguments[i + 1]:-}
            i=$((i + 2))
            ;;
        -w | --workspace | -R | --references | -P | --persist-dir | -n | --run-name | --workflow)
            i=$((i + 2))
            ;;
        --)
            sheet=${arguments[i + 1]:-}
            break
            ;;
        -*)
            i=$((i + 1))
            ;;
        *)
            sheet=${arguments[i]}
            break
            ;;
    esac
done
if [[ -z $profile ]]; then
    profile=$("$script_dir/get-default-profile.sh")
fi
if [[ -z $job_name ]]; then
    job_name="genesis-$(basename -- "${sheet:-run}" .tsv)"
fi

# Where the manager job goes, by site.
site_options=()
environment=""
case ",$profile," in
    *,nersc,*)
        manager_queue=${manager_queue:-shared}
        if [[ $manager_queue == workflow ]]; then
            site_options=(--qos=workflow --constraint=cron)
        else
            site_options=("--qos=$manager_queue" --constraint=cpu)
        fi
        # Pull here, on the login node: pulls on compute nodes run out of memory.
        "$script_dir/pull-nersc-images.sh"
        environment="export GENESIS_IMAGES_PULLED=1"
        ;;
    *,sherlock,*)
        site_options=("--partition=${manager_queue:-normal}")
        ;;
    *,slurm,*)
        site_options=("--partition=${manager_queue:-regular}")
        ;;
    *)
        printf 'submit-pipeline.sh needs a Slurm profile (nersc, sherlock or slurm), not %s;\n' \
            "$profile" >&2
        printf 'pass one with -p, or run scripts/run-pipeline.sh directly.\n' >&2
        exit 2
        ;;
esac
if ! command -v sbatch > /dev/null; then
    printf 'sbatch not found: submit-pipeline.sh runs on a Slurm login node.\n' >&2
    exit 2
fi

# Nextflow's JVM gets 3/4 of the job's memory, leaving room for the shell and pixi.
heap=""
if [[ $memory =~ ^([0-9]+)[Gg]$ ]]; then
    heap="-Xmx$((BASH_REMATCH[1] * 3 / 4 > 0 ? BASH_REMATCH[1] * 3 / 4 : 1))g"
elif [[ $memory =~ ^([0-9]+)[Mm]?$ ]]; then
    heap="-Xmx$((BASH_REMATCH[1] * 3 / 4))m"
fi

command_line=$(printf '%q ' "$script_dir/run-pipeline.sh" "$@")
script=$(
    cat << EOF
#!/usr/bin/env bash
set -euo pipefail
echo "Nextflow manager: job \${SLURM_JOB_ID:-?} on \$(hostname)"
# Task jobs that Nextflow submits from here would inherit this job's SLURM_* settings (CPUs,
# memory, partition). Keep only SLURM_CONF, which the Slurm commands need.
unset \$(compgen -v SLURM_ | grep -vx SLURM_CONF) || true
# The log is a file, so print progress lines rather than redrawing the terminal.
export NXF_ANSI_LOG=false
export NXF_OPTS="$heap\${NXF_OPTS:+ \$NXF_OPTS}"
$environment
# Slurm sends this shell SIGTERM 2 minutes before the walltime (--signal); pass it on to
# Nextflow, which then cancels its task jobs rather than leaving them in the queue. pixi forwards
# it to Nextflow and waits for it to finish.
$command_line&
pid=\$!
trap 'kill -TERM "\$pid" 2> /dev/null || true' TERM
status=0
while :; do
    wait "\$pid" && status=0 || status=\$?
    kill -0 "\$pid" 2> /dev/null || break
done
exit "\$status"
EOF
)

mkdir -p -- "$log_dir"
# A requeued job would start the run again from scratch; resubmit instead, which resumes it.
job_id=$(printf '%s\n' "$script" | sbatch --parsable \
    --job-name="$job_name" \
    --output="$log_dir/%x-%j.out" \
    --time="$time_limit" \
    --cpus-per-task="$cpus" \
    --mem="$memory" \
    --chdir="$PWD" \
    --no-requeue \
    --signal=B:TERM@120 \
    "${site_options[@]}")

cat << EOF
Submitted the Nextflow manager as job $job_id (${site_options[*]}, $time_limit, $cpus CPUs, $memory).
    log:    $log_dir/$job_name-$job_id.out
    status: squeue -j $job_id
EOF
