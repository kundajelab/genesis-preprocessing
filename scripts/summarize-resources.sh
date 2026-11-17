#!/usr/bin/env bash
# Summarize a run's resource use in PERSIST_RUN/resources.tsv, from its Nextflow trace files, its
# tasks' work folders, Slurm accounting and du. Light enough for a login node.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

function usage {
    cat << EOF
Usage: $0 PERSIST_RUN [WORK_DIR]

Write PERSIST_RUN/resources.tsv: one line per field, its name and value (%.3f) tab-separated:
    cpu-hours           Σ requested CPUs × realtime
    mem-gb-hours        Σ requested memory × realtime
    node-hours          Σ Slurm billing × elapsed / 256 (Perlmutter CPU node threads), over the
                        tasks and the Nextflow manager jobs; nan without Slurm accounting
    scratch-gb          du of WORK_DIR
    persistent-disk-gb  du of PERSIST_RUN
    gb-written          Σ wchar (every write, including to pipes and logs)
    avg-%cpu            Σ(%cpu ÷ cpus × cpu-hours) ÷ Σ cpu-hours
    avg-%mem            Σ(peak_rss ÷ memory × 100 × mem-gb-hours) ÷ Σ mem-gb-hours
GB is 2^30 bytes, as in Nextflow. Every task attempt with a measured realtime counts, from every
PERSIST_RUN/trace/execution_trace_*.txt; cached rows repeat earlier attempts and are skipped.
Requests come from the trace's cpus and memory columns, with any executor (Slurm, or local with
conda or a container). Traces without them, from older runs, need the tasks' work folders: the
requests are then Slurm's directives or the container's limits in each .command.run. Only
scripts/submit-pipeline.sh runs this when a run finishes; run it yourself after other runs.

PERSIST_RUN: the run's folder of published output and trace files (PERSIST_DIR/RUN_NAME).
WORK_DIR: the run's Nextflow work folder. Default: WORKSPACE/RUN_NAME/work, WORKSPACE from
    scripts/get-default-workspace.sh and RUN_NAME the basename of PERSIST_RUN.
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    usage
    exit 0
fi
if (($# < 1 || $# > 2)); then
    usage >&2
    exit 2
fi
persist=$(cd -- "$1" && pwd)
run_name=$(basename -- "$persist")
work=${2:-$("$script_dir/get-default-workspace.sh")/$run_name/work}
if [[ ! -d $work ]]; then
    printf 'No work folder %s; pass WORK_DIR\n' "$work" >&2
    exit 1
fi
work=$(cd -- "$work" && pwd)
traces=("$persist"/trace/execution_trace_*.txt)
if [[ ! -e ${traces[0]} ]]; then
    printf 'No trace files in %s/trace\n' "$persist" >&2
    exit 1
fi

temporary=$(mktemp -d)
trap 'rm -rf -- "$temporary"' EXIT

# Task attempts that ran, once each (a hash is one attempt): hash, native_id, realtime (hours),
# %cpu, peak_rss (bytes), wchar (bytes), and the requested CPUs and memory (MB) where the trace
# records them, with "-" where the trace has none.
awk -F '\t' -v OFS='\t' '
    function hours(text,    total, n, i, part, value) {
        if (text == "-" || text == "") return "-"
        n = split(text, part, " ")
        for (i = 1; i <= n; i++) {
            value = part[i] + 0
            if (part[i] ~ /ms$/) total += value / 3600000
            else if (part[i] ~ /s$/) total += value / 3600
            else if (part[i] ~ /m$/) total += value / 60
            else if (part[i] ~ /h$/) total += value
            else if (part[i] ~ /d$/) total += value * 24
        }
        return total
    }
    function bytes(text,    n, part, scale) {
        if (text == "-" || text == "") return "-"
        n = split(text, part, " ")
        scale = n < 2 ? 1 : part[2] == "KB" ? 2^10 : part[2] == "MB" ? 2^20 : \
            part[2] == "GB" ? 2^30 : part[2] == "TB" ? 2^40 : part[2] == "PB" ? 2^50 : 1
        return part[1] * scale
    }
    FNR == 1 {
        delete column
        for (i = 1; i <= NF; i++) column[$i] = i
        for (name in required) if (!(name in column)) {
            printf "%s has no %s column\n", FILENAME, name > "/dev/stderr"
            failed = 1
            exit 1
        }
        next
    }
    BEGIN { split("hash native_id status realtime %cpu peak_rss wchar", names, " ")
        for (i in names) required[names[i]] = 1 }
    $column["status"] != "CACHED" && $column["realtime"] != "-" {
        cpu = $column["%cpu"]; sub(/%$/, "", cpu)
        cpus = ("cpus" in column ? $column["cpus"] : "-")
        memory = ("memory" in column ? bytes($column["memory"]) : "-")
        if (cpus !~ /^[0-9]+$/ || memory == "-") cpus = memory = "-"
        else memory = memory / 2^20
        row[$column["hash"]] = $column["native_id"] OFS hours($column["realtime"]) OFS cpu OFS \
            bytes($column["peak_rss"]) OFS bytes($column["wchar"]) OFS cpus OFS memory
    }
    END { if (failed) exit 1
        for (hash in row) print hash, row[hash] }
' "${traces[@]}" > "$temporary/tasks.tsv"
if [[ ! -s $temporary/tasks.tsv ]]; then
    printf 'No task in %s ran\n' "$persist/trace" >&2
    exit 1
fi

# Each attempt's work folder (WORK/ab/cdef..., the trace's hash ab/cdef12 abbreviated), then the
# CPUs and memory (MB) it requested: Slurm's directives, or else the container's limits. Only an
# attempt whose trace row has no requests needs its folder; the others' show whether Slurm ran it.
find "$work" -mindepth 2 -maxdepth 2 -type d > "$temporary/folders.txt"
awk -F '\t' -v OFS='\t' -v work="$work/" '
    NR == FNR { key = substr($1, 1, 9); folder[key] = ""; needed[key] = ($7 == "-"); next }
    { relative = substr($0, length(work) + 1); key = substr(relative, 1, 9)
      if (key in folder) folder[key] = $0 }
    END { for (key in folder) {
        if (folder[key] != "") print key, folder[key]
        else if (needed[key]) { printf "No work folder for %s in %s\n", key, work > "/dev/stderr"
            failed = 1 } }
      exit failed }
' "$temporary/tasks.tsv" "$temporary/folders.txt" > "$temporary/paths.tsv"
# The $ fields are awk's, after xargs, not the shell's.
# shellcheck disable=SC2016
cut -f 2 "$temporary/paths.tsv" | sed 's|$|/.command.run|' | tr '\n' '\0' \
    | xargs -0 awk -v OFS='\t' '
        function emit() {
            if (prev == "") return
            cpus = (slurm_cpus != "" ? slurm_cpus : (shares != "" ? shares / 1024 : ""))
            memory = (slurm_memory != "" ? slurm_memory : container_memory)
            print prev, cpus, memory, (slurm_cpus != "" ? "slurm" : "-")
        }
        FNR == 1 { emit(); prev = FILENAME; slurm_cpus = slurm_memory = shares = container_memory = "" }
        /^#SBATCH -c [0-9]+$/ { slurm_cpus = $3 }
        /^#SBATCH --mem [0-9]+M$/ { slurm_memory = $3 + 0 }
        match($0, /--cpu-shares [0-9]+/) { shares = substr($0, RSTART + 13, RLENGTH - 13) }
        match($0, /--memory [0-9]+m/) { container_memory = substr($0, RSTART + 9, RLENGTH - 10) }
        END { emit() }
    ' > "$temporary/requests.tsv"

# Slurm's billing for the tasks and for the Nextflow manager jobs whose logs name this run.
node_hours=nan
if command -v sacct > /dev/null && grep -q $'\tslurm$' "$temporary/requests.tsv"; then
    cut -f 2 "$temporary/tasks.tsv" > "$temporary/jobs.txt"
    for log in "$HOME"/logs/genesis/*.out; do
        [[ -e $log ]] || continue
        if awk -v persist="$persist" '$1 == "publish_base" && $3 == persist {found = 1; exit}
            END {exit !found}' "$log"; then
            awk '/^Nextflow manager: job [0-9]+/ {print $4; exit}' "$log"
        fi
    done >> "$temporary/jobs.txt"
    sort -u "$temporary/jobs.txt" | { grep -E '^[0-9]+$' || true; } | xargs -n 500 \
        | while read -r jobs; do
            sacct -X -n -P -o ElapsedRaw,AllocTRES -j "${jobs// /,}"
        done > "$temporary/billing.txt"
    node_hours=$(awk -F '|' '
        match($2, /billing=[0-9]+/) { total += substr($2, RSTART + 8, RLENGTH - 8) * $1 }
        END { printf "%.3f", total / 3600 / 256 }' "$temporary/billing.txt")
fi

scratch_kb=$(du -sk -- "$work" | cut -f 1)
persistent_kb=$(du -sk -- "$persist" | cut -f 1)

awk -F '\t' -v OFS='\t' -v node_hours="$node_hours" -v scratch_kb="$scratch_kb" \
    -v persistent_kb="$persistent_kb" -v work="$work/" '
    FILENAME ~ /paths.tsv$/ { path[$2 "/.command.run"] = $1; next }
    FILENAME ~ /requests.tsv$/ {
        if ($2 != "" && $3 != "") { cpus[path[$1]] = $2; memory[path[$1]] = $3 / 1024 }
        else folder[path[$1]] = $1
        next
    }
    {
        key = substr($1, 1, 9); realtime = $3
        # Requests from the trace, or else from the work folder.
        if ($7 != "-") { cpus[key] = $7; memory[key] = $8 / 1024 }
        else if (!(key in cpus)) { printf "No CPU or memory request in %s\n", folder[key] > "/dev/stderr"
            failed = 1; exit 1 }
        cpu_hours = cpus[key] * realtime; memory_hours = memory[key] * realtime
        total_cpu += cpu_hours; total_memory += memory_hours
        if ($6 != "-") written += $6
        if ($4 != "-" && $4 != "") { cpu_sum += $4 / cpus[key] * cpu_hours; cpu_weight += cpu_hours }
        if ($5 != "-") { memory_sum += $5 / 2^30 / memory[key] * 100 * memory_hours
            memory_weight += memory_hours }
    }
    END {
        if (failed) exit 1
        printf "cpu-hours\t%.3f\n", total_cpu
        printf "mem-gb-hours\t%.3f\n", total_memory
        printf "node-hours\t%s\n", node_hours
        printf "scratch-gb\t%.3f\n", scratch_kb / 2^20
        printf "persistent-disk-gb\t%.3f\n", persistent_kb / 2^20
        printf "gb-written\t%.3f\n", written / 2^30
        printf "avg-%%cpu\t%s\n", (cpu_weight ? sprintf("%.3f", cpu_sum / cpu_weight) : "nan")
        printf "avg-%%mem\t%s\n", (memory_weight ? sprintf("%.3f", memory_sum / memory_weight) : "nan")
    }
' "$temporary/paths.tsv" "$temporary/requests.tsv" "$temporary/tasks.tsv" > "$temporary/resources.tsv"
mv -- "$temporary/resources.tsv" "$persist/resources.tsv"
cat -- "$persist/resources.tsv"
