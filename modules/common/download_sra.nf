include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from './retry'
include { batchTag ; countsPath ; readUrls ; sraRun ; sraTime } from './download_support'

// Each sample's SRA run's .sra file from AWS, checked against NCBI's MD5 and converted to FASTQs
// named and formatted as ENA's (header "@RUN.N NAME/MATE", bare "+"): the same bytes, decompressed.
// A sample whose run is unavailable (an HTTP failure, a stalled transfer, an MD5 mismatch) or
// whose converted reads are not the FASTQs its URLs name (another layout, broken records) is left
// out, its partial files deleted, and recorded in status.tsv as "fallback" with the reason, for
// DOWNLOAD_URLS to fetch instead; the task still succeeds. Any other failure (fasterq-dump, pigz,
// a full disk) stops the task. status.tsv has one line per sample: ID, "ok" or "fallback", reason.
process DOWNLOAD_SRA {
    // AWS's copy takes many connections; NCBI's locator is asked once per run.
    maxForks 32
    // fasterq-dump and pigz use them all.
    cpus 8
    memory { oomMemoryOf(task, 2.GB) }
    time { retryTimeOf(task, sraTime(samples)) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${batchTag(samples)}"
    publishDir "${params.publish_base}/output", mode: 'copy', pattern: '*.fastq.lines',
        saveAs: { name -> countsPath(samples, name) }
    input:
    val samples
    output:
    tuple val(samples), path('status.tsv'), path('*.read*.fastq.gz', arity: '0..*'), emit: downloads
    path '*.fastq.lines', arity: '0..*', emit: counts
    script:
    // One line per sample: ID, SRA run, number of mates.
    def plan = samples.collect { meta -> [meta.id, sraRun(meta), readUrls(meta).size()].join('\t') }
        .join('\n')
    """
    set -o pipefail
    # sra-tools keeps its settings outside the home folder.
    export NCBI_SETTINGS="\$PWD/ncbi-settings.mkfg"
    # Why the last run was unavailable, for status.tsv: set before fetch_sra returns 2.
    reason=

    # Fetch and convert one run. Returns 2 when the run is unavailable, for the sample to be
    # fetched from its URLs instead, and 1 on any other failure.
    fetch_sra() {
        local id=\$1 run=\$2 mates=\$3 located link md5 mate expected file try status
        located=\$(curl --fail --silent --show-error --location --connect-timeout 60 \\
            --retry 3 --retry-all-errors --retry-delay 10 '${params.sra_locate_url}'"\$run") || {
            reason="NCBI's locator failed"
            return 2
        }
        # The first file of type sra: its MD5, then its first location's link.
        located=\$(grep -o '"type": *"sra".*' <<< "\$located") || {
            reason="NCBI's locator lists no .sra file"
            return 2
        }
        md5=\$(grep -o '"md5": *"[0-9a-f]\\{32\\}"' <<< "\$located" | head -n 1 | grep -o '[0-9a-f]\\{32\\}')
        link=\$(grep -o '"link": *"[^"]*"' <<< "\$located" | head -n 1 | sed 's/^"link": *"//; s/"\$//')
        [[ -n \$md5 && -n \$link ]] || {
            reason="NCBI's locator gave no MD5 or link"
            return 2
        }
        # Resume after a failed try, up to 5 tries. A transfer slower than 2 MB/s for 5 minutes
        # gives up at once: a run that keeps up with that fits the second retry's time.
        for try in 1 2 3 4 5; do
            status=0
            curl --fail --silent --show-error --location --connect-timeout 60 \\
                --speed-limit 2000000 --speed-time 300 \\
                --continue-at - --output "\$run.sra" "\$link" || status=\$?
            case \$status in
                0) break ;;
                # Not the server's doing: the file could not be written, or curl ran out of memory.
                23 | 27)
                    echo "\$run.sra: curl failed locally (exit \$status)" >&2
                    return 1
                    ;;
                28)
                    reason="the .sra transfer timed out or stayed under 2 MB/s for 5 minutes"
                    return 2
                    ;;
            esac
            reason="the .sra transfer failed (curl exit \$status)"
            (( try < 5 )) || return 2
            sleep 30
        done
        [[ \$(md5sum < "\$run.sra" | cut -d ' ' -f 1) == "\$md5" ]] || {
            echo "\$run.sra: MD5 differs from NCBI's \$md5" >&2
            reason="the .sra file's MD5 differs from NCBI's"
            return 2
        }
        fasterq-dump --split-3 --threads ${task.cpus} --temp "tmp-\$run" --outdir "sra-\$run" \\
            --seq-defline '@\$ac.\$si \$sn/\$ri' --qual-defline '+' "\$run.sra" >&2 || {
            echo "\$run: fasterq-dump failed" >&2
            return 1
        }
        rm "\$run.sra"
        # Exactly the files ENA serves for this layout: RUN.fastq, or RUN_1 and RUN_2.
        if (( mates == 1 )); then expected="\$run.fastq"; else expected="\${run}_1.fastq \${run}_2.fastq"; fi
        [[ \$(cd "sra-\$run" && ls | sort | tr '\\n' ' ') == "\$expected " ]] || {
            echo "\$run: fasterq-dump wrote \$(ls "sra-\$run"), not \$expected" >&2
            reason="its .sra file is not of the sample's layout"
            return 2
        }
        for mate in \$(seq 1 "\$mates"); do
            if (( mates == 1 )); then file="sra-\$run/\$run.fastq"; else file="sra-\$run/\${run}_\$mate.fastq"; fi
            wc -l < "\$file" | tr -d ' ' > "\$id.read\$mate.lines" || return 1
            # Level 4 compresses 5.5 times as fast as the default (6), into files 10% larger.
            pigz -p ${task.cpus} -4 -c "\$file" > "\$id.read\$mate.fastq.gz" || return 1
            rm "\$file"
        done
        rm -r "sra-\$run" "tmp-\$run"
        # Whole FASTQ records, at least one, and as many in each mate.
        awk '{if (\$1 == 0 || \$1 % 4) exit 1; print \$1}' "\$id.read1.lines" > "\$id.fastq.lines" &&
            { (( mates == 1 )) || cmp -s "\$id.read1.lines" "\$id.read2.lines"; } || {
            reason="its converted reads are not whole, matching FASTQ records"
            return 2
        }
        rm "\$id".read*.lines
    }

    cat > plan.tsv << 'END_OF_PLAN'
${plan}
END_OF_PLAN
    : > status.tsv
    # On its own descriptor, so that no command in the loop can read the plan as its input.
    while IFS=\$'\\t' read -r -u 4 id run mates; do
        status=0
        fetch_sra "\$id" "\$run" "\$mates" || status=\$?
        if (( status == 2 )); then
            echo "WARNING: \$id: SRA run \$run unavailable (\$reason); its FASTQs will be downloaded instead" >&2
            rm -rf "\$run.sra" "sra-\$run" "tmp-\$run" "\$id".read* "\$id.fastq.lines"
            printf '%s\\tfallback\\t%s\\n' "\$id" "\$reason" >> status.tsv
        elif (( status == 0 )); then
            printf '%s\\tok\\t\\n' "\$id" >> status.tsv
        else
            exit "\$status"
        fi
    done 4< plan.tsv
    """
    stub:
    def fakes = samples.collect { meta ->
        def second = meta.layout == 'PE' ? "\ncp '${meta.id}.read1.fastq.gz' '${meta.id}.read2.fastq.gz'" : ''
        "printf '@read1\\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip > '${meta.id}.read1.fastq.gz'\n" +
            "printf '4\\n' > '${meta.id}.fastq.lines'\n" +
            "printf '%s\\tok\\t\\n' '${meta.id}' >> status.tsv" + second
    }.join('\n')
    """
    : > status.tsv
${fakes}
    """
}
