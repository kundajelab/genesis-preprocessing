include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'
include { sraReadsTime } from '../common/download_support'

// One SRA run of a library: its .sra file from AWS, checked against NCBI's MD5 (sra-fetch), then
// its two genomic mates and its cell-barcode read, streamed into gzipped FASTQs (sra-split). NCBI
// keeps the 10x barcode read as a technical read; ENA's FASTQs drop it, so there is no fallback
// to them as DOWNLOAD_SRA has. Fetching and converting in one task waits in the queue once, and
// the .sra file is deleted once converted rather than kept in work/.
process DOWNLOAD_SCATAC_READS {
    // As DOWNLOAD_SRA: AWS's copy takes many connections; NCBI's locator is asked once per run.
    // Queued jobs count too, and on Perlmutter most of a task's time is in the queue. Each task
    // holds at most its .sra file and compressed FASTQs (up to about 150 GB for an 80 GB run).
    maxForks 32
    // fastq-dump slices and their pigz use them all.
    cpus 8
    // Each fastq-dump slice takes about 100-200 MB.
    memory { oomMemoryOf(task, 4.GB) }
    time { retryTimeOf(task, sraReadsTime(run)) }
    // Retry an unavailable run (exit 75: a failed or stalled transfer, a corrupt file) with the
    // same resources, and one that ran out of memory or time with twice that
    // (modules/common/retry.nf). Any other failure (e.g. a run without the expected reads, which
    // converting again cannot fix) stops the run.
    errorStrategy { task.exitStatus == 75 ? 'retry' : retryErrorStrategy(task) }
    maxRetries 2
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${meta.id}:${run}"
    input:
    tuple val(meta), val(run)
    output:
    tuple val(meta), val(run), path("${run}.R1.fastq.gz"), path("${run}.R2.fastq.gz"),
        path("${run}.barcode.fastq.gz"), emit: reads
    script:
    // Reads are told apart by length, as submitters order them differently (e.g. 8/50/16/49,
    // 50/16/50, 51/51/16): the barcode read is meta.barcode_length long, the genomic mates are the
    // two longest, and an 8-base sample index read is dropped. Each read's length is its longest
    // among the first 1,000 spots, so a trimmed first read does not decide it.
    """
    set -o pipefail
    # sra-tools keeps its settings outside the home folder.
    export NCBI_SETTINGS="\$PWD/ncbi-settings.mkfg"
    # sra-fetch and sra-split are this repository's workflows/bin/ scripts, not the image's.
    status=0
    reason=\$(sra-fetch '${params.sra_locate_url}' '${run}') || status=\$?
    if (( status == 2 )); then
        printf '%s: SRA run unavailable (%s)\\n' '${run}' "\$reason" >&2
        exit 75
    elif (( status != 0 )); then
        exit "\$status"
    fi
    # A run loaded from BAM stores aligned reads, and no barcode read.
    if vdb-dump --table_enum '${run}.sra' | grep -q '^tbl *: *PRIMARY_ALIGNMENT\$'; then
        printf '%s stores aligned reads (loaded from BAM), with no barcode read\\n' '${run}' >&2
        exit 1
    fi
    # reads.tsv: read number, longest and shortest read among the first 1,000 spots, by longest.
    # Headers end in /READ, as sra-split numbers the reads.
    fastq-dump -X 1000 --split-spot -Z --defline-seq '@\$ac.\$si \$sn/\$ri' --defline-qual '+' \\
        '${run}.sra' | awk -v OFS='\\t' '
        NR % 4 == 1 {read = \$0; sub(/.*\\//, "", read)}
        NR % 4 == 2 {
            n = length(\$0)
            if (!(read in max) || n > max[read]) max[read] = n
            if (!(read in min) || n < min[read]) min[read] = n
        }
        END {for (r in max) print r, max[r], min[r]}' | sort -t \$'\\t' -k2,2n -k1,1n > reads.tsv
    cat reads.tsv >&2
    barcode=\$(awk -F '\\t' '\$2 == ${meta.barcode_length} {print \$1}' reads.tsv)
    if [ "\$(printf '%s' "\$barcode" | grep -c .)" != 1 ]; then
        printf '%s has no single %s-base barcode read; is the chemistry right?\\n' '${run}' \\
            '${meta.barcode_length}' >&2
        exit 1
    fi
    # Trimmed barcode reads would give some cells shorter barcodes.
    shortest=\$(awk -F '\\t' -v read="\$barcode" '\$1 == read {print \$3}' reads.tsv)
    if [ "\$shortest" != ${meta.barcode_length} ]; then
        printf '%s has barcode reads of %s to %s bases; expected all %s (trimmed reads?)\\n' \\
            '${run}' "\$shortest" '${meta.barcode_length}' '${meta.barcode_length}' >&2
        exit 1
    fi
    mates=\$(awk -F '\\t' '\$2 >= 30' reads.tsv | tail -n 2 | cut -f1 | sort -n)
    if [ "\$(printf '%s\\n' "\$mates" | grep -c .)" != 2 ]; then
        printf '%s has no two genomic reads of 30 bases or more\\n' '${run}' >&2
        exit 1
    fi
    # The lower-numbered genomic read is R1.
    first=\$(printf '%s\\n' "\$mates" | head -n 1)
    second=\$(printf '%s\\n' "\$mates" | tail -n 1)
    sra-split '${run}.sra' ${task.cpus} 1 "\$first=${run}.R1.fastq.gz" \\
        "\$second=${run}.R2.fastq.gz" "\$barcode=${run}.barcode.fastq.gz" > counts.tsv
    rm '${run}.sra'
    cat counts.tsv >&2
    # Every spot must have each of the three reads.
    awk -F '\\t' -v reads="\$first \$second \$barcode" -v run='${run}' '
        \$1 == "spots" {spots = \$2; next}
        {records[\$1] = \$2}
        END {
            n = split(reads, chosen, " ")
            for (i = 1; i <= n; i++) if (records[chosen[i]] + 0 != spots) {
                printf "%s: read %s in %d of %d spots\\n", run, chosen[i], records[chosen[i]], spots > "/dev/stderr"
                bad = 1
            }
            exit bad || !spots
        }' counts.tsv
    rm reads.tsv counts.tsv
    """
    stub:
    """
    printf '@r1\\nACGT\\n+\\nIIII\\n' | gzip > '${run}.R1.fastq.gz'
    printf '@r1\\nACGT\\n+\\nIIII\\n' | gzip > '${run}.R2.fastq.gz'
    printf '@r1\\nACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIII\\n' | gzip > '${run}.barcode.fastq.gz'
    """
}
