include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from './retry'
include { readUrls ; sraRun ; sraTime } from './download_support'

// A sample's SRA run's .sra file from AWS, checked against NCBI's MD5 and converted to FASTQs
// named and formatted as ENA's (header "@RUN.N NAME/MATE", bare "+"): the same bytes, decompressed.
// SRA_STATUS is "ok", or "fallback" with the reason in SRA_REASON when the run is unavailable (an
// HTTP failure, a stalled transfer, an MD5 mismatch), stores aligned reads (loaded from BAM; 1 of
// the 127 runs in the test sheets), or its converted reads are not the FASTQs the sample's URLs
// name (another layout, a spot missing a mate): its partial files are deleted, and the task still
// succeeds, for DOWNLOAD_URLS to fetch the sample instead. Any other failure (fastq-dump, pigz, a
// full disk) stops the task.
process DOWNLOAD_SRA {
    // AWS's copy takes many connections; NCBI's locator is asked once per run.
    maxForks 32
    // fastq-dump slices and their pigz use them all.
    cpus 8
    // Each fastq-dump slice takes about 100-200 MB.
    memory { oomMemoryOf(task, 4.GB) }
    time { retryTimeOf(task, sraTime(meta)) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.publish_base}/output", mode: 'copy', pattern: '*.fastq.lines',
        saveAs: { name -> "${meta.species}/${meta.id}/${name}" }
    input:
    val meta
    output:
    tuple val(meta), env('SRA_STATUS'), env('SRA_REASON'), path("${meta.id}.read*.fastq.gz", arity: '0..*'), emit: downloads
    path "${meta.id}.fastq.lines", optional: true, emit: counts
    script:
    def id = meta.id
    def run = sraRun(meta)
    def mates = readUrls(meta).size()
    def outputs = (1..mates).collect { mate -> "${mate}=${id}.read${mate}.fastq.gz" }.join(' ')
    """
    set -o pipefail
    # sra-tools keeps its settings outside the home folder.
    export NCBI_SETTINGS="\$PWD/ncbi-settings.mkfg"
    SRA_STATUS=ok
    # Why the run is unavailable: set before convert returns 2.
    SRA_REASON=

    # Fetch and convert the run. Returns 2 when it is unavailable, for the sample to be fetched
    # from its URLs instead, and 1 on any other failure.
    convert() {
        local status=0
        # sra-fetch and sra-split are this repository's workflows/bin/ scripts, not the image's.
        SRA_REASON=\$(sra-fetch '${params.sra_locate_url}' '${run}') || status=\$?
        (( status == 0 )) || return "\$status"
        # Biological reads only, as ENA serves them; no file but the compressed FASTQs.
        sra-split --skip-technical '${run}.sra' ${task.cpus} 1 ${outputs} > counts.tsv || status=\$?
        if (( status == 3 )); then
            SRA_REASON="its .sra file stores aligned reads, too slow to convert"
            return 2
        elif (( status != 0 )); then
            echo "${run}: sra-split failed" >&2
            return 1
        fi
        rm '${run}.sra'
        cat counts.tsv >&2
        # Exactly the reads ENA serves for this layout: read 1, or reads 1 and 2, in every spot.
        awk -F '\\t' -v mates=${mates} '
            \$1 == "spots" {spots = \$2; next}
            {n++; if (\$1 > mates || \$2 != spots) bad = 1}
            END {exit !(spots > 0 && n == mates && !bad)}' counts.tsv || {
            echo "${run}: reads per spot are not this ${meta.layout} sample's ${mates}" >&2
            SRA_REASON="its .sra file is not of the sample's layout"
            return 2
        }
        awk -F '\\t' '\$1 == 1 {print 4 * \$2}' counts.tsv > '${id}.fastq.lines'
        rm counts.tsv
    }

    status=0
    convert || status=\$?
    if (( status == 2 )); then
        echo "WARNING: ${id}: SRA run ${run} unavailable (\$SRA_REASON); its FASTQs will be downloaded instead" >&2
        rm -rf '${run}.sra' counts.tsv sra-split.* '${id}'.read* '${id}.fastq.lines'
        SRA_STATUS=fallback
    elif (( status != 0 )); then
        exit "\$status"
    fi
    """
    stub:
    def second = meta.layout == 'PE' ? "cp '${meta.id}.read1.fastq.gz' '${meta.id}.read2.fastq.gz'" : ''
    """
    printf '@read1\\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip > '${meta.id}.read1.fastq.gz'
    ${second}
    printf '4\\n' > '${meta.id}.fastq.lines'
    SRA_STATUS=ok
    SRA_REASON=
    """
}
