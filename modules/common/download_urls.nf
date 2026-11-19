include { dotenv } from 'plugin/nf-dotenv'
include { batchTag ; countsPath ; readUrls ; urlTime } from './download_support'

// Each sample's FASTQs from their URLs as given; decompressing them all also checks their gzip
// integrity. For samples with no SRA run, and those whose run DOWNLOAD_SRA could not use.
process DOWNLOAD_URLS {
    // ENA caps each connection and refuses (HTTP 403) clients with too many.
    maxForks 8
    // A transfer can end early or arrive corrupt (the checks below catch it); fetch the batch
    // again rather than stop the whole run.
    errorStrategy 'retry'
    maxRetries 2
    // curl, and pigz decompressing to check the files.
    cpus 2
    memory '1 GB'
    time { urlTime(samples, task.attempt) }
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${batchTag(samples)}"
    publishDir "${params.publish_base}/output", mode: 'copy', pattern: '*.fastq.lines',
        saveAs: { name -> countsPath(samples, name) }
    input:
    val samples
    output:
    tuple val(samples), path('*.read*.fastq.gz'), emit: downloads
    path '*.fastq.lines', emit: counts
    script:
    // One line per sample: ID, number of mates, then its URLs.
    def plan = samples.collect { meta -> ([meta.id, readUrls(meta).size()] + readUrls(meta)).join('\t') }
        .join('\n')
    """
    set -o pipefail
    retrying_curl() {
        curl --fail --silent --show-error --location --retry 20 --retry-all-errors \\
            --retry-delay 30 --connect-timeout 60 "\$@"
    }

    cat > plan.tsv << 'END_OF_PLAN'
${plan}
END_OF_PLAN
    # On its own descriptor, so that no command in the loop can read the plan as its input.
    while IFS=\$'\\t' read -r -u 4 id mates urls; do
        read -r -a urls <<< "\$urls"
        for mate in \$(seq 1 "\$mates"); do
            retrying_curl --continue-at - --output "\$id.read\$mate.fastq.gz" "\${urls[mate - 1]}"
            pigz -dc "\$id.read\$mate.fastq.gz" | wc -l | tr -d ' ' > "\$id.read\$mate.lines"
        done
        # Whole FASTQ records, at least one, and as many in each mate.
        awk '{if (\$1 == 0 || \$1 % 4) exit 1; print \$1}' "\$id.read1.lines" > "\$id.fastq.lines"
        if (( mates == 2 )); then
            cmp -s "\$id.read1.lines" "\$id.read2.lines"
        fi
        rm "\$id".read*.lines
    done 4< plan.tsv
    """
    stub:
    def fakes = samples.collect { meta ->
        def second = meta.layout == 'PE' ? "\ncp '${meta.id}.read1.fastq.gz' '${meta.id}.read2.fastq.gz'" : ''
        "printf '@read1\\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip > '${meta.id}.read1.fastq.gz'\n" +
            "printf '4\\n' > '${meta.id}.fastq.lines'" + second
    }.join('\n')
    """
${fakes}
    """
}
