include { dotenv } from 'plugin/nf-dotenv'
include { readUrls ; urlTime } from './download_support'

// A sample's FASTQs from its URLs as given; decompressing them also checks their gzip integrity.
// For samples with no SRA run, and those whose run DOWNLOAD_SRA could not use.
process DOWNLOAD_URLS {
    // ENA caps each connection and refuses (HTTP 403) clients with too many.
    maxForks 8
    // A transfer can end early or arrive corrupt (the checks below catch it); fetch the sample
    // again rather than stop the whole run.
    errorStrategy 'retry'
    maxRetries 2
    // curl, and pigz decompressing to check the files.
    cpus 2
    memory '1 GB'
    time { urlTime(meta, task.attempt) }
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.publish_base}/output", mode: 'copy', pattern: '*.fastq.lines',
        saveAs: { name -> "${meta.species}/${meta.id}/${name}" }
    input:
    val meta
    output:
    tuple val(meta), path("${meta.id}.read*.fastq.gz", arity: '1..2'), emit: downloads
    path "${meta.id}.fastq.lines", emit: counts
    script:
    def id = meta.id
    def urls = readUrls(meta)
    """
    set -o pipefail
    retrying_curl() {
        curl --fail --silent --show-error --location --retry 20 --retry-all-errors \\
            --retry-delay 30 --connect-timeout 60 "\$@"
    }

    # The mates' URLs, read from a here-document so that no URL is parsed as shell code.
    read -r -a urls << 'END_OF_URLS'
${urls.join(' ')}
END_OF_URLS
    for mate in \$(seq 1 ${urls.size()}); do
        retrying_curl --continue-at - --output "${id}.read\$mate.fastq.gz" "\${urls[mate - 1]}"
        pigz -dc "${id}.read\$mate.fastq.gz" | wc -l | tr -d ' ' > "${id}.read\$mate.lines"
    done
    # Whole FASTQ records, at least one, and as many in each mate.
    awk '{if (\$1 == 0 || \$1 % 4) exit 1; print \$1}' '${id}.read1.lines' > '${id}.fastq.lines'
    if (( ${urls.size()} == 2 )); then
        cmp -s '${id}.read1.lines' '${id}.read2.lines'
    fi
    rm '${id}'.read*.lines
    """
    stub:
    def second = meta.layout == 'PE' ? "cp '${meta.id}.read1.fastq.gz' '${meta.id}.read2.fastq.gz'" : ''
    """
    printf '@read1\\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip > '${meta.id}.read1.fastq.gz'
    ${second}
    printf '4\\n' > '${meta.id}.fastq.lines'
    """
}
