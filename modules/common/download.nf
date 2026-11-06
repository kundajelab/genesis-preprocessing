include { dotenv } from 'plugin/nf-dotenv'

// The size in bytes the server reports for a URL (HTTP HEAD), or -1 if it reports none.
def contentLength(url) {
    try {
        def connection = URI.create(url).toURL().openConnection() as java.net.HttpURLConnection
        connection.requestMethod = 'HEAD'
        connection.connectTimeout = 10000
        connection.readTimeout = 10000
        def length = connection.responseCode == 200 ? connection.contentLengthLong : -1L
        connection.disconnect()
        return length
    } catch (Exception _error) {
        return -1L
    }
}

// Time for a DOWNLOAD task: its FASTQs' sizes at 0.5 MB/s, below the slowest rate seen from ENA
// (0.67 MB/s for a whole file), plus 30 minutes for converting and counting reads. SRA runs come
// from AWS far faster (about 40 MB/s from Perlmutter), but a sample falls back to ENA when that
// fails. Where a server reports no size, a generous 15 hours per sample. The sizes are looked up
// here, on the manager, rather than passed in, so they do not change any task's cache key.
def downloadTime(samples, attempt) {
    def urls = samples.collectMany { meta ->
        meta.layout == 'PE' ? [meta.read1_url, meta.read2_url] : [meta.read1_url]
    }
    def sizes = urls.collect { url -> contentLength(url) }
    if (sizes.any { size -> size < 0 }) {
        return 15.h * samples.size() * attempt
    }
    def seconds = 1800L + (sizes.sum() / 0.5e6) as long
    return "${seconds * attempt}s"
}

// The SRA run whose reads a sample's URLs are, when they are ENA's FASTQs of one whole run:
// .../RUN/RUN.fastq.gz for single-end, .../RUN/RUN_1.fastq.gz and .../RUN/RUN_2.fastq.gz for
// paired-end. Otherwise null, and the sample is fetched from its URLs as given.
def sraRun(meta) {
    def first = meta.read1_url =~ /\/([SED]RR[0-9]+)\/([SED]RR[0-9]+)(_1)?\.fastq\.gz$/
    if (!first.find() || first.group(1) != first.group(2)) {
        return null
    }
    def run = first.group(1)
    if (meta.layout == 'PE') {
        def second = meta.read2_url.endsWith("/${run}/${run}_2.fastq.gz")
        return first.group(3) && second ? run : null
    }
    return first.group(3) ? null : run
}

process DOWNLOAD {
    // Each task holds a sample's uncompressed reads while it converts them, and a sample falls
    // back to ENA, which caps each connection and refuses (HTTP 403) clients with too many.
    maxForks 8
    // A transfer can end early or arrive corrupt (the checks below catch it); fetch the batch
    // again rather than stop the whole run.
    errorStrategy 'retry'
    maxRetries 2
    // fasterq-dump and pigz use them all.
    cpus 4
    memory '2 GB'
    time { downloadTime(samples, task.attempt) }
    conda "${moduleDir}/../../environments/GENESIS_SRA.yaml"
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${samples.size() == 1 ? samples[0].id : "${samples[0].id} +${samples.size() - 1} more"}"
    publishDir "${params.publish_base}/output", mode: 'copy', pattern: '*.fastq.lines',
        saveAs: { name ->
            def meta = samples.find { sample -> "${sample.id}.fastq.lines" == name }
            "${meta.species}/${meta.id}/${name}"
        }
    input:
    val samples
    output:
    tuple val(samples), path('*.read*.fastq.gz'), path('*.fastq.lines'), emit: downloads
    script:
    // One line per sample: ID, SRA run (or -), number of mates, then its URLs.
    def plan = samples.collect { meta ->
        def urls = meta.layout == 'PE' ? [meta.read1_url, meta.read2_url] : [meta.read1_url]
        ([meta.id, sraRun(meta) ?: '-', urls.size()] + urls).join('\t')
    }.join('\n')
    """
    set -o pipefail
    # sra-tools keeps its settings outside the home folder.
    export NCBI_SETTINGS="\$PWD/ncbi-settings.mkfg"
    retrying_curl() {
        curl --fail --silent --show-error --location --retry 20 --retry-all-errors \\
            --retry-delay 30 --connect-timeout 60 "\$@"
    }

    # An SRA run's .sra file from AWS, checked against NCBI's MD5, and converted to FASTQs named and
    # formatted as ENA's (header "@RUN.N NAME/MATE", bare "+"): the same bytes, decompressed.
    # Fails, for the caller to fall back to the URLs, if any step does.
    fetch_sra() {
        local id=\$1 run=\$2 mates=\$3 located link md5 mate expected file
        located=\$(retrying_curl --retry 3 --retry-delay 10 '${params.sra_locate_url}'"\$run") || return 1
        # The first file of type sra: its MD5, then its first location's link.
        located=\$(grep -o '"type": *"sra".*' <<< "\$located") || return 1
        md5=\$(grep -o '"md5": *"[0-9a-f]\\{32\\}"' <<< "\$located" | head -n 1 | grep -o '[0-9a-f]\\{32\\}')
        link=\$(grep -o '"link": *"[^"]*"' <<< "\$located" | head -n 1 | sed 's/^"link": *"//; s/"\$//')
        [[ -n \$md5 && -n \$link ]] || return 1
        retrying_curl --continue-at - --output "\$run.sra" "\$link" || return 1
        [[ \$(md5sum < "\$run.sra" | cut -d ' ' -f 1) == "\$md5" ]] || {
            echo "\$run.sra: MD5 differs from NCBI's \$md5" >&2
            return 1
        }
        fasterq-dump --split-3 --threads ${task.cpus} --temp "tmp-\$run" --outdir "sra-\$run" \\
            --seq-defline '@\$ac.\$si \$sn/\$ri' --qual-defline '+' "\$run.sra" >&2 || return 1
        rm "\$run.sra"
        # Exactly the files ENA serves for this layout: RUN.fastq, or RUN_1 and RUN_2.
        if (( mates == 1 )); then expected="\$run.fastq"; else expected="\${run}_1.fastq \${run}_2.fastq"; fi
        [[ \$(cd "sra-\$run" && ls | sort | tr '\\n' ' ') == "\$expected " ]] || {
            echo "\$run: fasterq-dump wrote \$(ls "sra-\$run"), not \$expected" >&2
            return 1
        }
        for mate in \$(seq 1 "\$mates"); do
            if (( mates == 1 )); then file="sra-\$run/\$run.fastq"; else file="sra-\$run/\${run}_\$mate.fastq"; fi
            wc -l < "\$file" | tr -d ' ' > "\$id.read\$mate.lines" || return 1
            pigz -p ${task.cpus} -c "\$file" > "\$id.read\$mate.fastq.gz" || return 1
            rm "\$file"
        done
        rm -r "sra-\$run" "tmp-\$run"
    }

    # FASTQs from their URLs as given; decompressing them all also checks their gzip integrity.
    fetch_urls() {
        local id=\$1 mates=\$2 mate
        shift 2
        for mate in \$(seq 1 "\$mates"); do
            retrying_curl --continue-at - --output "\$id.read\$mate.fastq.gz" "\$1"
            pigz -dc "\$id.read\$mate.fastq.gz" | wc -l | tr -d ' ' > "\$id.read\$mate.lines"
            shift
        done
    }

    cat > plan.tsv << 'END_OF_PLAN'
${plan}
END_OF_PLAN
    # On its own descriptor, so that no command in the loop can read the plan as its input.
    while IFS=\$'\\t' read -r -u 4 id run mates urls; do
        read -r -a urls <<< "\$urls"
        if [[ \$run != - ]]; then
            if ! fetch_sra "\$id" "\$run" "\$mates"; then
                echo "WARNING: \$id: SRA run \$run failed; downloading its FASTQs instead" >&2
                rm -rf "\$run.sra" "sra-\$run" "tmp-\$run" "\$id".read*
                run=-
            fi
        fi
        if [[ \$run == - ]]; then
            fetch_urls "\$id" "\$mates" "\${urls[@]}"
        fi
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
