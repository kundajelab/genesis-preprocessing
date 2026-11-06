include { dotenv } from 'plugin/nf-dotenv'

process DOWNLOAD {
    // ENA caps each connection (about 1 MB/s from Perlmutter) and refuses (HTTP 403) clients with
    // too many, so throughput comes from many files in flight over few connections: 8 tasks, each
    // fetching one file at a time over 2 connections, 16 in all.
    maxForks 8
    // A transfer can end early or arrive corrupt (the MD5 and gzip checks below catch it); fetch
    // the batch again rather than stop the whole run.
    errorStrategy 'retry'
    maxRetries 2
    // One CPU per sample: their reads are counted in parallel, each by a single-threaded gzip
    // decompression, as ENA's files are plain gzip rather than BGZF.
    cpus { samples.size() }
    memory '1 GB'
    time { 1.h * samples.size() * task.attempt }
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
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
    def reads = samples.collectMany { meta ->
        def urls = meta.layout == 'PE' ? [meta.read1_url, meta.read2_url] : [meta.read1_url]
        urls.withIndex().collect { url, index -> [url, "${meta.id}.read${index + 1}.fastq.gz"] }
    }
    def downloads = reads.collect { url, name -> "${url}\t${name}" }.join('\n')
    // Each sample's reads are counted in the background; a PE sample's mates must agree.
    def checks = samples.collect { meta ->
        def count = "count_lines '${meta.id}.read1.fastq.gz' | " +
            "awk '{if (\$1 == 0 || \$1 % 4) exit 1; print \$1}' > '${meta.id}.fastq.lines'"
        def mates = meta.layout == 'PE' ?
            "\n    test \"\$(count_lines '${meta.id}.read2.fastq.gz')\" -eq \"\$(cat '${meta.id}.fastq.lines')\"" :
            ''
        "(\n    set -euo pipefail\n    ${count}${mates}\n) &"
    }.join('\n')
    """
    cat > downloads.tsv << 'END_OF_URLS'
${downloads}
END_OF_URLS
    # ENA publishes each file's MD5. Its FASTQ URLs end in RUN/FILE, so look up the run's files and
    # checksums; elsewhere (other archives, or a failed lookup) the gzip check below still applies.
    ena_md5() {
        local path=\${1#*://} run report query
        [[ \$path == ftp.sra.ebi.ac.uk/vol1/fastq/* ]] || return 0
        run=\$(basename -- "\$(dirname -- "\$path")")
        report="ena-\${run}.tsv"
        query="accession=\${run}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
        if [[ ! -s \$report ]]; then
            aria2c -q --max-tries=5 --retry-wait=10 --allow-overwrite=true --out="\$report" \\
                "https://www.ebi.ac.uk/ena/portal/api/filereport?\${query}" || return 0
        fi
        awk -F '\\t' -v path="\$path" '
            NR == 1 { for (i = 1; i <= NF; i++) column[\$i] = i; next }
            {
                n = split(\$column["fastq_ftp"], files, ";")
                split(\$column["fastq_md5"], sums, ";")
                for (i = 1; i <= n; i++) if (files[i] == path) print sums[i]
            }' "\$report"
    }
    while IFS=\$'\\t' read -r url name; do
        printf '%s\\n  out=%s\\n' "\$url" "\$name"
        md5=\$(ena_md5 "\$url")
        if [[ -n \$md5 ]]; then printf '  checksum=md5=%s\\n' "\$md5"; fi
    done < downloads.tsv > aria2c-input.txt
    # One file at a time over 2 connections: with maxForks, this bounds the connections to the
    # host (see above). Servers also drop connections now and then: keep retrying, resuming the
    # partial file. (No --lowest-speed-limit: a download slower than it is aborted, not retried.)
    aria2c \\
        --input-file=aria2c-input.txt \\
        --max-connection-per-server=2 \\
        --split=2 \\
        --max-concurrent-downloads=1 \\
        --max-tries=20 \\
        --retry-wait=30 \\
        --timeout=120 \\
        --rlimit-nofile \$(ulimit -n)
    # Decompressing the whole file also checks its gzip integrity.
    count_lines() {
        bgzip --decompress --stdout "\$1" | wc -l | tr -d ' '
    }
${checks}
    for pid in \$(jobs -p); do
        wait "\$pid"
    done
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
