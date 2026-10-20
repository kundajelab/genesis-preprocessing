include { dotenv } from 'plugin/nf-dotenv'

process DOWNLOAD {
    // Each task fetches a batch of samples with parallel aria2c connections, so only one task
    // runs at a time to keep the load on the remote servers bounded.
    maxForks 1
    cpus 2
    memory '1 GB'
    time { 1.h * samples.size() }
    conda "environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${samples.size() == 1 ? samples[0].id : "${samples[0].id} +${samples.size() - 1} more"}"
    publishDir "${params.run_dir}/output", mode: 'copy', pattern: '*.fastq.lines',
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
    def aria2cInput = reads.collect { url, name -> "${url}\n  out=${name}" }.join('\n')
    def checks = samples.collect { meta ->
        def count = "count_lines '${meta.id}.read1.fastq.gz' | " +
            "awk '{if (\$1 == 0 || \$1 % 4) exit 1; print \$1}' > '${meta.id}.fastq.lines'"
        meta.layout == 'PE' ?
            count + "\ntest \"\$(count_lines '${meta.id}.read2.fastq.gz')\" -eq \"\$(cat '${meta.id}.fastq.lines')\"" :
            count
    }.join('\n')
    """
    aria2c_input_file=aria2c-input.txt
    cat > "\${aria2c_input_file}" << 'END_OF_URLS'
${aria2cInput}
END_OF_URLS
    # Batches usually hit a single host, which refuses connections beyond a per-client cap, so keep
    # the total (concurrent downloads x connections) modest and back off on refusals.
    aria2c \\
        --input-file="\${aria2c_input_file}" \\
        --max-connection-per-server=4 \\
        --split=4 \\
        --max-concurrent-downloads=8 \\
        --retry-wait=5 \\
        --rlimit-nofile \$(ulimit -n)
    # Decompressing the whole file also checks its gzip integrity.
    count_lines() {
        bgzip --decompress --stdout --threads ${task.cpus} "\$1" | wc -l | tr -d ' '
    }
${checks}
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
