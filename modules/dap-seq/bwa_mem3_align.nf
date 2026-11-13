include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process BWA_MEM3_ALIGN {
    cpus 16
    // bwa-mem3 loads ~3.5 bytes per reference base (it reads no .0123); 7 leaves headroom. samtools
    // sort uses 512 MB per thread, and each bwa-mem3 thread adds a little more.
    memory { oomMemoryOf(task, 1.5.GB + 7.B * meta.genome_size + (512.MB + 256.MB) * task.cpus) }
    // On Perlmutter with 16 CPUs (7 samples of 1.3-4.7 GB of compressed reads), bwa-mem3's three
    // alignments took about a minute plus 3.2-6.2 minutes per GB, and 9.2 for the repeat-rich
    // AL2G32420 (25 minutes). Allow twice 2 minutes plus 9.2 per GB, for 16 CPUs.
    time {
        def gigabytes = [reads].flatten().sum { read -> read.size() } / 1e9
        retryTimeOf(task, (2.min + 9.2.min * gigabytes) * 2 * 16 / task.cpus)
    }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
        pattern: '*.alignment.tsv'
    input:
    tuple val(meta), path(reads), path(metadata), path(index)
    val bwa_options
    val spp_read_length
    output:
    tuple val(meta), path("${meta.id}.primary.bam"), path("${meta.id}.primary.bam.bai"),
        path("${meta.id}.qc.bam"), path("${meta.id}.qc.bam.bai"), emit: aligned
    tuple val(meta), path("${meta.id}.spp.bam"), emit: spp
    tuple val(meta), path("${meta.id}.alignment.tsv"), emit: provenance
    script:
    def files = reads instanceof List ? reads.sort { read -> read.name } : [reads]
    def inputs = files.collect { read -> "'${read}'" }.join(' ')
    def options = "-K ${bwa_options.batch_size} -k ${bwa_options.seed_length} " +
        "-c ${bwa_options.max_seed_occurrences} -T ${bwa_options.score_threshold}"
    """
    # main alignment
    bwa-mem3 mem \\
        -t ${task.cpus} \\
        ${options} \\
        -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA' \\
        "${index}/genome" \\
        ${inputs} \\
        2> main.bwa.log \\
    | samtools view -u -F 2308 - \\
    | samtools sort -@ ${task.cpus} --write-index -m 512M \\
        -o '${meta.id}.primary.bam##idx##${meta.id}.primary.bam.bai' -

    # QC uses read 1 aligned on its own. For a single-end sample that is the main alignment (same
    # reads, same options), so link it rather than align again; its read group is the main one.
    if [ ${files.size()} -eq 1 ]; then
        ln '${meta.id}.primary.bam' '${meta.id}.qc.bam'
        ln '${meta.id}.primary.bam.bai' '${meta.id}.qc.bam.bai'
    else
        bwa-mem3 mem \\
            -t ${task.cpus} \\
            ${options} \\
            -R '@RG\\tID:${meta.id}.qc\\tSM:${meta.id}\\tPL:ILLUMINA' \\
            "${index}/genome" \\
            '${files[0]}' \\
            2> qc.bwa.log \\
        | samtools view -u -F 2308 - \\
        | samtools sort -@ ${task.cpus} --write-index -m 512M \\
            -o '${meta.id}.qc.bam##idx##${meta.id}.qc.bam.bai' -
    fi

    # SPP gets read 1 cut to its first bases before alignment: with full-length reads as long as the
    # fragments, the read-length phantom peak hides the fragment peak and SPP cannot estimate it.
    bwa-mem3 mem \\
        -t ${task.cpus} \\
        ${options} \\
        -R '@RG\\tID:${meta.id}.spp\\tSM:${meta.id}\\tPL:ILLUMINA' \\
        "${index}/genome" \\
        <(\\
            gzip -dc '${files[0]}' \\
            | awk -v bases=${spp_read_length} 'NR % 2 == 0 {\$0 = substr(\$0, 1, bases)} 1'\\
        ) \\
        2> spp.bwa.log \\
    | samtools view -u -F 2308 - \\
    | samtools sort -@ ${task.cpus} -m 512M -o '${meta.id}.spp.bam' -

    printf 'sample_id\\taligner\\tversion\\tlayout\\tread_length\\tmapq_filter\\tflag_exclude\\tbwa_options\\tqc_reads\\tspp_reads\\n' > '${meta.id}.alignment.tsv'
    printf '%s\\tbwa-mem3\\t%s\\t%s\\t%s\\tnone\\t2308\\t%s\\tfull_read1_unpaired\\tread1_first_${spp_read_length}bp_unpaired\\n' \\
        '${meta.id}' "\$(bwa-mem3 version | sed -n 1p)" '${meta.layout}' '${meta.read_length}' '${options}' >> '${meta.id}.alignment.tsv'
    """
    stub:
    """
    touch '${meta.id}.primary.bam' '${meta.id}.primary.bam.bai' '${meta.id}.qc.bam' '${meta.id}.qc.bam.bai' \\
        '${meta.id}.spp.bam'
    printf 'stub\\n' > '${meta.id}.alignment.tsv'
    """
}
