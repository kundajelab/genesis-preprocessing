include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomErrorStrategy } from '../common/retry'

process ATAC_ALIGN {
    cpus 8
    // bwa-mem3 loads ~3.5 bytes per reference base (it reads no .0123); 7 leaves headroom. samtools
    // sort uses 512 MB per thread, and each bwa-mem3 thread adds a little more.
    memory { oomMemoryOf(task, 1.5.GB + 7.B * meta.genome_size + (512.MB + 256.MB) * task.cpus) }
    time { 1.h * (1 + reads.sum { read -> read.size() } / 1e9) * task.attempt }
    // Retry with twice the memory only when an attempt ran out of memory (modules/common/retry.nf).
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}/qc" }, mode: 'copy',
        pattern: '*.marked.*.txt'
    input:
    tuple val(meta), path(reads), path(index)
    val bwa_options
    output:
    tuple val(meta), path("${meta.id}.marked.bam"), path("${meta.id}.marked.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.marked.*.txt"), emit: reports
    script:
    def files = reads.sort { read -> read.name }
    def options = "-K ${bwa_options.batch_size} -k ${bwa_options.seed_length} " +
        "-c ${bwa_options.max_seed_occurrences} -T ${bwa_options.score_threshold}"
    // bwa-mem3 writes mates together, as fixmate needs. Duplicates are marked here and removed with
    // the other filters, so the complexity metrics can still count them. Every stage gets all the
    // CPUs: the samtools stages need little while bwa-mem3 runs, and sort's merge and markdup come
    // after it.
    """
    bwa-mem3 mem -t ${task.cpus} ${options} -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA' \\
        "${index}/genome" '${files[0]}' '${files[1]}' 2> bwa.log | \\
        samtools fixmate -@ ${task.cpus} -m -u - - | \\
        samtools sort -@ ${task.cpus} -u -m 512M - | \\
        samtools markdup -@ ${task.cpus} -f '${meta.id}.marked.markdup.txt' - '${meta.id}.marked.bam'
    samtools index -@ ${task.cpus} '${meta.id}.marked.bam'
    samtools flagstat -@ ${task.cpus} '${meta.id}.marked.bam' > '${meta.id}.marked.flagstat.txt'
    samtools idxstats -@ ${task.cpus} '${meta.id}.marked.bam' > '${meta.id}.marked.idxstats.txt'
    samtools stats -@ ${task.cpus} '${meta.id}.marked.bam' > '${meta.id}.marked.stats.txt'
    """
    stub:
    """
    touch '${meta.id}.marked.bam' '${meta.id}.marked.bam.bai'
    for report in markdup flagstat idxstats stats; do
        printf 'stub\\n' > '${meta.id}.marked.'"\$report"'.txt'
    done
    """
}
