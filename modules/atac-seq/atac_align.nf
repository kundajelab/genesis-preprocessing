include { dotenv } from 'plugin/nf-dotenv'

process ATAC_ALIGN {
    cpus 8
    // The loaded index is ~5.5 bytes per reference base; samtools sort uses 512 MB, and each bwa-mem2
    // thread adds a little more.
    memory { (1.5.GB + 7.B * meta.genome_size + 256.MB * task.cpus) * task.attempt }
    time { 1.h * (1 + reads.sum { read -> read.size() } / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.marked.*.txt'
    input:
    tuple val(meta), path(reads), path(index)
    val bwa_options
    output:
    tuple val(meta), path("${meta.id}.marked.bam"), path("${meta.id}.marked.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.marked.*.txt"), emit: reports
    script:
    def files = reads.sort { it.name }
    def threads = Math.max(1, task.cpus - 1)
    def options = "-K ${bwa_options.batch_size} -k ${bwa_options.seed_length} " +
        "-c ${bwa_options.max_seed_occurrences} -T ${bwa_options.score_threshold}"
    // bwa-mem2 writes mates together, as fixmate needs. Duplicates are marked here and removed with
    // the other filters, so the complexity metrics can still count them.
    """
    bwa-mem2 mem -t ${threads} ${options} -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA' \\
        "${index}/genome" '${files[0]}' '${files[1]}' 2> bwa.log | \\
        samtools fixmate -m -u - - | \\
        samtools sort -u -m 512M - | \\
        samtools markdup -f '${meta.id}.marked.markdup.txt' - '${meta.id}.marked.bam'
    samtools index '${meta.id}.marked.bam'
    samtools flagstat '${meta.id}.marked.bam' > '${meta.id}.marked.flagstat.txt'
    samtools idxstats '${meta.id}.marked.bam' > '${meta.id}.marked.idxstats.txt'
    samtools stats '${meta.id}.marked.bam' > '${meta.id}.marked.stats.txt'
    """
    stub:
    """
    touch '${meta.id}.marked.bam' '${meta.id}.marked.bam.bai'
    for report in markdup flagstat idxstats stats; do
        printf 'stub\\n' > '${meta.id}.marked.'"\$report"'.txt'
    done
    """
}
