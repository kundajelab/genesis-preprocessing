include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomErrorStrategy } from '../common/retry'

process ATAC_TRIM {
    cpus 4
    memory { oomMemoryOf(task, 2.GB) }
    time { 1.h * (1 + reads.sum { read -> read.size() } / 1e9) * task.attempt }
    // Retry with twice the memory only when an attempt ran out of memory (modules/common/retry.nf).
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_CUTADAPT.yaml"
    container "${dotenv('GENESIS_CUTADAPT_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}/qc" }, mode: 'copy',
        pattern: '*.cutadapt.log'
    input:
    tuple val(meta), path(reads)
    val adapter
    output:
    tuple val(meta), path("${meta.id}.trimmed.read{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.cutadapt.log"), emit: report
    script:
    def files = reads.sort { read -> read.name }
    // Tn5 inserts the same Nextera adapter at both ends. ENCODE ATAC drops pairs where either mate
    // is shorter than 5 bases after trimming.
    """
    cutadapt -j ${task.cpus} -a '${adapter}' -A '${adapter}' --minimum-length 5 \\
        -o '${meta.id}.trimmed.read1.fastq.gz' -p '${meta.id}.trimmed.read2.fastq.gz' \\
        '${files[0]}' '${files[1]}' > '${meta.id}.cutadapt.log'
    """
    stub:
    def files = reads.sort { read -> read.name }
    """
    cp '${files[0]}' '${meta.id}.trimmed.read1.fastq.gz'
    cp '${files[1]}' '${meta.id}.trimmed.read2.fastq.gz'
    printf 'stub\\n' > '${meta.id}.cutadapt.log'
    """
}
