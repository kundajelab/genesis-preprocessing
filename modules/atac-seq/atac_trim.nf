include { dotenv } from 'plugin/nf-dotenv'

process ATAC_TRIM {
    cpus 4
    memory '2 GB'
    time { 1.h * (1 + reads.sum { read -> read.size() } / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/ATAC_TRIMMING.yaml"
    container "${dotenv('ATAC_TRIMMING_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.cutadapt.log'
    input:
    tuple val(meta), path(reads)
    val adapter
    output:
    tuple val(meta), path("${meta.id}.trimmed.read{1,2}.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.cutadapt.log"), emit: report
    script:
    def files = reads.sort { it.name }
    // Tn5 inserts the same Nextera adapter at both ends. ENCODE ATAC drops pairs where either mate
    // is shorter than 5 bases after trimming.
    """
    cutadapt -j ${task.cpus} -a '${adapter}' -A '${adapter}' --minimum-length 5 \\
        -o '${meta.id}.trimmed.read1.fastq.gz' -p '${meta.id}.trimmed.read2.fastq.gz' \\
        '${files[0]}' '${files[1]}' > '${meta.id}.cutadapt.log'
    """
    stub:
    def files = reads.sort { it.name }
    """
    cp '${files[0]}' '${meta.id}.trimmed.read1.fastq.gz'
    cp '${files[1]}' '${meta.id}.trimmed.read2.fastq.gz'
    printf 'stub\\n' > '${meta.id}.cutadapt.log'
    """
}
