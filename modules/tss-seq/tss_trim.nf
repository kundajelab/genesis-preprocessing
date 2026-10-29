include { dotenv } from 'plugin/nf-dotenv'

process TSS_TRIM {
    cpus 4
    memory '2 GB'
    time { 1.h * (1 + reads.size() / 1e9) * task.attempt }
    // retry with more time when killed (walltime or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_CUTADAPT.yaml"
    container "${dotenv('GENESIS_CUTADAPT_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.cutadapt.log'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.id}.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.cutadapt.log"), emit: log
    script:
    // The read's 5' end is the TSS: remove only the sheet's fixed barcode or linker there, and
    // quality-trim the 3' end alone.
    def cut = meta.trim_5p == '-' ? '' : "--cut ${meta.trim_5p}"
    def adapter = meta.adapter_3p == '-' ? '' : "--adapter ${meta.adapter_3p}"
    """
    cutadapt --cores ${task.cpus} ${cut} ${adapter} --quality-cutoff 0,20 \\
        --minimum-length ${meta.min_length} -o '${meta.id}.trimmed.fastq.gz' '${reads}' \\
        > '${meta.id}.cutadapt.log'
    """
    stub:
    """
    printf '@r\\nACGT\\n+\\nIIII\\n' | gzip > '${meta.id}.trimmed.fastq.gz'
    printf 'Total reads processed: 1\\nReads written (passing filters): 1 (100.0%%)\\n' \\
        > '${meta.id}.cutadapt.log'
    """
}
