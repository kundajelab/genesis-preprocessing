include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process TSS_TRIM {
    cpus 4
    // On Perlmutter, 24 Arabidopsis libraries (0.5-9.6 GB) peaked at 0.24 GB and took at most
    // 0.019 hours per GB (5.5 minutes for the largest). Allow 10 minutes plus 0.05 hours per GB.
    memory { oomMemoryOf(task, 2.GB) }
    time { retryTimeOf(task, 10.min + 1.h * (0.05 * reads.size() / 1e9)) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    container "${dotenv('GENESIS_CUTADAPT_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}/qc" }, mode: 'copy',
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
