include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process ATAC_PEAKS {
    // As CALL_PEAKS: MACS3 bounds its threads and peak workers to --threads.
    cpus 4
    // As CALL_PEAKS, from 73 DAP-seq samples (all within 1.06 GB plus 0.30 per GB of BAMs, and 8.6
    // seconds per GB). ATAC-seq calls more peaks, with no control; retries cover more memory.
    memory { oomMemoryOf(task, 1.5.GB + 1.B * (0.30 * bam.size()).toLong()) }
    time { retryTimeOf(task, 10.min + 30.s * bam.size() / 1e9) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_MACS3.yaml"
    container "${dotenv('GENESIS_MACS3_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
        pattern: '*.macs3*'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path("${meta.id}.macs3_peaks.narrowPeak.gz"), emit: peaks
    tuple val(meta), path("${meta.id}.macs3*"), emit: results
    script:
    // BAMPE peaks over whole unshifted fragments. Duplicates were already removed, so keep the rest.
    """
    macs3 callpeak -t '${bam}' -f BAMPE -g ${meta.nuclear_size} -q 0.01 --keep-dup all \\
        -n '${meta.id}.macs3' --threads ${task.cpus}
    gzip -f '${meta.id}.macs3_peaks.narrowPeak'
    """
    stub:
    """
    printf 'chr1\\t100\\t200\\tpeak1\\t100\\t.\\t1\\t2\\t3\\t50\\n' | gzip > '${meta.id}.macs3_peaks.narrowPeak.gz'
    touch '${meta.id}.macs3_peaks.xls'
    """
}
