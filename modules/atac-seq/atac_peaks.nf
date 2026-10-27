include { dotenv } from 'plugin/nf-dotenv'

process ATAC_PEAKS {
    // As CALL_PEAKS: MACS3 bounds its threads and peak workers to --threads.
    cpus 4
    memory { (1.GB + 1.B * (0.25 * bam.size()).toLong()) * task.attempt }
    time { 1.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_MACS3.yaml"
    container "${dotenv('GENESIS_MACS3_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
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
