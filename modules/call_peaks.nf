include { dotenv } from 'plugin/nf-dotenv'

process CALL_PEAKS {
    // MACS3 bounds its threads and peak workers to --threads (upstream PR #787). On a
    // Sorghum sample (52M treatment + 12M control pairs) 4 CPUs ran 5.7x faster than 1, and 8
    // CPUs no faster than 4.
    cpus 4
    // That sample needed ~1.1 GB at 1 CPU and ran within 2 GB at 4 CPUs, for 8.1 GB of BAMs.
    memory { (1.GB + 1.B * (0.25 * (bam.size() + control_bam.size())).toLong()) * task.attempt }
    time { 1.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../environments/DAP_SEQ_PEAKS.yaml"
    container "${dotenv('DAP_SEQ_PEAKS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.macs3*'
    input:
    tuple val(meta), path(bam), path(bai), path(control_bam), path(control_bai)
    output:
    tuple val(meta), path("${meta.id}.macs3_peaks.narrowPeak.gz"), emit: peaks
    tuple val(meta), path("${meta.id}.macs3*"), emit: results
    script:
    def format = meta.layout == 'PE' ? 'BAMPE' : 'BAM'
    """
    macs3 callpeak -t '${bam}' -c '${control_bam}' -f ${format} -g ${meta.genome_size} \\
        -n '${meta.id}.macs3' --threads ${task.cpus}
    gzip -f '${meta.id}.macs3_peaks.narrowPeak'
    """
    stub:
    """
    printf 'chr1\\t100\\t200\\tpeak1\\t100\\t.\\t1\\t2\\t3\\t50\\n' | gzip > '${meta.id}.macs3_peaks.narrowPeak.gz'
    touch '${meta.id}.macs3_peaks.xls'
    """
}
