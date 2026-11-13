include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process CALL_PEAKS {
    // MACS3 bounds its threads and peak workers to --threads (upstream PR #787). On a
    // Sorghum sample (52M treatment + 12M control pairs) 4 CPUs ran 5.7x faster than 1, and 8
    // CPUs no faster than 4.
    cpus 4
    // On Perlmutter, 73 Sorghum samples (2.0-30.6 GB of treatment plus control BAMs, all with the
    // same 1.8 GB control) peaked at about 0.34 GB plus 0.25 per GB of BAMs, and all within 1.06 GB
    // plus 0.30 per GB; the largest peaked at 9.0 GB. Allow 1.5 GB plus 0.30 per GB.
    memory { oomMemoryOf(task, 1.5.GB + 1.B * (0.30 * (bam.size() + control_bam.size())).toLong()) }
    // They took 8-78 seconds, at most 8.6 seconds per GB of BAMs. Allow 10 minutes plus 30 seconds
    // per GB.
    time { retryTimeOf(task, 10.min + 30.s * (bam.size() + control_bam.size()) / 1e9) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_MACS3.yaml"
    container "${dotenv('GENESIS_MACS3_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
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
