include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

// MACS3 peaks of a treatment against its control. A single-end sample (MACS3's BAM format; BAMPE
// takes fragment sizes from the pairs) whose treatment has too few paired + and - strand peaks for
// MACS3 to build its fragment-size model gets no peaks: the task still succeeds, publishing only
// SAMPLE.failure-reason.txt with MACS3's message, and nothing downstream (QUANTIFY) runs for it.
// Any other MACS3 failure stops the run.
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
    container "${dotenv('GENESIS_MACS3_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
        pattern: '{*.macs3*,*.failure-reason.txt}'
    input:
    tuple val(meta), path(bam), path(bai), path(control_bam), path(control_bai)
    output:
    tuple val(meta), path("${meta.id}.macs3_peaks.narrowPeak.gz"), optional: true, emit: peaks
    tuple val(meta), path("${meta.id}.macs3*"), optional: true, emit: results
    tuple val(meta), path("${meta.id}.failure-reason.txt"), optional: true, emit: failure
    script:
    def format = meta.layout == 'PE' ? 'BAMPE' : 'BAM'
    """
    status=0
    macs3 callpeak -t '${bam}' -c '${control_bam}' -f ${format} -g ${meta.genome_size} \\
        -n '${meta.id}.macs3' --threads ${task.cpus} 2> callpeak.log || status=\$?
    cat callpeak.log >&2
    if (( status != 0 )); then
        # MACS3 warns, then exits 1:
        #   MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can
        #   only find 88! Please make your MFOLD range broader and try again. ...
        #   Process for pairing-model is terminated!
        grep -q 'Process for pairing-model is terminated!' callpeak.log || exit "\$status"
        reason=\$(grep -o -m 1 'MACS3 needs at least [0-9]* paired peaks.*can only find [0-9]*!' \\
            callpeak.log || echo 'MACS3 could not build its peak model')
        rm -f '${meta.id}'.macs3*
        printf '%s\\n' "\$reason" > '${meta.id}.failure-reason.txt'
        printf 'WARNING: ${meta.id}: no peaks called: %s\\n' "\$reason" >&2
    else
        gzip -f '${meta.id}.macs3_peaks.narrowPeak'
    fi
    rm callpeak.log
    """
    stub:
    """
    printf 'chr1\\t100\\t200\\tpeak1\\t100\\t.\\t1\\t2\\t3\\t50\\n' | gzip > '${meta.id}.macs3_peaks.narrowPeak.gz'
    touch '${meta.id}.macs3_peaks.xls'
    """
}
