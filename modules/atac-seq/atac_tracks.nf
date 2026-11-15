include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process ATAC_TRACKS {
    // As DAP-seq's TRACKS: bamCoverage runs 1.7x as fast on 4 CPUs as on 2.
    cpus 4
    // On Perlmutter, 24 Sorghum samples peaked at 2.9 GB of 3, whatever the BAM size: samtools
    // sort's 512 MB per thread, plus bamCoverage. Allow 1.5 GB more than the sort's buffers.
    memory { oomMemoryOf(task, 1.5.GB + 512.MB * task.cpus) }
    // On Perlmutter with 4 CPUs, 24 Sorghum BAMs (0.13-11.9 GB; 0.73 Gbp genome) took at most 0.14
    // hours per GB (50 minutes for 6.0 GB; 64 for the largest). Allow 10 minutes plus 0.25 hours
    // per Gbp and 0.2 per GB of BAM, at 4 CPUs' speed.
    time {
        def hours = (0.25 * meta.genome_size + 0.2 * bam.size()) / 1e9
        retryTimeOf(task, (10.min + 1.h * hours) * 4 / task.cpus)
    }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_DEEPTOOLS.yaml"
    container "${dotenv('GENESIS_DEEPTOOLS_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
        pattern: '*.bw'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path("${meta.id}.cuts.counts.bw"), emit: cuts
    script:
    // Tn5 cut sites: each read's 5' base after the usual +4/-5 shift (alignmentSieve --ATACshift),
    // i.e. BED0 start + 4 on the forward strand and end - 6 on the reverse. Raw counts per base.
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    # bamCoverage's temporary files go in the task's folder, not a /tmp that may be held in memory.
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$MPLCONFIGDIR" "\$TMPDIR"
    alignmentSieve -b '${bam}' -o shifted.unsorted.bam --ATACshift -p ${task.cpus}
    samtools sort -@ ${task.cpus} -m 512M -o shifted.bam shifted.unsorted.bam
    samtools index -@ ${task.cpus} shifted.bam
    bamCoverage -b shifted.bam -o '${meta.id}.cuts.counts.bw' --binSize 1 --Offset 1 \\
        --normalizeUsing None -p ${task.cpus}
    """
    stub:
    """
    touch '${meta.id}.cuts.counts.bw'
    """
}
