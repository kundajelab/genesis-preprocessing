include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process ATAC_TRACKS {
    // As DAP-seq's TRACKS: bamCoverage runs 1.7x as fast on 4 CPUs as on 2.
    cpus 4
    // DAP-seq's bamCoverage processes peaked at 1.58 GB on 4 CPUs whatever the BAM size, with their
    // temporary files in the task's folder; samtools sort uses 512 MB per thread.
    memory { oomMemoryOf(task, 1.GB + 512.MB * task.cpus) }
    // Not yet measured for ATAC-seq. DAP-seq's 1-bp tracks took about 0.15 hours per Gbp of genome
    // plus 0.013 per GB of BAM for each bamCoverage run on 4 CPUs; this has alignmentSieve, a sort
    // and one bamCoverage run. Allow 30 minutes plus 1 hour per Gbp and 0.1 hours per GB of BAM.
    time {
        retryTimeOf(task, (30.min + 1.h * (meta.genome_size + 0.1 * bam.size()) / 1e9) * 4 / task.cpus)
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
