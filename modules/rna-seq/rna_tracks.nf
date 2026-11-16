include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process RNA_TRACKS {
    // bamCoverage kept 2 CPUs busy; as for ATAC and DAP-seq tracks, it runs 1.7x as fast on 4.
    cpus 4
    // On Perlmutter with 2 CPUs, 20 Sorghum BAMs peaked at 0.78 GB; bamCoverage's workers each
    // hold a chunk of the genome.
    memory { oomMemoryOf(task, 1.GB + 512.MB * task.cpus) }
    // On Perlmutter with 2 CPUs, 20 Sorghum BAMs (0.17-24 GB; 0.72 Gbp genome) took 4 minutes for
    // the smallest and at most 4.9 minutes per GB (1.9 hours for 24 GB). Allow 10 minutes plus
    // 0.25 hours per Gbp and 0.1 per GB of BAM, at 2 CPUs' speed: 4 CPUs' speed-up is margin.
    time {
        def hours = (0.25 * meta.genome_size + 0.1 * bam.size()) / 1e9
        retryTimeOf(task, 10.min + 1.h * hours)
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
    tuple val(meta), path("${meta.id}.cpm.bw"), emit: cpm
    script:
    // Uniquely aligned reads (STAR gives them MAPQ 255), per base and in CPM; spliced reads cover
    // their exons only. Both strands together.
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    # bamCoverage's temporary files go in the task's folder, not a /tmp that may be held in memory.
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$MPLCONFIGDIR" "\$TMPDIR"
    bamCoverage -b '${bam}' -o '${meta.id}.cpm.bw' --binSize 1 --minMappingQuality 255 \\
        --normalizeUsing CPM -p ${task.cpus}
    """
    stub:
    """
    touch '${meta.id}.cpm.bw'
    """
}
