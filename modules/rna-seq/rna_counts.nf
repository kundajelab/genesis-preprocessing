include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process RNA_COUNTS {
    cpus 4
    // On Perlmutter, 20 Sorghum BAMs (0.17-24 GB) peaked at 1.4 GB, at most 2.9 minutes per GB of
    // BAM (43 minutes for 14.7 GB; multimapping libraries are the slowest, not the largest).
    // Allow 10 minutes plus 0.1 hours per GB.
    memory { oomMemoryOf(task, 2.GB) }
    time { retryTimeOf(task, 10.min + 0.1.h * bam.size() / 1e9) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_RNA_COUNTS.yaml"
    container "${dotenv('GENESIS_RNA_COUNTS_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
        pattern: '*.featureCounts.tsv'
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}/qc" }, mode: 'copy',
        pattern: '*.featureCounts.tsv.summary'
    input:
    tuple val(meta), path(bam), path(bai), path(gtf)
    output:
    tuple val(meta), path("${meta.id}.featureCounts.tsv"), emit: counts
    tuple val(meta), path("${meta.id}.featureCounts.tsv.summary"), emit: summary
    script:
    def strand = [unstranded: 0, forward: 1, reverse: 2][meta.strandedness]
    // Pairs count once. Reads on several loci (NH > 1) or overlapping several genes are left
    // out, as featureCounts does by default; the summary records how many.
    def pairs = meta.layout == 'PE' ? '-p --countReadPairs' : ''
    """
    featureCounts -T ${task.cpus} -a '${gtf}' -t exon -g gene_id -s ${strand} ${pairs} \\
        -o '${meta.id}.featureCounts.tsv' '${bam}'
    """
    stub:
    """
    printf 'Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${meta.id}.bam\\n' > '${meta.id}.featureCounts.tsv'
    printf 'Status\\t${meta.id}.bam\\nAssigned\\t0\\n' > '${meta.id}.featureCounts.tsv.summary'
    """
}
