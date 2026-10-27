include { dotenv } from 'plugin/nf-dotenv'

process RNA_COUNTS {
    cpus 4
    memory '2 GB'
    time { 1.h * (1 + bam.size() / 1e9) }
    conda "${moduleDir}/../../environments/GENESIS_RNA_COUNTS.yaml"
    container "${dotenv('GENESIS_RNA_COUNTS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.featureCounts.tsv'
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
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
