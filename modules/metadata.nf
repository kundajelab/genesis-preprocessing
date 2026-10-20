include { dotenv } from 'plugin/nf-dotenv'

process METADATA {
    cpus 1
    memory '1 GB'
    time '1h'
    conda "environments/GENESIS_TOOLS.yaml"
    container "${dotenv('GENESIS_TOOLS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.metadata.tsv'
    input:
    tuple val(meta), path(reads), path(chrom_sizes)
    output:
    tuple val(meta), path("${meta.id}.metadata.tsv"), emit: metadata
    script:
    def files = reads instanceof List ? reads.sort { it.name } : [reads]
    def mate = meta.layout == 'PE' ? "--read2 '${files[1]}'" : ''
    """
    genesis-tools metadata --sample-id '${meta.id}' \
        --species '${meta.species}' --reference-fasta '${meta.reference_fasta}' \
        --chrom-sizes '${chrom_sizes}' \
        --read1 '${files[0]}' ${mate} --output '${meta.id}.metadata.tsv'
    """
    stub:
    def files = reads instanceof List ? reads.sort { it.name } : [reads]
    def mate = meta.layout == 'PE' ? "--read2 '${files[1]}'" : ''
    """
    genesis-tools metadata --sample-id '${meta.id}' \
        --species '${meta.species}' --reference-fasta '${meta.reference_fasta}' \
        --chrom-sizes '${chrom_sizes}' \
        --read1 '${files[0]}' ${mate} --output '${meta.id}.metadata.tsv'
    """
}
