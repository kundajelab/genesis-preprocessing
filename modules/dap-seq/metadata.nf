include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from '../common/retry'

process METADATA {
    cpus 1
    memory { oomMemoryOf(task, 1.GB) }
    time { retryTimeOf(task, 1.h) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_TOOLS.yaml"
    container "${dotenv('GENESIS_TOOLS_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
        pattern: '*.metadata.tsv'
    input:
    tuple val(meta), path(reads), path(chrom_sizes)
    output:
    tuple val(meta), path("${meta.id}.metadata.tsv"), emit: metadata
    script:
    def files = reads instanceof List ? reads.sort { read -> read.name } : [reads]
    def mate = meta.layout == 'PE' ? "--read2 '${files[1]}'" : ''
    """
    genesis-tools metadata --sample-id '${meta.id}' \
        --species '${meta.species}' --reference-fasta '${meta.reference_fasta}' \
        --chrom-sizes '${chrom_sizes}' \
        --read1 '${files[0]}' ${mate} --output '${meta.id}.metadata.tsv'
    """
    stub:
    def files = reads instanceof List ? reads.sort { read -> read.name } : [reads]
    def mate = meta.layout == 'PE' ? "--read2 '${files[1]}'" : ''
    """
    genesis-tools metadata --sample-id '${meta.id}' \
        --species '${meta.species}' --reference-fasta '${meta.reference_fasta}' \
        --chrom-sizes '${chrom_sizes}' \
        --read1 '${files[0]}' ${mate} --output '${meta.id}.metadata.tsv'
    """
}
