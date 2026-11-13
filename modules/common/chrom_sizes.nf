include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; retryTimeOf ; retryErrorStrategy } from './retry'

process CHROM_SIZES {
    cpus 1
    memory { oomMemoryOf(task, 1.GB) }
    time { retryTimeOf(task, 1.h) }
    // Retry with twice the memory or time when an attempt ran out of it (modules/common/retry.nf).
    errorStrategy { retryErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${ref.id}"
    // Stub runs never publish into the shared references folder: their fake products would be reused.
    // Runs resume by default (nextflow.config), and publishDir then keeps existing files;
    // overwrite, so a task taken from the cache still replaces the stale products that sent it
    // here.
    publishDir params.references, mode: 'copy', pattern: '*.chrom.sizes', overwrite: true,
        enabled: !workflow.stubRun
    input:
    tuple val(ref), path(fasta)
    output:
    tuple val(ref), path("${ref.id}.chrom.sizes"), emit: sizes
    script:
    """
    gzip -dc -- '${fasta}' > reference.fa
    samtools faidx reference.fa
    cut -f1,2 reference.fa.fai > '${ref.id}.chrom.sizes'
    """
    stub:
    """
    printf 'chr1\\t1000\\n' > '${ref.id}.chrom.sizes'
    """
}
