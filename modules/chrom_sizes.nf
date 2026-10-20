include { dotenv } from 'plugin/nf-dotenv'

process CHROM_SIZES {
    cpus 1
    memory '1 GB'
    time '1h'
    conda "environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${ref.id}"
    publishDir params.references, mode: 'copy', pattern: '*.chrom.sizes'
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
