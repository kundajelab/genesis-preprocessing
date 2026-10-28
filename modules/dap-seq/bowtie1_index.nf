include { dotenv } from 'plugin/nf-dotenv'

// Available for future comparison workflows; not called by workflows/dap-seq.nf.
process BOWTIE1_INDEX {
    cpus { resources.cpus }
    memory { resources.memory }
    time { resources.time }
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${ref.id}"
    publishDir params.references, mode: 'copy', pattern: '*.bowtie1'
    input:
    tuple val(ref), path(fasta)
    val resources
    output:
    tuple val(ref), path("${ref.id}.bowtie1"), emit: index
    script:
    """
    gzip -dc -- '${fasta}' > reference.fa
    mkdir '${ref.id}.bowtie1'
    bowtie-build --threads ${task.cpus} reference.fa '${ref.id}.bowtie1/genome'
    """
    stub:
    """
    mkdir '${ref.id}.bowtie1'
    for suffix in 1 2 3 4 rev.1 rev.2; do
        printf 'stub\\n' > '${ref.id}.bowtie1/genome.'"\$suffix"'.ebwt'
    done
    """
}
