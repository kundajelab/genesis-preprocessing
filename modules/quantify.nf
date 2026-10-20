include { dotenv } from 'plugin/nf-dotenv'

process QUANTIFY {
    tag "${meta.id}"
    // One streaming pass over the BAM; the second CPU decompresses BGZF blocks.
    cpus 2
    memory '2 GB'
    time { 1.h * (1 + bam.size() / 1e9) }
    conda "environments/GENESIS_TOOLS.yaml"
    container "${dotenv('GENESIS_TOOLS_IMAGE')}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.peaks.*.tsv'
    input:
    tuple val(meta), path(peaks), path(bam), path(coverage), path(chrom_sizes)
    output:
    tuple val(meta), path("${meta.id}.peaks.RPM.tsv"),
        path("${meta.id}.peaks.mean_RPKM.tsv"), emit: quantified
    script:
    """
    genesis-tools quantify --bed '${peaks}' --bam '${bam}' --coverage '${coverage}' \\
        --chrom-sizes '${chrom_sizes}' --weighting primary --threads ${task.cpus} \\
        --rpm-output '${meta.id}.peaks.RPM.tsv' --rpkm-output '${meta.id}.peaks.mean_RPKM.tsv'
    """
    stub:
    """
    gzip -dc '${peaks}' | awk -v OFS='\\t' '{print \$0, 0}' > '${meta.id}.peaks.RPM.tsv'
    cp '${meta.id}.peaks.RPM.tsv' '${meta.id}.peaks.mean_RPKM.tsv'
    """
}
