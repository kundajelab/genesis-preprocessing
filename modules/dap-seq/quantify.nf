include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomErrorStrategy } from '../common/retry'

process QUANTIFY {
    tag "${meta.id}"
    // One streaming pass over the BAM; the second CPU decompresses BGZF blocks. On Perlmutter, 19
    // Sorghum samples used 131-146% CPU.
    cpus 2
    // They peaked at 56-137 MB whatever the BAM size (1.3-11.1 GB).
    memory { oomMemoryOf(task, 512.MB) }
    // They took 52 seconds per GB of BAM (at most 58), 9.4 minutes for the largest. Allow twice
    // that, at least 10 minutes.
    time { 10.min + 2.min * bam.size() / 1e9 }
    // Retry with twice the memory only when an attempt ran out of memory (modules/common/retry.nf).
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_TOOLS.yaml"
    container "${dotenv('GENESIS_TOOLS_IMAGE')}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
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
