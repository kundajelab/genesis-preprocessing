include { dotenv } from 'plugin/nf-dotenv'

process ATAC_TRACKS {
    cpus { resources.cpus }
    memory { resources.memory }
    time { resources.time }
    conda "${moduleDir}/../../environments/DAP_SEQ_TRACKS.yaml"
    container "${dotenv('DAP_SEQ_TRACKS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.bw'
    input:
    tuple val(meta), path(bam), path(bai)
    val resources
    output:
    tuple val(meta), path("${meta.id}.cuts.counts.bw"), emit: cuts
    script:
    // Tn5 cut sites: each read's 5' base after the usual +4/-5 shift (alignmentSieve --ATACshift),
    // i.e. BED0 start + 4 on the forward strand and end - 6 on the reverse. Raw counts per base.
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"
    alignmentSieve -b '${bam}' -o shifted.unsorted.bam --ATACshift -p ${task.cpus}
    samtools sort -m 512M -o shifted.bam shifted.unsorted.bam
    samtools index shifted.bam
    bamCoverage -b shifted.bam -o '${meta.id}.cuts.counts.bw' --binSize 1 --Offset 1 \\
        --normalizeUsing None -p ${task.cpus}
    """
    stub:
    """
    touch '${meta.id}.cuts.counts.bw'
    """
}
