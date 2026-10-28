include { dotenv } from 'plugin/nf-dotenv'

process RNA_TRACKS {
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
    tuple val(meta), path("${meta.id}.cpm.bw"), emit: cpm
    script:
    // Uniquely aligned reads (STAR gives them MAPQ 255), per base and in CPM; spliced reads cover
    // their exons only. Both strands together.
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"
    bamCoverage -b '${bam}' -o '${meta.id}.cpm.bw' --binSize 1 --minMappingQuality 255 \\
        --normalizeUsing CPM -p ${task.cpus}
    """
    stub:
    """
    touch '${meta.id}.cpm.bw'
    """
}
