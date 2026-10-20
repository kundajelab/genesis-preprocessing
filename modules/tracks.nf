include { dotenv } from 'plugin/nf-dotenv'

process TRACKS {
    cpus { resources.cpus }
    memory { resources.memory }
    time { resources.time }
    conda "environments/DAP_SEQ_TRACKS.yaml"
    container "${dotenv('DAP_SEQ_TRACKS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.bw'
    input:
    tuple val(meta), path(bam), path(bai), path(qc_bam), path(qc_bai)
    val resources
    output:
    tuple val(meta), path("${meta.id}.tracks*.bw"), emit: tracks
    tuple val(meta), path("${meta.id}.tracks.RPKM.bedGraph"), emit: rpkm
    script:
    def qcTracks = meta.layout == 'PE' ?
        "make_tracks '${qc_bam}' '${meta.id}.tracks.read1'" : ''
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"
    make_tracks() {
        local input=\$1 stem=\$2 count factor
        count=\$(samtools view -c "\$input")
        test "\$count" -gt 0
        factor=\$(awk -v count="\$count" 'BEGIN {printf "%.17g", 1000000/count}')
        bamCoverage -b "\$input" -o "\$stem.CPM.bw" --binSize 1 \
            --normalizeUsing CPM --exactScaling -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.plus.CPM.bw" --binSize 1 \
            --normalizeUsing None --scaleFactor "\$factor" --samFlagExclude 16 -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.minus.CPM.bw" --binSize 1 \
            --normalizeUsing None --scaleFactor "\$factor" --samFlagInclude 16 -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.plus.5p.counts.bw" --binSize 1 \
            --normalizeUsing None --Offset 1 --samFlagExclude 16 -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.minus.5p.counts.bw" --binSize 1 \
            --normalizeUsing None --Offset 1 --samFlagInclude 16 -p ${task.cpus}
    }
    make_tracks '${bam}' '${meta.id}.tracks'
    ${qcTracks}
    bamCoverage -b '${bam}' -o '${meta.id}.tracks.RPKM.bedGraph' \
        --outFileFormat bedgraph --binSize 1 --normalizeUsing RPKM --exactScaling -p ${task.cpus}
    """
    stub:
    """
    touch '${meta.id}.tracks.CPM.bw' '${meta.id}.tracks.plus.CPM.bw' '${meta.id}.tracks.minus.CPM.bw'
    touch '${meta.id}.tracks.plus.5p.counts.bw' '${meta.id}.tracks.minus.5p.counts.bw'
    printf 'chr1\\t0\\t1000\\t0\\n' > '${meta.id}.tracks.RPKM.bedGraph'
    """
}
