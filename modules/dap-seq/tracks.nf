include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomErrorStrategy } from '../common/retry'

process TRACKS {
    // bamCoverage runs 1.72x as fast on 4 CPUs as on 2 (on a 1 GB BAM), and 2.7x on 8 for 48%
    // more CPU time; 4 is the better trade. On 2 CPUs, 57 Sorghum samples used 183-196% CPU.
    cpus 4
    // On Perlmutter, 142 Sorghum samples' processes peaked at 0.85 GB on 2 CPUs and 1.58 GB on 4,
    // whatever the BAM size (up to 28.8 GB). Yet 17 large ones had bamCoverage workers killed for
    // memory at 5.5 GB, apparently for bamCoverage's temporary files, which a compute node's /tmp
    // keeps in memory; the script puts them in the task's folder instead. bamCoverage waits for
    // killed workers until the time limit, which modules/common/retry.nf retries as out of memory.
    memory { oomMemoryOf(task, 2.0.GB) }
    // On 2 CPUs, 1 hour for the Sorghum genome (718 Mb, which the 1-bp tracks span) plus 0.24
    // hours per GB of the BAMs tracked covered every sample; the largest, 17.9 GB, took 3.64
    // hours. Allow twice that at 4 CPUs' speed (proportionally longer on fewer), plus 10 minutes.
    time {
        def tracked = bam.size() + (meta.layout == 'PE' ? qc_bam.size() : 0)
        def hours = (1.4 * meta.genome_size + 0.24 * tracked) / 1e9 * 2 / 1.72 * 4 / task.cpus
        10.min + 1.h * hours
    }
    // Retry with twice the memory only when an attempt ran out of memory (modules/common/retry.nf).
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_DEEPTOOLS.yaml"
    container "${dotenv('GENESIS_DEEPTOOLS_IMAGE')}"
    tag "${meta.id}"
    publishDir { "${params.publish_base}/output/${meta.species}/${meta.id}" }, mode: 'copy',
        pattern: '*.bw'
    input:
    tuple val(meta), path(bam), path(bai), path(qc_bam), path(qc_bai)
    output:
    tuple val(meta), path("${meta.id}.tracks*.bw"), emit: tracks
    tuple val(meta), path("${meta.id}.tracks.RPKM.bedGraph"), emit: rpkm
    script:
    def qcTracks = meta.layout == 'PE' ?
        "make_tracks '${qc_bam}' '${meta.id}.tracks.read1'" : ''
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    # bamCoverage's temporary files go in the task's folder, not a /tmp that may be held in memory.
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$MPLCONFIGDIR" "\$TMPDIR"
    make_tracks() {
        local input=\$1 stem=\$2 count factor
        count=\$(samtools view -c "\$input")
        test "\$count" -gt 0
        factor=\$(awk -v count="\$count" 'BEGIN {printf "%.17g", 1000000/count}')
        bamCoverage -b "\$input" -o "\$stem.CPM.bw" --binSize 1 \\
            --normalizeUsing CPM --exactScaling -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.plus.CPM.bw" --binSize 1 \\
            --normalizeUsing None --scaleFactor "\$factor" --samFlagExclude 16 -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.minus.CPM.bw" --binSize 1 \\
            --normalizeUsing None --scaleFactor "\$factor" --samFlagInclude 16 -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.plus.5p.counts.bw" --binSize 1 \\
            --normalizeUsing None --Offset 1 --samFlagExclude 16 -p ${task.cpus}
        bamCoverage -b "\$input" -o "\$stem.minus.5p.counts.bw" --binSize 1 \\
            --normalizeUsing None --Offset 1 --samFlagInclude 16 -p ${task.cpus}
    }
    make_tracks '${bam}' '${meta.id}.tracks'
    ${qcTracks}
    bamCoverage -b '${bam}' -o '${meta.id}.tracks.RPKM.bedGraph' \\
        --outFileFormat bedgraph --binSize 1 --normalizeUsing RPKM --exactScaling -p ${task.cpus}
    """
    stub:
    """
    touch '${meta.id}.tracks.CPM.bw' '${meta.id}.tracks.plus.CPM.bw' '${meta.id}.tracks.minus.CPM.bw'
    touch '${meta.id}.tracks.plus.5p.counts.bw' '${meta.id}.tracks.minus.5p.counts.bw'
    printf 'chr1\\t0\\t1000\\t0\\n' > '${meta.id}.tracks.RPKM.bedGraph'
    """
}
