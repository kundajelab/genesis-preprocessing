include { dotenv } from 'plugin/nf-dotenv'

process SCATAC_TRACKS {
    cpus { resources.cpus }
    memory { resources.memory }
    time { resources.time }
    conda "${moduleDir}/../../environments/DAP_SEQ_TRACKS.yaml"
    container "${dotenv('DAP_SEQ_TRACKS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.library}/groups/${meta.group}",
        mode: 'copy', pattern: '*.bw'
    input:
    tuple val(meta), path(fragments), path(chrom_sizes)
    val resources
    output:
    tuple val(meta), path("${meta.id}.cuts.counts.bw"), emit: cuts
    tuple val(meta), path("${meta.id}.cuts.bam"), path("${meta.id}.cuts.bam.bai"), emit: bam
    script:
    // Tn5 cut sites: the first and last base of each fragment, whose ends already carry the Tn5
    // offsets. These are the bases bulk ATAC counts after its +4/-5 shift. Each cut is a 1-base
    // alignment, so deepTools and samtools count cuts exactly as they count bulk ATAC reads.
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"
    {
        awk -v OFS='\\t' '{print "@SQ", "SN:" \$1, "LN:" \$2}' '${chrom_sizes}'
        awk -v OFS='\\t' '{
            print "c" NR "s", 0, \$1, \$2 + 1, 255, "1M", "*", 0, 0, "*", "*"
            print "c" NR "e", 0, \$1, \$3, 255, "1M", "*", 0, 0, "*", "*"
        }' '${fragments}'
    } | samtools sort -m 512M -o '${meta.id}.cuts.bam' -
    samtools index '${meta.id}.cuts.bam'
    bamCoverage -b '${meta.id}.cuts.bam' -o '${meta.id}.cuts.counts.bw' --binSize 1 \\
        --normalizeUsing None -p ${task.cpus}
    """
    stub:
    """
    touch '${meta.id}.cuts.counts.bw' '${meta.id}.cuts.bam' '${meta.id}.cuts.bam.bai'
    """
}
