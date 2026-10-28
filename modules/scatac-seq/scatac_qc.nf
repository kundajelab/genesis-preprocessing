include { dotenv } from 'plugin/nf-dotenv'

process SCATAC_QC {
    cpus { resources.cpus }
    memory { resources.memory }
    time { resources.time }
    conda "${moduleDir}/../../environments/GENESIS_DEEPTOOLS.yaml"
    container "${dotenv('GENESIS_DEEPTOOLS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.library}/groups/${meta.group}",
        mode: 'copy', pattern: '{*_mqc.tsv,*.tss.*}'
    input:
    // genes is an empty list when the sheet gives no gene annotation.
    tuple val(meta), path(fragments), path(cuts_bam), path(cuts_bai), path(peaks), path(cuts), path(genes)
    val resources
    output:
    tuple val(meta), path("${meta.id}.scatac_qc_mqc.tsv"), emit: summary
    tuple val(meta), path("${meta.id}.tss.*"), optional: true, emit: tss
    script:
    // TSS enrichment exactly as ATAC_QC computes it, on this pseudobulk's cut sites.
    def tss = genes ? """
    computeMatrix reference-point -S '${cuts}' -R '${genes}' --referencePoint TSS -b 2000 -a 2000 \\
        --binSize 1 --missingDataAsZero -p ${task.cpus} -o tss.matrix.gz
    plotProfile -m tss.matrix.gz -o '${meta.id}.tss.png' --outFileNameData '${meta.id}.tss.profile.tab'
    tss=\$(awk 'END {
        for (i = 3; i < 103; i++) flank += \$i
        for (i = NF - 99; i <= NF; i++) flank += \$i
        print (flank ? \$2003 / (flank / 200) : "NA")
    }' '${meta.id}.tss.profile.tab')
    """ : 'tss=NA'
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"
    gzip -dc '${peaks}' > peaks.bed
    peaks=\$(wc -l < peaks.bed | tr -d ' ')
    cells=\$(cut -f4 '${fragments}' | sort -u | wc -l | tr -d ' ')
    fragments=\$(wc -l < '${fragments}' | tr -d ' ')
    insertions=\$(samtools view -c '${cuts_bam}')
    in_peaks=\$(samtools view -c -L peaks.bed '${cuts_bam}')
    ${tss}
    {
        printf '# id: scatac_qc\\n# section_name: scATAC pseudobulks\\n# plot_type: table\\n'
        printf 'pseudobulk\\tcells\\tfragments\\tinsertions\\tinsertions_in_peaks\\tFRiP\\tpeaks\\tTSS_enrichment\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' '${meta.id}' "\$cells" "\$fragments" \\
            "\$insertions" "\$in_peaks" \\
            "\$(awk -v n="\$in_peaks" -v d="\$insertions" 'BEGIN {print (d ? n / d : "NA")}')" \\
            "\$peaks" "\$tss"
    } > '${meta.id}.scatac_qc_mqc.tsv'
    """
    stub:
    """
    printf 'pseudobulk\\n${meta.id}\\n' > '${meta.id}.scatac_qc_mqc.tsv'
    """
}
