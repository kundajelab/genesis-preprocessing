include { dotenv } from 'plugin/nf-dotenv'

process ATAC_QC {
    cpus { resources.cpus }
    memory { resources.memory }
    time { resources.time }
    conda "${moduleDir}/../../environments/DAP_SEQ_TRACKS.yaml"
    container "${dotenv('DAP_SEQ_TRACKS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '{*_mqc.tsv,*.tss.*}'
    input:
    // genes is an empty list when the sheet gives no gene annotation.
    tuple val(meta), path(bam), path(bai), path(peaks), path(cuts), path(complexity), path(genes)
    val resources
    output:
    tuple val(meta), path("${meta.id}.atac_qc_mqc.tsv"), emit: summary
    tuple val(meta), path("${meta.id}.tss.*"), optional: true, emit: tss
    script:
    // TSS enrichment as ENCODE reports it: mean cut count at the TSS divided by the mean over the
    // outer 100 bases at each end of the +/-2 kb window. Gene starts come from the strand in BED6.
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
    reads=\$(samtools view -c '${bam}')
    in_peaks=\$(samtools view -c -L peaks.bed '${bam}')
    ${tss}
    {
        printf '# id: atac_qc\\n# section_name: ATAC QC\\n# plot_type: table\\n'
        printf 'sample\\tusable_reads\\treads_in_peaks\\tFRiP\\tpeaks\\tTSS_enrichment\\t'
        head -n 1 '${complexity}'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t' '${meta.id}' "\$reads" "\$in_peaks" \\
            "\$(awk -v n="\$in_peaks" -v d="\$reads" 'BEGIN {print (d ? n / d : "NA")}')" "\$peaks" "\$tss"
        tail -n 1 '${complexity}'
    } > '${meta.id}.atac_qc_mqc.tsv'
    """
    stub:
    """
    printf 'sample\\n${meta.id}\\n' > '${meta.id}.atac_qc_mqc.tsv'
    """
}
