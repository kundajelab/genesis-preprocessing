include { dotenv } from 'plugin/nf-dotenv'

process QC {
    cpus { resources.cpus }
    memory { resources.memory }
    time { resources.time }
    conda "environments/DAP_SEQ_QC.yaml"
    container "${dotenv('DAP_SEQ_QC_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.qc.report.*'
    input:
    tuple val(meta), path(bam), path(bai), path(qc_bam), path(qc_bai), path(spp_bam)
    val resources
    path spp_script
    output:
    tuple val(meta), path("${meta.id}.qc.report.*"), emit: reports
    script:
    def insert = meta.layout == 'PE' ? """
        awk 'BEGIN {OFS="\\t"; print "insert_length","count","fraction"}
            \$1 == "IS" {n++; lengths[n]=\$2; counts[n]=\$3; total+=\$3}
            END {for (i=1;i<=n;i++) print lengths[i],counts[i],(total ? counts[i]/total : 0)}' \
            '${meta.id}.qc.report.main.stats.txt' > '${meta.id}.qc.report.insert_lengths.tsv'
        """ : ''
    """
    samtools flagstat '${bam}' > '${meta.id}.qc.report.flagstat.txt'
    samtools idxstats '${bam}' > '${meta.id}.qc.report.idxstats.tsv'
    samtools stats '${bam}' > '${meta.id}.qc.report.main.stats.txt'
    samtools stats '${qc_bam}' > '${meta.id}.qc.report.read1.stats.txt'
    ${insert}
    mkdir qc-bin
    ln -s "\$(command -v gawk)" qc-bin/awk
    export PATH="\$PWD/qc-bin:\$PATH"
    # run_spp.R crashes when no fragment-length peak stands out (e.g. little enrichment). Record an
    # NA result from what it logged rather than failing the run; any other failure is still fatal.
    if ! Rscript "${spp_script}" -c='${spp_bam}' \
        -p=${task.cpus} -rf -s=-0:2:400 \
        '-savp=${meta.id}.qc.report.spp.pdf' '-out=${meta.id}.qc.report.spp.tsv' > spp.log 2>&1; then
        cat spp.log >&2
        grep -q '^Top 3 estimates for fragment length NA' spp.log || exit 1
        awk -v OFS='\t' -v name='${spp_bam}' '
            /^done[.] read [0-9]+ fragments/ {reads = \$3}
            /^Phantom peak location/ {phantom = \$NF}
            /^Phantom peak Correlation/ {phantom_cc = \$NF}
            /^Minimum cross-correlation shift/ {min_shift = \$NF}
            /^Minimum cross-correlation value/ {min_cc = \$NF}
            END {print name, reads, "NA", "NA", phantom, phantom_cc, min_shift, min_cc, "NA", "NA", "NA"}
        ' spp.log > '${meta.id}.qc.report.spp.tsv'
    fi
    cat spp.log
    """
    stub:
    """
    printf 'stub\\n' > '${meta.id}.qc.report.spp.tsv'
    touch '${meta.id}.qc.report.spp.pdf' '${meta.id}.qc.report.flagstat.txt'
    """
}
