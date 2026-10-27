include { dotenv } from 'plugin/nf-dotenv'

process RNA_QC {
    cpus 4
    memory '2 GB'
    time { 1.h * (1 + bam.size() / 1e9) }
    conda "${moduleDir}/../../environments/GENESIS_RNA_COUNTS.yaml"
    container "${dotenv('GENESIS_RNA_COUNTS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*_mqc.tsv'
    input:
    tuple val(meta), path(bam), path(bai), path(gtf), path(rrna_genes), path(counts), path(summary),
        path(star_log)
    output:
    tuple val(meta), path("${meta.id}.rna_qc_mqc.tsv"), emit: summary
    script:
    def organelles = meta.organelles.join(' ')
    def pairs = meta.layout == 'PE' ? '-p --countReadPairs' : ''
    """
    # A misspelled organelle contig would silently report no organelle reads.
    for contig in ${organelles}; do
        samtools view -H '${bam}' | awk -F '\\t' '\$1 == "@SQ" {print \$2}' | grep -qx "SN:\$contig" || {
            printf 'Organelle contig %s is not in %s\\n' "\$contig" '${meta.reference_fasta}' >&2
            exit 1
        }
    done
    star() { awk -F '|' -v key="\$1" 'index(\$1, key) {gsub(/[\\t %]/, "", \$2); print \$2}' '${star_log}'; }
    # Libraries are often submitted with the wrong strandedness: count the reads that would be
    # assigned with each orientation. Near 1 is forward, near 0 reverse, near 0.5 unstranded.
    for strand in 1 2; do
        featureCounts -T ${task.cpus} -a '${gtf}' -t exon -g gene_id -s "\$strand" ${pairs} \\
            -o "strand\$strand.tsv" '${bam}' 2> /dev/null
    done
    forward=\$(awk '\$1 == "Assigned" {a = \$2} END {print a + 0}' strand1.tsv.summary)
    reverse=\$(awk '\$1 == "Assigned" {a = \$2} END {print a + 0}' strand2.tsv.summary)
    # Unique alignments on organelle contigs, of all unique alignments.
    unique=\$(samtools view -c -q 255 -F 0x900 '${bam}')
    organelle=\$(if [ -n '${organelles}' ]; then samtools view -c -q 255 -F 0x900 '${bam}' ${organelles}; else echo 0; fi)
    # Assigned reads in rRNA genes, of all assigned reads; NA when the annotation names none
    # (MaizeGDB and Phytozome gene sets often have no rRNA genes).
    rrna=NA
    if [ -s '${rrna_genes}' ]; then
        rrna=\$(awk -F '\\t' 'NR == FNR {r[\$1] = 1; next} FNR > 2 {all += \$7; if (\$1 in r) n += \$7}
            END {print (all ? n / all : "NA")}' '${rrna_genes}' '${counts}')
    fi
    assigned=\$(awk -F '\\t' 'NR > 1 {all += \$2; if (\$1 == "Assigned") a = \$2}
        END {print (all ? a / all : "NA")}' '${summary}')
    {
        printf '# id: rna_qc\\n# section_name: RNA-seq libraries\\n# plot_type: table\\n'
        printf 'sample\\tstrandedness\\tforward_fraction\\tinput_reads\\tunique\\tmultimapped\\tassigned\\trRNA\\torganelle\\tmax_intron\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' '${meta.id}' '${meta.strandedness}' \\
            "\$(awk -v f="\$forward" -v r="\$reverse" 'BEGIN {print (f + r ? f / (f + r) : "NA")}')" \\
            "\$(star 'Number of input reads')" \\
            "\$(star 'Uniquely mapped reads %' | awk '{print \$1 / 100}')" \\
            "\$(star '% of reads mapped to multiple loci' | awk '{print \$1 / 100}')" \\
            "\$assigned" "\$rrna" \\
            "\$(awk -v n="\$organelle" -v d="\$unique" 'BEGIN {print (d ? n / d : "NA")}')" \\
            '${meta.max_intron}'
    } > '${meta.id}.rna_qc_mqc.tsv'
    """
    stub:
    """
    printf 'sample\\n${meta.id}\\n' > '${meta.id}.rna_qc_mqc.tsv'
    """
}
