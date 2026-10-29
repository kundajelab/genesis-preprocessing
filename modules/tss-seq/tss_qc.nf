include { dotenv } from 'plugin/nf-dotenv'

process TSS_QC {
    cpus 1
    memory { 4.GB * task.attempt }
    time { 1.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_TSS.yaml"
    container "${dotenv('GENESIS_TSS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*_mqc.tsv'
    input:
    tuple val(meta), path(ctss), path(clusters), path(gtf), path(rrna_genes), path(trim_log),
        path(star_log), path(chrom_sizes)
    output:
    tuple val(meta), path("${meta.id}.tss_qc_mqc.tsv"), emit: summary
    script:
    def organelles = meta.organelles.join(',')
    // Tags near annotated transcript starts (the GTF's transcripts, strand-aware), on the
    // opposite strand nearby, on organelle contigs, in rRNA genes and in clusters. These are
    // diagnostics: no plant thresholds are applied.
    """
    # A misspelled organelle contig would silently report no organelle tags.
    for contig in \$(printf '%s' '${organelles}' | tr ',' ' '); do
        awk -v c="\$contig" '\$1 == c {found = 1} END {exit !found}' '${chrom_sizes}' || {
            printf 'Organelle contig %s is not in %s\\n' "\$contig" '${meta.reference_fasta}' >&2
            exit 1
        }
    done
    gzip -dc '${ctss}' | awk -v OFS='\\t' '{print \$1, \$3 - 1, \$3, ".", \$4, \$2}' \\
        | LC_ALL=C sort -k1,1 -k2,2n > ctss.bed
    awk -F '\\t' -v OFS='\\t' '\$3 == "transcript" {
            start = (\$7 == "+") ? \$4 : \$5
            print \$1, start - 1, start, ".", 0, \$7
        }' '${gtf}' | LC_ALL=C sort -u -k1,1 -k2,2n -k6,6 > starts.bed
    sum() { awk '{n += \$5} END {print n + 0}'; }
    tags=\$(sum < ctss.bed)
    near100=\$(bedtools window -a ctss.bed -b starts.bed -w 100 -sm -u | sum)
    near500=\$(bedtools window -a ctss.bed -b starts.bed -w 500 -sm -u | sum)
    antisense500=\$(bedtools window -a ctss.bed -b starts.bed -w 500 -Sm -u | sum)
    organelle=\$(awk -v list='${organelles}' 'BEGIN {n = split(list, o, ","); for (i = 1; i <= n; i++) org[o[i]] = 1}
        \$1 in org {t += \$5} END {print t + 0}' ctss.bed)
    # rRNA genes (RNA_ANNOTATION's list) over their exons, strand-aware; NA without any.
    if [ -s '${rrna_genes}' ]; then
        awk -F '\\t' -v OFS='\\t' 'NR == FNR {rrna[\$1] = 1; next}
            \$3 == "exon" && match(\$9, /gene_id "[^"]+"/) && substr(\$9, RSTART + 9, RLENGTH - 10) in rrna {
                print \$1, \$4 - 1, \$5, ".", 0, \$7
            }' '${rrna_genes}' '${gtf}' | LC_ALL=C sort -k1,1 -k2,2n > rrna.bed
        rrna=\$(bedtools intersect -a ctss.bed -b rrna.bed -s -u | sum)
    else
        rrna=NA
    fi
    # The cluster's exact tags (column 10); the BED score is capped at 1000.
    clustered=\$(awk '{n += \$10} END {print n + 0}' '${clusters}')
    cluster_count=\$(grep -c . '${clusters}' || true)
    input=\$(awk -F ':' '/^Total reads processed/ {gsub(/[^0-9]/, "", \$2); print \$2}' '${trim_log}')
    trimmed=\$(awk -F ':' '/^Reads written/ {sub(/\\(.*/, "", \$2); gsub(/[^0-9]/, "", \$2); print \$2}' '${trim_log}')
    unique=\$(awk -F '|' 'index(\$1, "Uniquely mapped reads number") {gsub(/[ \\t]/, "", \$2); print \$2}' '${star_log}')
    {
        printf '# id: tss_qc\\n# section_name: TSS libraries\\n# plot_type: table\\n'
        printf 'sample\\tmethod\\tinput_reads\\ttrimmed_reads\\tunique_reads\\ttags\\tnear_TSS_100\\tnear_TSS_500\\tantisense_500\\torganelle\\trRNA\\tclusters\\tin_clusters\\n'
        awk -v OFS='\\t' -v id='${meta.id}' -v method='${meta.method}' -v input="\$input" \\
            -v trimmed="\$trimmed" -v unique="\$unique" -v tags="\$tags" -v n100="\$near100" \\
            -v n500="\$near500" -v anti="\$antisense500" -v org="\$organelle" -v rrna="\$rrna" \\
            -v clusters="\$cluster_count" -v clustered="\$clustered" 'BEGIN {
                f = "%.4f"
                print id, method, input, trimmed, unique, tags, sprintf(f, n100 / tags),
                    sprintf(f, n500 / tags), sprintf(f, anti / tags), sprintf(f, org / tags),
                    (rrna == "NA" ? "NA" : sprintf(f, rrna / tags)), clusters, sprintf(f, clustered / tags)
            }'
    } > '${meta.id}.tss_qc_mqc.tsv'
    """
    stub:
    """
    printf 'sample\\n${meta.id}\\n' > '${meta.id}.tss_qc_mqc.tsv'
    """
}
