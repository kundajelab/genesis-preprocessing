include { dotenv } from 'plugin/nf-dotenv'

process TSS_CLUSTERS {
    cpus 1
    memory { 4.GB * task.attempt }
    time { 1.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_TSS.yaml"
    container "${dotenv('GENESIS_TSS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.tss_clusters.bed'
    input:
    tuple val(meta), path(ctss)
    output:
    tuple val(meta), path("${meta.id}.tss_clusters.bed"), emit: clusters
    script:
    // paraclu clusters the tags of each strand; paraclu-cut keeps, as FANTOM does, clusters at
    // most --tss_cluster_max_length long whose density rises at least 2-fold, and drops clusters
    // inside larger ones. Clusters need --tss_cluster_min_tpm tags per million (at least 2 tags).
    """
    gzip -dc '${ctss}' > ctss.tsv
    total=\$(awk '{n += \$4} END {print n}' ctss.tsv)
    min=\$(awk -v t="\$total" -v tpm=${params.tss_cluster_min_tpm} \\
        'BEGIN {m = tpm * t / 1e6; m = (m == int(m)) ? m : int(m) + 1; print (m < 2 ? 2 : m)}')
    paraclu "\$min" ctss.tsv > paraclu.txt
    paraclu-cut -l ${params.tss_cluster_max_length} -d 2 paraclu.txt > clusters.txt
    # BED6+4: score is the cluster's tags, capped at BED's 1000; then its TPM, its dominant TSS
    # (0-based) and that TSS's tags, and the cluster's exact tags.
    awk -v OFS='\\t' '!/^#/ {print \$1, \$3 - 1, \$4, \$2, \$6}' clusters.txt \\
        | LC_ALL=C sort -k1,1 -k2,2n > clusters.bed
    awk -v OFS='\\t' '{print \$1, \$3 - 1, \$3, \$2, \$4}' ctss.tsv | LC_ALL=C sort -k1,1 -k2,2n > ctss.bed
    bedtools intersect -sorted -wa -wb -a clusters.bed -b ctss.bed \\
        | awk -v OFS='\\t' '\$4 == \$9 {
            key = \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5
            if (!(key in best) || \$10 > best[key]) {best[key] = \$10; at[key] = \$7}
        }
        END {for (key in best) print key, at[key], best[key]}' \\
        | LC_ALL=C sort -k1,1 -k2,2n -k6,6 \\
        | awk -v OFS='\\t' -v id='${meta.id}' -v total="\$total" '{
            printf "%s\\t%d\\t%d\\t%s_%d\\t%d\\t%s\\t%.4f\\t%d\\t%d\\t%d\\n",
                \$1, \$2, \$3, id, NR, (\$5 > 1000 ? 1000 : \$5), \$4, \$5 * 1e6 / total, \$6, \$7, \$5
        }' > '${meta.id}.tss_clusters.bed'
    """
    stub:
    """
    printf 'chr1\\t90\\t110\\t${meta.id}_1\\t5\\t+\\t1.0\\t99\\t5\\t5\\n' > '${meta.id}.tss_clusters.bed'
    """
}
