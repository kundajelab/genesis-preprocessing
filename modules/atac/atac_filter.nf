include { dotenv } from 'plugin/nf-dotenv'

process ATAC_FILTER {
    cpus 2
    memory { 2.GB * task.attempt }
    time { 1.h * (1 + bam.size() / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.fragments.bed.gz*'
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '{*.usable.flagstat.txt,*.complexity.tsv}'
    input:
    tuple val(meta), path(bam), path(bai), path(chrom_sizes)
    val mapq
    output:
    tuple val(meta), path("${meta.id}.usable.bam"), path("${meta.id}.usable.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.fragments.bed.gz"), path("${meta.id}.fragments.bed.gz.tbi"),
        emit: fragments
    tuple val(meta), path("${meta.id}.complexity.tsv"), emit: complexity
    tuple val(meta), path("${meta.id}.usable.flagstat.txt"), emit: reports
    script:
    def organelles = meta.organelles.join(',')
    // Each proper pair's outer span, BED0 half-open, from read 1: the leftmost mate's start plus the
    // template length.
    def spans = "awk -v OFS='\\t' '{start = (\$4 < \$8 ? \$4 : \$8) - 1; " +
        "print \$3, start, start + (\$9 < 0 ? -\$9 : \$9)}' | LC_ALL=C sort -S 512M -T . -k1,1 -k2,2n -k3,3n"
    """
    # A misspelled organelle contig would silently count organelle reads as nuclear.
    for contig in \$(printf '%s' '${organelles}' | tr ',' ' '); do
        awk -v contig="\$contig" '\$1 == contig {found = 1} END {exit !found}' '${chrom_sizes}' || {
            printf 'Organelle contig %s is not in %s\\n' "\$contig" '${meta.reference_fasta}' >&2
            exit 1
        }
    done
    awk -v OFS='\\t' -v organelles='${organelles}' '
        BEGIN {n = split(organelles, names, ","); for (i = 1; i <= n; i++) skip[names[i]] = 1}
        !(\$1 in skip) {print \$1, 0, \$2}' '${chrom_sizes}' > nuclear.bed
    # Keep pairs where both mates are primary, properly paired, not QC-failed, on a nuclear contig and
    # at least MAPQ ${mapq}. A read whose mate fails loses its mate here; fixmate then clears its
    # proper-pair flag, so the second filter drops it. Duplicates stay until complexity is counted.
    samtools view -u -f 2 -F 2828 -q ${mapq} -L nuclear.bed '${bam}' | \\
        samtools collate -O -u - | \\
        samtools fixmate -u - - | \\
        samtools view -u -f 2 -F 2828 - | \\
        samtools sort -m 512M -o filtered.bam -
    # Library complexity as ENCODE defines it, over fragments before duplicate removal.
    samtools view -f 66 filtered.bam | ${spans} | uniq -c | awk -v OFS='\\t' '
        {total += \$1; distinct++; one += (\$1 == 1); two += (\$1 == 2)}
        END {
            print "total_fragments", "distinct_fragments", "one_read_pair", "two_read_pairs", "NRF", "PBC1", "PBC2"
            print total, distinct, one, two, (total ? distinct / total : "NA"),
                (distinct ? one / distinct : "NA"), (two ? one / two : "NA")
        }' > '${meta.id}.complexity.tsv'
    samtools view -b -F 1024 -o '${meta.id}.usable.bam' filtered.bam
    samtools index '${meta.id}.usable.bam'
    samtools flagstat '${meta.id}.usable.bam' > '${meta.id}.usable.flagstat.txt'
    samtools view -f 66 '${meta.id}.usable.bam' | ${spans} | bgzip > '${meta.id}.fragments.bed.gz'
    tabix -p bed '${meta.id}.fragments.bed.gz'
    """
    stub:
    """
    touch '${meta.id}.usable.bam' '${meta.id}.usable.bam.bai' '${meta.id}.fragments.bed.gz.tbi'
    printf 'chr1\\t100\\t300\\n' | gzip > '${meta.id}.fragments.bed.gz'
    printf 'stub\\n' > '${meta.id}.complexity.tsv'
    printf 'stub\\n' > '${meta.id}.usable.flagstat.txt'
    """
}
