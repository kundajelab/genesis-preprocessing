include { dotenv } from 'plugin/nf-dotenv'

process TSS_CTSS {
    cpus 2
    memory { 4.GB * task.attempt }
    time { 1.h * (1 + bam.size() / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_TSS.yaml"
    container "${dotenv('GENESIS_TSS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '{*.ctss.tsv.gz,*.bw}'
    input:
    tuple val(meta), path(bam), path(bai), path(chrom_sizes)
    output:
    tuple val(meta), path("${meta.id}.ctss.tsv.gz"), emit: ctss
    tuple val(meta), path("${meta.id}.*.bw"), emit: tracks
    script:
    // One tag per uniquely aligned read (STAR gives them MAPQ 255) at its 5' end: the alignment
    // start on +, its last aligned base on -. Cap-based methods add a non-templated G before the
    // first base; when the read starts with a G that the genome lacks (MD tag), the TSS is the
    // next base, as in FANTOM.
    """
    samtools view -q 255 -F 0x904 '${bam}' | gawk -v OFS='\\t' -v correct=${meta.g_correction ? 1 : 0} '
        {
            md = ""
            for (i = 12; i <= NF; i++) if (\$i ~ /^MD:Z:/) md = substr(\$i, 6)
            if (and(\$2, 16)) {
                span = 0; cigar = \$6
                while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
                    n = substr(cigar, 1, RLENGTH - 1) + 0; op = substr(cigar, RLENGTH, 1)
                    if (op ~ /[MDN=X]/) span += n
                    cigar = substr(cigar, RLENGTH + 1)
                }
                strand = "-"; pos = \$4 + span - 1
                if (correct && substr(\$10, length(\$10), 1) == "C" && md ~ /[ACGTN]0\$/) pos--
            } else {
                strand = "+"; pos = \$4
                if (correct && substr(\$10, 1, 1) == "G" && md ~ /^0[ACGTN]/) pos++
            }
            count[\$3 "\\t" strand "\\t" pos]++
        }
        END {for (key in count) print key, count[key]}' |
        LC_ALL=C sort -t \$'\\t' -k1,1 -k2,2 -k3,3n > ctss.tsv
    test -s ctss.tsv || { printf 'No uniquely aligned reads\\n' >&2; exit 1; }
    gzip -n -c ctss.tsv > '${meta.id}.ctss.tsv.gz'
    # Per strand: raw tag counts and counts per million tags, as bigWigs (bedGraph is 0-based).
    total=\$(awk '{n += \$4} END {print n}' ctss.tsv)
    for strand in plus minus; do
        sign=\$( [ "\$strand" = plus ] && echo + || echo - )
        awk -F '\\t' -v OFS='\\t' -v s="\$sign" '\$2 == s {print \$1, \$3 - 1, \$3, \$4}' ctss.tsv \\
            | LC_ALL=C sort -k1,1 -k2,2n > "\$strand.bedGraph"
        awk -v OFS='\\t' -v total="\$total" '{print \$1, \$2, \$3, \$4 * 1e6 / total}' \\
            "\$strand.bedGraph" > "\$strand.cpm.bedGraph"
        bedGraphToBigWig "\$strand.bedGraph" '${chrom_sizes}' '${meta.id}.ctss.'"\$strand"'.bw'
        bedGraphToBigWig "\$strand.cpm.bedGraph" '${chrom_sizes}' '${meta.id}.ctss.'"\$strand"'.cpm.bw'
    done
    """
    stub:
    """
    printf 'chr1\\t+\\t100\\t5\\n' | gzip > '${meta.id}.ctss.tsv.gz'
    for f in plus minus plus.cpm minus.cpm; do touch '${meta.id}.ctss.'"\$f"'.bw'; done
    """
}
