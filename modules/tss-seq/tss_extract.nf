include { dotenv } from 'plugin/nf-dotenv'

process TSS_EXTRACT {
    cpus 2
    memory '2 GB'
    time { 1.h * (1 + reads.sum { read -> read.size() } / 2e9) * task.attempt }
    // retry with more time when killed (walltime or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_TSS.yaml"
    container "${dotenv('GENESIS_TSS_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.tso.log'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.id}.tss_reads.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.tso.log"), emit: log
    script:
    // Smar2C2: a read from the cap carries the Smart-seq2 primer (AAGCAGTGGTATCAACGCAGAGTAC), the
    // UMI and AT, then the template-switch Gs and the RNA. Either mate can carry it, and the primer
    // can sit anywhere in the read: Tn5 also cuts inside the oligo upstream of it (in GSE197144
    // most primers start at base 6, many at 50-110). The rest are internal Tn5 fragments. Reads
    // with Ts after the primer instead come from the oligo-dT primer (3' ends) and are dropped.
    // The UMI goes to the end of the read name for deduplication; the Gs stay for alignment to
    // soft-clip.
    def umi = meta.umi_length
    """
    for read in ${reads.collect { read -> "'${read}'" }.join(' ')}; do gzip -dc "\$read"; done | gawk -v umi=${umi} '
        BEGIN {
            pattern = "AAGCAGTGGTATCAACGCAGAGTAC[ACGTN]{" umi "}AT"
        }
        NR % 4 == 1 {name = \$1}
        NR % 4 == 2 {sequence = \$0}
        NR % 4 == 0 {
            reads++
            if (match(sequence, pattern)) {
                start = RSTART + RLENGTH
                code = substr(sequence, RSTART + RLENGTH - umi - 2, umi)
                kept++
                print name "_" code "\\n" substr(sequence, start) "\\n+\\n" substr(\$0, start)
            } else if (sequence ~ /AAGCAGTGGTATCAACGCAGAGTACT{10}/) {
                oligo_dt++
            }
        }
        END {
            printf "reads\\t%d\\ntss_reads\\t%d\\noligo_dt_reads\\t%d\\n", reads, kept, oligo_dt > "counts.tsv"
        }' | gzip -c > '${meta.id}.tss_reads.fastq.gz'
    mv counts.tsv '${meta.id}.tso.log'
    test "\$(awk '\$1 == "tss_reads" {print \$2}' '${meta.id}.tso.log')" -gt 0 || {
        printf 'No read has the Smar2C2 primer, UMI and AT\\n' >&2
        exit 1
    }
    """
    stub:
    """
    printf '@r_AAAAAAAAAAAA\\nGGGACGT\\n+\\nIIIIIII\\n' | gzip > '${meta.id}.tss_reads.fastq.gz'
    printf 'reads\\t2\\ntss_reads\\t1\\noligo_dt_reads\\t0\\n' > '${meta.id}.tso.log'
    """
}
