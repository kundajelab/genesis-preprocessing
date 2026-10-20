include { dotenv } from 'plugin/nf-dotenv'

process BWA_MEM2_ALIGN {
    cpus 8
    // The loaded index is ~5.5 bytes per reference base; samtools sort uses 512 MB, and each bwa-mem2
    // thread adds a little more.
    memory { (1.5.GB + 7.B * meta.genome_size + 256.MB * task.cpus) * task.attempt }
    // Each sample is aligned three times (all reads, read 1 for QC and tracks, and short read 1 for
    // SPP). Measured throughput gives a wide margin at one hour per GB of compressed reads.
    time { 1.h * (1 + [reads].flatten().sum { read -> read.size() } / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.alignment.tsv'
    input:
    tuple val(meta), path(reads), path(metadata), path(index)
    val bwa_options
    val spp_read_length
    output:
    tuple val(meta), path("${meta.id}.primary.bam"), path("${meta.id}.primary.bam.bai"),
        path("${meta.id}.qc.bam"), path("${meta.id}.qc.bam.bai"), emit: aligned
    tuple val(meta), path("${meta.id}.spp.bam"), emit: spp
    tuple val(meta), path("${meta.id}.alignment.tsv"), emit: provenance
    script:
    def files = reads instanceof List ? reads.sort { it.name } : [reads]
    def inputs = files.collect { "'${it}'" }.join(' ')
    def threads = Math.max(1, task.cpus - 1)
    def options = "-K ${bwa_options.batch_size} -k ${bwa_options.seed_length} " +
        "-c ${bwa_options.max_seed_occurrences} -T ${bwa_options.score_threshold}"
    """
    bwa-mem2 mem -t ${threads} ${options} -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA' \
        "${index}/genome" ${inputs} \
        2> main.bwa.log | samtools view -u -F 2308 - | \
        samtools sort -m 512M -o '${meta.id}.primary.bam' -
    samtools index '${meta.id}.primary.bam'
    bwa-mem2 mem -t ${threads} ${options} -R '@RG\\tID:${meta.id}.qc\\tSM:${meta.id}\\tPL:ILLUMINA' \
        "${index}/genome" '${files[0]}' \
        2> qc.bwa.log | samtools view -u -F 2308 - | \
        samtools sort -m 512M -o '${meta.id}.qc.bam' -
    samtools index '${meta.id}.qc.bam'
    # SPP gets read 1 cut to its first bases before alignment: with full-length reads as long as the
    # fragments, the read-length phantom peak hides the fragment peak and SPP cannot estimate it.
    bwa-mem2 mem -t ${threads} ${options} -R '@RG\\tID:${meta.id}.spp\\tSM:${meta.id}\\tPL:ILLUMINA' \
        "${index}/genome" <(gzip -dc '${files[0]}' | \
            awk -v bases=${spp_read_length} 'NR % 2 == 0 {\$0 = substr(\$0, 1, bases)} 1') \
        2> spp.bwa.log | samtools view -u -F 2308 - | \
        samtools sort -m 512M -o '${meta.id}.spp.bam' -
    printf 'sample_id\\taligner\\tversion\\tlayout\\tread_length\\tmapq_filter\\tflag_exclude\\tbwa_options\\tqc_reads\\tspp_reads\\n' > '${meta.id}.alignment.tsv'
    printf '%s\\tbwa-mem2\\t%s\\t%s\\t%s\\tnone\\t2308\\t%s\\tfull_read1_unpaired\\tread1_first_${spp_read_length}bp_unpaired\\n' \
        '${meta.id}' "\$(bwa-mem2 version)" '${meta.layout}' '${meta.read_length}' '${options}' >> '${meta.id}.alignment.tsv'
    """
    stub:
    """
    touch '${meta.id}.primary.bam' '${meta.id}.primary.bam.bai' '${meta.id}.qc.bam' '${meta.id}.qc.bam.bai' \\
        '${meta.id}.spp.bam'
    printf 'stub\\n' > '${meta.id}.alignment.tsv'
    """
}
