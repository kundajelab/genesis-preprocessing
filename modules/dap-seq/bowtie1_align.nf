include { dotenv } from 'plugin/nf-dotenv'

// Stock Bowtie 1 comparison module; no NH weighting is inferred here.
process BOWTIE1_ALIGN {
    cpus 2
    memory '8 GB'
    time '4h'
    conda "${moduleDir}/../../environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    input:
    tuple val(meta), path(reads), path(metadata), path(index)
    output:
    tuple val(meta), path("${meta.id}.bowtie1.bam"), path("${meta.id}.bowtie1.bam.bai"),
        path("${meta.id}.bowtie1.qc.bam"), path("${meta.id}.bowtie1.qc.bam.bai"), emit: aligned
    script:
    def files = reads instanceof List ? reads.sort { it.name } : [reads]
    def inputs = meta.layout == 'PE' ?
        "-1 '${files[0]}' -2 '${files[1]}' -X 1000" :
        "'${files[0]}'"
    """
    bowtie -p ${Math.max(1, task.cpus - 1)} -v 1 -a --best --strata --sam \
        "${index}/genome" ${inputs} | \
        samtools view -u -F 4 - | samtools sort -n -O SAM -o grouped.sam -
    samtools sort -m 512M -o '${meta.id}.bowtie1.bam' grouped.sam
    samtools index '${meta.id}.bowtie1.bam'
    bowtie -p ${Math.max(1, task.cpus - 1)} -v 2 -k 2 -m 1 --best --strata --sam \
        "${index}/genome" '${files[0]}' | \
        samtools view -u -F 4 - | samtools sort -m 512M -o '${meta.id}.bowtie1.qc.bam' -
    samtools index '${meta.id}.bowtie1.qc.bam'
    """
    stub:
    """
    touch '${meta.id}.bowtie1.bam' '${meta.id}.bowtie1.bam.bai'
    touch '${meta.id}.bowtie1.qc.bam' '${meta.id}.bowtie1.qc.bam.bai'
    """
}
