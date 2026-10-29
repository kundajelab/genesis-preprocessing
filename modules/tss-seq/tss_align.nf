include { dotenv } from 'plugin/nf-dotenv'

process TSS_ALIGN {
    cpus 8
    // As RNA_ALIGN: the loaded index is ~10 bytes per reference base, then samtools sorts.
    memory { (4.GB + 11.B * meta.genome_size) * task.attempt }
    time { 1.h * (1 + reads.size() / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_STAR.yaml"
    container "${dotenv('GENESIS_STAR_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.star.Log.final.out'
    input:
    tuple val(meta), path(reads), path(index)
    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.star.Log.final.out"), emit: log
    script:
    def sortMemory = Math.max(256L, (task.memory.toMega() * 3 / 4 / task.cpus) as long)
    // Extend5pOfRead1: the 5' end is never soft-clipped, so a mismatched cap G stays visible
    // (MD tag) for TSS_CTSS to correct. Template-switching methods align locally instead, which
    // soft-clips the added Gs that the genome lacks. Introns no longer than the longest annotated.
    def ends = meta.soft_clip_5p ? 'Local' : 'Extend5pOfRead1'
    """
    STAR --runThreadN ${task.cpus} --genomeDir '${index}' --readFilesCommand gunzip -c \\
        --readFilesIn '${reads}' --alignEndsType ${ends} \\
        --alignIntronMax ${meta.max_intron} --outFilterMultimapNmax 20 \\
        --outSAMattributes NH HI AS nM MD --outSAMtype BAM Unsorted \\
        --outFileNamePrefix star. --outTmpDir star.tmp
    samtools sort -@ ${task.cpus} -m ${sortMemory}M -T sort.tmp -o '${meta.id}.bam' \\
        star.Aligned.out.bam
    rm star.Aligned.out.bam
    samtools index '${meta.id}.bam'
    mv star.Log.final.out '${meta.id}.star.Log.final.out'
    """
    stub:
    """
    touch '${meta.id}.bam' '${meta.id}.bam.bai'
    printf 'stub\\n' > '${meta.id}.star.Log.final.out'
    """
}
