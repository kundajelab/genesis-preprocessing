include { dotenv } from 'plugin/nf-dotenv'

process RNA_ALIGN {
    cpus 8
    // The loaded index is ~10 bytes per reference base; sorting gets its own 2 GB.
    memory { (4.GB + 11.B * meta.genome_size) * task.attempt }
    time { 1.h * (1 + [reads].flatten().sum { read -> read.size() } / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_STAR.yaml"
    container "${dotenv('GENESIS_STAR_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.star.*'
    input:
    // reads is one FASTQ for single-end libraries and both mates for paired-end ones.
    tuple val(meta), path(reads), path(index)
    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.star.Log.final.out"), emit: log
    tuple val(meta), path("${meta.id}.star.SJ.out.tab"), emit: junctions
    script:
    def files = [reads].flatten().sort { it.name }
    // No intron longer than the longest annotated intron (see RNA_ANNOTATION). The gap between
    // mates can hold such an intron plus the unread middle of the insert, up to 1 kb.
    // Up to 20 loci per read, as ENCODE does; gene counts later keep only unique alignments.
    """
    STAR --runThreadN ${task.cpus} --genomeDir '${index}' --readFilesCommand zcat \\
        --readFilesIn ${files.collect { "'${it}'" }.join(' ')} \\
        --alignIntronMax ${meta.max_intron} --alignMatesGapMax ${meta.max_intron + 1000} \\
        --outFilterMultimapNmax 20 --outSAMattributes NH HI AS nM \\
        --outSAMattrRGline 'ID:${meta.id}' 'SM:${meta.id}' 'PL:ILLUMINA' \\
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 2000000000 \\
        --outFileNamePrefix star. --outTmpDir star.tmp
    mv star.Aligned.sortedByCoord.out.bam '${meta.id}.bam'
    samtools index '${meta.id}.bam'
    # MultiQC recognises the STAR log by its name.
    mv star.Log.final.out '${meta.id}.star.Log.final.out'
    mv star.SJ.out.tab '${meta.id}.star.SJ.out.tab'
    """
    stub:
    """
    touch '${meta.id}.bam' '${meta.id}.bam.bai' '${meta.id}.star.SJ.out.tab'
    printf 'stub\\n' > '${meta.id}.star.Log.final.out'
    """
}
