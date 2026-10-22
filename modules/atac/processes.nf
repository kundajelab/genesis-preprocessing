process ATAC_ACQUIRE {
    tag "${meta.library_id}"
    container params.images.tools
    publishDir "${params.outdir}/${meta.library_id}/provenance", mode: 'copy', pattern: 'acquisition.json'
    cache 'deep'
    maxRetries 2
    errorStrategy { task.exitStatus == 75 && task.attempt <= 2 ? 'retry' : 'finish' }
    input:
    tuple val(meta), path(manifest), path(localReads, stageAs: 'inputs??/*')
    path code
    val sourceHash
    output:
    tuple val(meta), path('raw_R1.fastq.gz'), path('raw_R2.fastq.gz'), emit: reads
    tuple val(meta), path('acquisition.json'), emit: provenance
    script:
    """
    PYTHONPATH='${code}' python -m genesis_tools.atac.stages acquire --manifest '${manifest}'
    """
}

process ATAC_REFERENCE {
    tag "${ref.reference_id}"
    container params.images.alignment
    cpus 1
    memory { 1.GB + 32.B * ref.total_bases }
    cache 'deep'
    input:
    tuple val(ref), path(fasta, stageAs: 'input.fa.gz')
    output:
    tuple val(ref.reference_id), path('index'), emit: index
    script:
    """
    mkdir index
    gzip -dc input.fa.gz > index/genome.fa
    samtools faidx index/genome.fa
    cut -f1,2 index/genome.fa.fai > index/chrom.sizes
    bwa-mem2 index index/genome.fa
    bwa-mem2 version > index/bwa.version.txt
    (cd index; sha256sum genome.fa* chrom.sizes bwa.version.txt > checksums.sha256)
    """
}

process ATAC_FASTQC {
    tag "${meta.library_id}:${mate}"
    container params.images.fastqc
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}/${meta.library_id}/qc/fastqc", mode: 'copy'
    input:
    tuple val(meta), val(mate), path(read)
    output:
    tuple val(meta), path('*_fastqc.zip'), path('*_fastqc.html'), emit: reports
    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    fastqc --threads 1 --memory 512 --noextract '${read}'
    test -s raw_${mate}_fastqc.zip
    test -s raw_${mate}_fastqc.html
    """
}

process ATAC_ADAPTERS {
    tag "${meta.library_id}"
    container params.images.cutadapt
    publishDir "${params.outdir}/${meta.library_id}/qc", mode: 'copy', pattern: 'adapters.json'
    input:
    tuple val(meta), path(read1), path(read2)
    output:
    tuple val(meta), path('analysis_R1.fastq.gz'), path('analysis_R2.fastq.gz'), emit: reads
    path 'adapters.json'
    script:
    def command = meta.adapter_r1 == '-' ?
        "ln -s '${read1}' analysis_R1.fastq.gz; ln -s '${read2}' analysis_R2.fastq.gz; printf '%s\\n' '{\"policy\":\"none\"}' > adapters.json" :
        "cutadapt -j ${task.cpus} -a '${meta.adapter_r1}' -A '${meta.adapter_r2}' --json adapters.json -o analysis_R1.fastq.gz -p analysis_R2.fastq.gz '${read1}' '${read2}'"
    """
    ${command}
    """
}

process ATAC_ALIGN {
    tag "${meta.library_id}"
    container params.images.alignment
    publishDir "${params.outdir}/${meta.library_id}/qc", mode: 'copy', pattern: '{raw.*.txt,duplicate-marking.txt,alignment.versions.txt}'
    memory { 4.GB + 8.B * meta.reference.total_bases }
    input:
    tuple val(meta), path(read1), path(read2), path(index)
    output:
    tuple val(meta), path('marked.bam'), emit: bam
    tuple val(meta), path('raw.*.txt'), emit: qc
    path 'duplicate-marking.txt'
    path 'alignment.versions.txt'
    script:
    """
    (cd '${index}'; sha256sum -c checksums.sha256)
    bwa-mem2 mem -t ${task.cpus} -R '@RG\\tID:${meta.library_id}\\tSM:${meta.sample_id}\\tLB:${meta.library_id}' '${index}/genome.fa' '${read1}' '${read2}' |
      samtools sort -n -m 128M -o query.bam
    samtools fixmate -m query.bam fixed.bam
    samtools sort -m 128M -o coordinate.bam fixed.bam
    samtools markdup -s coordinate.bam marked.bam 2> duplicate-marking.txt
    samtools index marked.bam
    samtools flagstat marked.bam > raw.flagstat.txt
    samtools stats marked.bam > raw.stats.txt
    samtools idxstats marked.bam > raw.idxstats.txt
    { bwa-mem2 version; samtools --version; } > alignment.versions.txt
    """
}

process ATAC_FRAGMENTS {
    tag "${meta.library_id}"
    container params.images.tools
    publishDir "${params.outdir}/${meta.library_id}/fragments", mode: 'copy', pattern: '*.{gz,tbi,json,bedgraph,txt,sizes}'
    input:
    tuple val(meta), path(bam), path(manifest)
    path code
    val sourceHash
    output:
    tuple val(meta), path('usable.bam'), emit: bam
    tuple val(meta), path('fragments.bed.gz'), path('fragments.bed.gz.tbi'), path('chrom.sizes'), emit: fragments
    tuple val(meta), path('usable.*.txt'), emit: qc
    tuple val(meta), path('metrics.json'), emit: metrics
    tuple val(meta), path('cuts.bedgraph'), path('chrom.sizes'), emit: cuts
    script:
    """
    PYTHONPATH='${code}' python -m genesis_tools.atac.stages fragments --manifest '${manifest}'
    """
}

process ATAC_PEAKS {
    tag "${meta.library_id}"
    container params.images.peaks
    publishDir "${params.outdir}/${meta.library_id}/peaks", mode: 'copy', pattern: '{atac*,macs3.version.txt}'
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path('atac_peaks.narrowPeak'), emit: peaks
    path 'macs3.version.txt'
    script:
    """
    macs3 callpeak -t '${bam}' -f BAMPE -g ${meta.reference.genome_size} -q 0.01 --keep-dup all -n atac \\
        --threads ${task.cpus}
    macs3 --version > macs3.version.txt
    """
}

process ATAC_ENRICHMENT {
    tag "${meta.library_id}"
    container params.images.tools
    publishDir "${params.outdir}/${meta.library_id}/qc", mode: 'copy', pattern: 'enrichment.json'
    input:
    tuple val(meta), path(bam), path(fragments), path(sizes), path(peaks), path(tss, stageAs: 'tss.bed')
    path code
    val sourceHash
    output:
    tuple val(meta), path('enrichment.json'), emit: metrics
    script:
    """
    PYTHONPATH='${code}' python -m genesis_tools.atac.stages enrichment
    """
}

process ATAC_REPORT {
    tag "${meta.library_id}"
    container params.images.multiqc
    publishDir "${params.outdir}/${meta.library_id}/qc", mode: 'copy'
    input:
    tuple val(meta), path(reports)
    path code
    val sourceHash
    output:
    tuple val(meta), path('multiqc_report.html'), path('multiqc_data'), emit: report
    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    printf '%s\\n' 'data_dir_name: multiqc_data' > multiqc_config.yaml
    multiqc . --force --require-logs -m fastqc -m samtools --filename multiqc_report.html --data-dir --data-format json --config multiqc_config.yaml
    PYTHONPATH='${code}' python -m genesis_tools.atac.stages report
    """
}

process ATAC_TRACKS {
    tag "${meta.library_id}"
    container params.images.tracks
    publishDir "${params.outdir}/${meta.library_id}/tracks", mode: 'copy'
    input:
    tuple val(meta), path(cuts), path(sizes)
    output:
    tuple val(meta), path('cuts.counts.bw'), emit: track
    script:
    """
    python - <<'PY'
    import pyBigWig
    sizes = [line.split() for line in open('${sizes}')]
    with pyBigWig.open('cuts.counts.bw', 'w') as target:
        target.addHeader([(chrom, int(n)) for chrom, n in sorted(sizes)])
        batch = []
        def flush():
            if batch:
                target.addEntries([r[0] for r in batch], [int(r[1]) for r in batch],
                                  ends=[int(r[2]) for r in batch], values=[float(r[3]) for r in batch])
                batch.clear()
        for line in open('${cuts}'):
            batch.append(line.split())
            if len(batch) >= 16384:
                flush()
        flush()
    PY
    """
}
