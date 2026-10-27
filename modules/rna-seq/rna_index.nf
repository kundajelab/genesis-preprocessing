include { dotenv } from 'plugin/nf-dotenv'

process RNA_INDEX {
    cpus 8
    // STAR needs about 10 bytes per reference base to build its suffix array (32 GB for a 3 Gb
    // genome); maize-size genomes need ~30 GB. ref.length is the summed chrom sizes.
    memory { (4.GB + 14.B * ref.length) * task.attempt }
    time { 1.h * (1 + ref.length / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_STAR.yaml"
    container "${dotenv('GENESIS_STAR_IMAGE')}"
    tag "${ref.id}"
    input:
    tuple val(ref), path(fasta), path(gtf)
    output:
    tuple val(ref), path("${ref.id}.star"), emit: index
    script:
    // Index sizes from the STAR manual, for small genomes and assemblies with many scaffolds.
    def log2 = { double value -> Math.log(value) / Math.log(2) }
    def saIndexBases = Math.min(14, (int) Math.floor(log2(ref.length) / 2 - 1))
    def chrBinBits = Math.min(18, (int) Math.floor(log2(Math.max(ref.length / ref.contigs as double, 100d))))
    """
    gzip -dc '${fasta}' > genome.fa
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir '${ref.id}.star' \\
        --genomeFastaFiles genome.fa --sjdbGTFfile '${gtf}' --sjdbOverhang 100 \\
        --genomeSAindexNbases ${saIndexBases} --genomeChrBinNbits ${chrBinBits} \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()} --outTmpDir star.tmp
    rm genome.fa
    """
    stub:
    """
    mkdir '${ref.id}.star'
    touch '${ref.id}.star/SA'
    """
}
