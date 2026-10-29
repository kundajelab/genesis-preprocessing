include { dotenv } from 'plugin/nf-dotenv'

// Shell that writes <index>.sha256 in the format GET_STAR_INDEX checks: the FASTA and generated GTF
// hashes, STAR's version, the index parameters, then every index file's size in bytes. Run after
// the index exists, so the recorded sizes are its final sizes.
def recordIndex(ref) {
    // Spliced into the script as is, so these escapes are the shell's.
    def parameters = ref.star.collect { name, value -> "printf '%s  %s\\n' '${value}' '${name}'" }
    return """\
{
    printf '%s  %s\\n' '${ref.sha256}' '${ref.filename}'
    printf '%s  %s\\n' '${ref.gtf_sha256}' '${ref.id}.genes.gtf'
    printf '%s  STAR\\n' "\$(STAR --version)"
    ${parameters.join('\n    ')}
    for file in '${ref.star_index}'/*; do
        if [ -f "\$file" ]; then
            printf '%s  %s\\n' "\$(wc -c < "\$file" | tr -d ' ')" "\$file"
        fi
    done
} > '${ref.star_index}.sha256'"""
}

process RNA_INDEX {
    cpus 8
    // STAR builds its suffix array in about 14 bytes per reference base (about 36 GB for a 2.3 Gb
    // maize genome). ref.length is the summed chrom sizes.
    memory { (4.GB + 14.B * ref.length) * task.attempt }
    time { 1.h * (1 + ref.length / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_STAR.yaml"
    container "${dotenv('GENESIS_STAR_IMAGE')}"
    tag "${ref.id}"
    // GET_STAR_INDEX reuses a published index only while its record matches; stub runs never
    // publish into the shared references folder.
    // Runs resume by default (nextflow.config), and publishDir then keeps existing files;
    // overwrite, so a task taken from the cache still replaces the stale products that sent it
    // here.
    publishDir params.references, mode: 'copy', pattern: '{*.star,*.star.sha256}',
        overwrite: true, enabled: !workflow.stubRun
    input:
    // ref.star holds the index parameters and ref.star_index the folder name (GET_STAR_INDEX).
    tuple val(ref), path(fasta), path(gtf)
    output:
    tuple val(ref), path("${ref.star_index}"), emit: index
    path "${ref.star_index}.sha256", emit: record
    script:
    """
    gzip -dc '${fasta}' > genome.fa
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir '${ref.star_index}' \\
        --genomeFastaFiles genome.fa --sjdbGTFfile '${gtf}' --sjdbOverhang ${ref.star.sjdbOverhang} \\
        --genomeSAindexNbases ${ref.star.genomeSAindexNbases} \\
        --genomeChrBinNbits ${ref.star.genomeChrBinNbits} \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()} --outTmpDir star.tmp
    rm genome.fa
${recordIndex(ref)}
    """
    stub:
    """
    mkdir '${ref.star_index}'
    printf 'stub\\n' > '${ref.star_index}/SA'
    STAR() { echo STUB; }
${recordIndex(ref)}
    """
}
