include { dotenv } from 'plugin/nf-dotenv'

// Shell that writes <ref_id>.chromap.sha256 in the format GET_CHROMAP_INDEX checks: the FASTA
// hash, chromap's release, then every index file's size in bytes. Run after the index exists, so
// the recorded sizes are its final sizes.
def recordIndex(ref) {
    // Spliced into the script as is, so these escapes are the shell's.
    return """\
{
    printf '%s  %s\\n' '${ref.sha256}' '${ref.filename}'
    # chromap --version prints the release and build, e.g. 0.3.2-r518; the pin names the release.
    printf '%s  chromap\\n' "\$(chromap --version 2>&1 | sed 's/-r[0-9]*\$//')"
    for file in '${ref.id}.chromap'/*; do
        if [ -f "\$file" ]; then
            printf '%s  %s\\n' "\$(wc -c < "\$file" | tr -d ' ')" "\$file"
        fi
    done
} > '${ref.id}.chromap.sha256'"""
}

process CHROMAP_INDEX {
    cpus 2
    // chromap peaks at about 10 bytes per reference base while indexing, for an index of about
    // 5 bytes per base (measured on 50 and 200 Mb genomes). ref.length is the summed chrom sizes.
    memory { (2.GB + 12.B * ref.length) * task.attempt }
    time { 1.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_CHROMAP.yaml"
    container "${dotenv('GENESIS_CHROMAP_IMAGE')}"
    tag "${ref.id}"
    // GET_CHROMAP_INDEX reuses a published index only while its record matches; stub runs never
    // publish into the shared references folder.
    // Runs resume by default (nextflow.config), and publishDir then keeps existing files;
    // overwrite, so a task taken from the cache still replaces the stale products that sent it
    // here.
    publishDir params.references, mode: 'copy', pattern: '{*.chromap,*.chromap.sha256}',
        overwrite: true, enabled: !workflow.stubRun
    input:
    // ref carries the FASTA's sha256 (GET_REFERENCES) and length (its chromosome sizes).
    tuple val(ref), path(fasta)
    output:
    tuple val(ref), path("${ref.id}.chromap"), emit: index
    path "${ref.id}.chromap.sha256", emit: record
    script:
    """
    gzip -dc '${fasta}' > genome.fa
    mkdir '${ref.id}.chromap'
    chromap -i -r genome.fa -o '${ref.id}.chromap/genome.index'
    rm genome.fa
${recordIndex(ref)}
    """
    stub:
    """
    mkdir '${ref.id}.chromap'
    printf 'stub\\n' > '${ref.id}.chromap/genome.index'
    chromap() { echo STUB; }
${recordIndex(ref)}
    """
}
