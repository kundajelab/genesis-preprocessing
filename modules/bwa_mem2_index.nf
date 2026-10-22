include { dotenv } from 'plugin/nf-dotenv'

// Shell that writes <fasta>.sha256: the FASTA hash, then each product's size in bytes.
// Run after the products exist, so the recorded sizes are their final sizes.
def recordProducts(ref, chromSizes) {
    return """\
{
    printf '%s  %s\\n' '${ref.sha256}' '${ref.filename}'
    printf '%s  %s\\n' "\$(wc -c < '${chromSizes}' | tr -d ' ')" '${ref.id}.chrom.sizes'
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        printf '%s  %s\\n' "\$(wc -c < '${ref.id}.bwa-mem2/genome'"\$suffix" | tr -d ' ')" \\
            '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
} > '${ref.filename}.sha256'"""
}

process BWA_MEM2_INDEX {
    // bwa-mem2 index is single-threaded. It peaks at ~28 bytes per reference base (measured 28.3,
    // and documented as 28N), so request 32 plus headroom. ref.length is the summed chrom sizes.
    cpus 1
    memory { (1.GB + 32.B * ref.length) * task.attempt }
    time { 1.h * (1 + ref.length / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${ref.id}"
    // <fasta>.sha256 records the FASTA hash and every product's size; GET_REFERENCES reuses products
    // only while both match. Stub runs never publish into the shared references folder.
    publishDir params.references, mode: 'copy', pattern: '{*.bwa-mem2,*.sha256}',
        enabled: !workflow.stubRun
    input:
    tuple val(ref), path(fasta), path(chrom_sizes)
    output:
    tuple val(ref), path("${ref.id}.bwa-mem2"), emit: index
    path "${ref.filename}.sha256", emit: sha256
    script:
    """
    mkdir '${ref.id}.bwa-mem2'
    # bwa-mem2 reads gzipped FASTA directly, producing an index identical to the uncompressed one.
    bwa-mem2 index -p '${ref.id}.bwa-mem2/genome' '${fasta}'
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        test -s '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
${recordProducts(ref, chrom_sizes)}
    """
    stub:
    """
    mkdir '${ref.id}.bwa-mem2'
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        printf 'stub\\n' > '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
${recordProducts(ref, chrom_sizes)}
    """
}
