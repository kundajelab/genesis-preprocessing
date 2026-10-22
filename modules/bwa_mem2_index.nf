include { dotenv } from 'plugin/nf-dotenv'

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
    // The FASTA hash is published with the index; GET_REFERENCES reuses products only while it matches.
    publishDir params.references, mode: 'copy', pattern: '{*.bwa-mem2,*.sha256}'
    input:
    tuple val(ref), path(fasta)
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
    # Written last, so a reference is recorded as built only once its index exists.
    printf '%s  %s\\n' '${ref.sha256}' '${ref.filename}' > '${ref.filename}.sha256'
    """
    stub:
    """
    mkdir '${ref.id}.bwa-mem2'
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        printf 'stub\\n' > '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
    printf '%s  %s\\n' '${ref.sha256}' '${ref.filename}' > '${ref.filename}.sha256'
    """
}
