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
    conda "environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${ref.id}"
    publishDir params.references, mode: 'copy', pattern: '*.bwa-mem2'
    input:
    tuple val(ref), path(fasta)
    output:
    tuple val(ref), path("${ref.id}.bwa-mem2"), emit: index
    script:
    """
    mkdir '${ref.id}.bwa-mem2'
    # bwa-mem2 reads gzipped FASTA directly, producing an index identical to the uncompressed one.
    bwa-mem2 index -p '${ref.id}.bwa-mem2/genome' '${fasta}'
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        test -s '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
    """
    stub:
    """
    mkdir '${ref.id}.bwa-mem2'
    for suffix in .0123 .amb .ann .bwt.2bit.64 .pac; do
        printf 'stub\\n' > '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
    """
}
