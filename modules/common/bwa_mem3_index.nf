include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomErrorStrategy } from './retry'

// The index files under <id>.bwa-mem2/genome that alignment needs. bwa-mem3 writes no .0123 (the
// unpacked reference bwa-mem2 also loaded); it reads bases from .pac instead.
def indexSuffixes() {
    return ['.amb', '.ann', '.bwt.2bit.64', '.pac']
}

// Shell that writes <fasta>.sha256: the FASTA hash, then each product's size in bytes.
// Run after the products exist, so the recorded sizes are their final sizes.
def recordProducts(ref, chromSizes) {
    return """\
{
    printf '%s  %s\\n' '${ref.sha256}' '${ref.filename}'
    printf '%s  %s\\n' "\$(wc -c < '${chromSizes}' | tr -d ' ')" '${ref.id}.chrom.sizes'
    for suffix in ${indexSuffixes().join(' ')}; do
        printf '%s  %s\\n' "\$(wc -c < '${ref.id}.bwa-mem2/genome'"\$suffix" | tr -d ' ')" \\
            '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
} > '${ref.filename}.sha256'"""
}

process BWA_MEM3_INDEX {
    // bwa-mem3 index builds the Sorghum reference (718 Mb) in 77 s on 1 thread, 46 s on 4, and 42 s
    // on 16: past 4, writing the FM-index takes most of the time. It peaks at ~14.5 bytes per
    // reference base whatever the threads; ref.length is the summed chrom sizes.
    cpus 4
    memory { oomMemoryOf(task, 1.GB + 16.B * ref.length) }
    time { 1.h * (1 + ref.length / 1e9) * task.attempt }
    // Retry with twice the memory only when an attempt ran out of memory (modules/common/retry.nf).
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${ref.id}"
    // <fasta>.sha256 records the FASTA hash and every product's size; GET_REFERENCES reuses products
    // only while both match. Stub runs never publish into the shared references folder.
    // Runs resume by default (nextflow.config), and publishDir then keeps existing files;
    // overwrite, so a task taken from the cache still replaces the stale products that sent it
    // here.
    publishDir params.references, mode: 'copy', pattern: '{*.bwa-mem2,*.sha256}',
        overwrite: true, enabled: !workflow.stubRun
    input:
    tuple val(ref), path(fasta), path(chrom_sizes)
    output:
    tuple val(ref), path("${ref.id}.bwa-mem2"), emit: index
    path "${ref.filename}.sha256", emit: sha256
    script:
    """
    mkdir '${ref.id}.bwa-mem2'
    # bwa-mem3 reads gzipped FASTA directly. Its thread count and memory budget default to what the
    # machine has, so give it the task's. The index is the same for any number of threads.
    bwa-mem3 index -t ${task.cpus} --max-memory ${task.memory.toMega()}M --tmp-dir . \\
        -p '${ref.id}.bwa-mem2/genome' '${fasta}'
    for suffix in ${indexSuffixes().join(' ')}; do
        test -s '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
${recordProducts(ref, chrom_sizes)}
    """
    stub:
    """
    mkdir '${ref.id}.bwa-mem2'
    for suffix in ${indexSuffixes().join(' ')}; do
        printf 'stub\\n' > '${ref.id}.bwa-mem2/genome'"\$suffix"
    done
${recordProducts(ref, chrom_sizes)}
    """
}
