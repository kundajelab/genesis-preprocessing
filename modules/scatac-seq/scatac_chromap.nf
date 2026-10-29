include { dotenv } from 'plugin/nf-dotenv'

process SCATAC_CHROMAP {
    cpus 16
    // The index is small (about 2 bytes per base); the barcode whitelist and per-cell duplicate
    // removal need most of the memory, which grows with the library's fragments.
    memory { (16.GB + 4.B * meta.genome_length) * task.attempt }
    time { 4.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_CHROMAP.yaml"
    container "${dotenv('GENESIS_CHROMAP_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.chromap.log'
    // Every barcode's fragments, for calling cells again with other tools.
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.source.fragments.gz',
        saveAs: { name -> name.replace('.source.fragments.gz', '.all_barcodes.fragments.tsv.gz') }
    input:
    // One library: the R1, R2 and barcode FASTQs of all its runs, in run order.
    tuple val(meta), path(r1, stageAs: 'r1/*'), path(r2, stageAs: 'r2/*'),
        path(barcodes, stageAs: 'barcodes/*'), path(fasta), path(whitelist)
    output:
    tuple val(meta), path("${meta.id}.source.fragments.gz"), emit: fragments
    tuple val(meta), path("${meta.id}.chromap.log"), emit: log
    script:
    def join = { files -> (files instanceof List ? files : [files]).join(',') }
    def corrected = meta.barcode_whitelist != '-'
    // The cell barcode is the last 16 bases of the barcode read: all of a 10x ATAC barcode read,
    // and after the 8-base spacer in a multiome ATAC one.
    def start = meta.barcode_length - 16
    """
    gzip -dc '${fasta}' > genome.fa
    chromap -i -r genome.fa -o genome.index
    orientation=+
    if ${corrected}; then
        gzip -dcf '${whitelist}' > known-barcodes.txt
        # Sequencers read the barcode on either strand; use the one that matches the whitelist.
        # Only the first 200,000 barcodes: the reader stops early, so zcat's SIGPIPE is expected.
        { zcat ${join(barcodes).replace(',', ' ')} || [ \$? -eq 141 ]; } | head -n 800000 |
            awk -v start=${start} '
            NR == FNR {list[\$1] = 1; next}
            FNR % 4 == 2 {
                code = substr(\$0, start + 1, 16); n++
                forward += code in list
                reverse = ""; for (i = 16; i >= 1; i--) reverse = reverse substr("TGCAN", index("ACGTN", substr(code, i, 1)), 1)
                backward += reverse in list
            }
            END {printf "%d\\t%d\\t%d\\n", n, forward, backward > "orientation.tsv"}' known-barcodes.txt -
        read -r n forward backward < orientation.tsv
        if [ "\$backward" -gt "\$forward" ]; then orientation=-; fi
        best=\$(( forward > backward ? forward : backward ))
        printf 'Barcode orientation %s: %s of %s barcodes in the whitelist\\n' "\$orientation" "\$best" "\$n" >&2
        if [ \$(( best * 2 )) -lt "\$n" ]; then
            printf 'Under half of the barcodes are in %s: wrong whitelist or chemistry?\\n' '${meta.barcode_whitelist}' >&2
            exit 1
        fi
    fi
    # --preset atac: adapter trimming, per-cell duplicate removal, and 10x fragments (BED0, Tn5
    # offsets applied: start + 4, end - 5), as Cell Ranger writes them.
    chromap --preset atac -t ${task.cpus} -x genome.index -r genome.fa \\
        -1 '${join(r1)}' -2 '${join(r2)}' -b '${join(barcodes)}' \\
        --read-format "r1:0:-1,bc:${start}:${meta.barcode_length - 1}:\$orientation,r2:0:-1" \\
        ${corrected ? '--barcode-whitelist known-barcodes.txt' : ''} \\
        -o fragments.tsv 2> '${meta.id}.chromap.log'
    rm -f genome.fa genome.index known-barcodes.txt
    test -s fragments.tsv || { cat '${meta.id}.chromap.log' >&2; exit 1; }
    bgzip --threads ${task.cpus} -c fragments.tsv > '${meta.id}.source.fragments.gz'
    rm fragments.tsv
    """
    stub:
    """
    printf 'chr1\\t100\\t300\\tACGTACGTACGTACGT\\t1\\n' | gzip > '${meta.id}.source.fragments.gz'
    touch '${meta.id}.chromap.log'
    """
}
