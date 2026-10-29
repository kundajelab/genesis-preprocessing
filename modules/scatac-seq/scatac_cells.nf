include { dotenv } from 'plugin/nf-dotenv'

process SCATAC_CELLS {
    cpus 2
    memory { 4.GB * task.attempt }
    time { 2.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '*.barcodes.tsv'
    input:
    tuple val(meta), path(fragments)
    output:
    tuple val(meta), path(fragments), path("${meta.id}.source.cells.tsv"), emit: cells
    tuple val(meta), path("${meta.id}.barcodes.tsv"), emit: barcodes
    script:
    def organelles = meta.organelles.join(',')
    """
    # Nuclear and organelle fragments per barcode, most nuclear fragments first.
    bgzip --decompress --stdout --threads ${task.cpus} '${fragments}' | awk -F '\\t' -v OFS='\\t' \\
        -v organelles='${organelles}' '
        BEGIN {n = split(organelles, names, ","); for (i = 1; i <= n; i++) organelle[names[i]] = 1}
        /^#/ {next}
        {if (\$1 in organelle) org[\$4]++; else nuclear[\$4]++; seen[\$4] = 1}
        END {for (b in seen) print b, nuclear[b] + 0, org[b] + 0}' |
        LC_ALL=C sort -t \$'\\t' -k2,2nr -k1,1 > counts.tsv
    {
        printf 'barcode\\tnuclear_fragments\\torganelle_fragments\\tcell\\n'
        awk -F '\\t' -v OFS='\\t' -v min=${params.scatac_min_fragments} -v max=${params.scatac_max_cells} '
            {cell = \$2 >= min && (max == 0 || NR <= max); print \$1, \$2, \$3, cell ? "yes" : "no"}' counts.tsv
    } > '${meta.id}.barcodes.tsv'
    if [ '${meta.cells_url}' = - ]; then
        # Called cells form one group, all; no TSS threshold is applied, as none is tuned for plants.
        {
            printf 'barcode\\tgroup\\n'
            awk -F '\\t' 'NR > 1 && \$4 == "yes" {print \$1 "\\tall"}' '${meta.id}.barcodes.tsv'
        } > '${meta.id}.source.cells.tsv'
    else
        aria2c --max-connection-per-server=4 --retry-wait=5 --out=authors.tsv '${meta.cells_url}'
        # The authors' barcodes may carry Cell Ranger's -1 suffix; chromap writes none.
        awk -F '\\t' -v OFS='\\t' 'NR == 1 {print; next} {sub(/-1\$/, "", \$1); print}' authors.tsv \\
            > '${meta.id}.source.cells.tsv'
    fi
    cells=\$(( \$(wc -l < '${meta.id}.source.cells.tsv') - 1 ))
    if [ "\$cells" -lt 1 ]; then
        printf 'No cells: no barcode has %s nuclear fragments\\n' '${params.scatac_min_fragments}' >&2
        exit 1
    fi
    printf '%s cells\\n' "\$cells" >&2
    """
    stub:
    """
    printf 'barcode\\tgroup\\nACGTACGTACGTACGT\\tall\\n' > '${meta.id}.source.cells.tsv'
    printf 'barcode\\tnuclear_fragments\\torganelle_fragments\\tcell\\n' > '${meta.id}.barcodes.tsv'
    """
}
