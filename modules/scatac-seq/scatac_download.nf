include { dotenv } from 'plugin/nf-dotenv'

process SCATAC_DOWNLOAD {
    // One library at a time, as DOWNLOAD does, to keep the load on the remote servers bounded.
    maxForks 1
    cpus 2
    memory '1 GB'
    time '4h'
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    input:
    val meta
    output:
    tuple val(meta), path("${meta.id}.source.fragments.gz"), path("${meta.id}.source.cells.tsv"),
        emit: downloads
    script:
    """
    cat > aria2c-input.txt << 'END_OF_URLS'
${meta.fragments_url}
  out=${meta.id}.source.fragments.gz
${meta.cells_url}
  out=${meta.id}.source.cells.tsv
END_OF_URLS
    aria2c --input-file=aria2c-input.txt --max-connection-per-server=4 --split=4 \\
        --retry-wait=5 --rlimit-nofile \$(ulimit -n)
    # Decompressing the whole file checks its gzip integrity.
    bgzip --decompress --stdout --threads ${task.cpus} '${meta.id}.source.fragments.gz' > /dev/null
    """
    stub:
    """
    printf 'chr1\\t100\\t300\\tAAA-1\\t1\\n' | gzip > '${meta.id}.source.fragments.gz'
    printf 'barcode\\tgroup\\nAAA-1\\tall\\n' > '${meta.id}.source.cells.tsv'
    """
}
