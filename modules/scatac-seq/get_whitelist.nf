include { dotenv } from 'plugin/nf-dotenv'

// A 10x barcode list missing from the references folder, downloaded from 10x's public source
// and checked against its pinned SHA256. The references folder is git-ignored and the list is
// used only to correct barcodes, never redistributed.
process GET_WHITELIST {
    cpus 1
    memory '1 GB'
    time { 30.min * task.attempt }
    // A transfer can stall or fail; fetch it again.
    errorStrategy 'retry'
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_ALIGNMENT.yaml"
    container "${dotenv('GENESIS_ALIGNMENT_IMAGE')}"
    tag "${name}"
    // Stub runs never publish into the shared references folder: their fake list would be reused.
    publishDir params.references, mode: 'copy', overwrite: true, enabled: !workflow.stubRun
    input:
    tuple val(name), val(url), val(sha256)
    output:
    tuple val(name), path(name), emit: whitelist
    script:
    """
    aria2c --max-connection-per-server=4 --retry-wait=5 --out=download.txt '${url}'
    printf '%s  download.txt\\n' '${sha256}' | sha256sum --check --quiet || {
        printf '%s from %s does not match its pinned SHA256\\n' '${name}' '${url}' >&2
        exit 1
    }
    mv download.txt '${name}'
    """
    stub:
    """
    printf 'ACGTACGTACGTACGT\\n' > '${name}'
    """
}
