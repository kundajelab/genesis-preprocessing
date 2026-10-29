include { dotenv } from 'plugin/nf-dotenv'

// NCBI SRA keeps the 10x cell-barcode read as a technical read; ENA's FASTQs drop it. Downloading
// and converting are separate tasks, so a failed conversion, or -resume after one, reuses the
// download instead of fetching the run again.
process SRA_PREFETCH {
    // A few runs at a time keep the load on NCBI bounded.
    maxForks 2
    cpus 1
    memory '2 GB'
    time { 4.h * task.attempt }
    // A transfer can stall or fail; fetch the run again, and on a timeout with more time.
    errorStrategy 'retry'
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_SRA.yaml"
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${meta.id}:${run}"
    input:
    tuple val(meta), val(run)
    output:
    tuple val(meta), val(run), path("${run}.sra"), emit: sra
    script:
    """
    # sra-tools keeps its settings outside the home folder.
    export NCBI_SETTINGS="\$PWD/ncbi-settings.mkfg"
    prefetch --max-size u --output-directory . '${run}'
    mv '${run}/${run}.sra' '${run}.sra'
    rm -r '${run}'
    """
    stub:
    """
    touch '${run}.sra'
    """
}
