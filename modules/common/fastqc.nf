include { dotenv } from 'plugin/nf-dotenv'

process FASTQC {
    cpus 1
    memory '1 GB'
    time '4h'
    conda "${moduleDir}/../../environments/DAP_SEQ_FASTQC.yaml"
    container "${dotenv('DAP_SEQ_FASTQC_IMAGE')}"
    tag "${meta.id}:${mate}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc/fastqc", mode: 'copy', pattern: '*fastqc*'
    input:
    tuple val(meta), val(mate), path(read)
    output:
    tuple val(meta), val(mate), path("${meta.id}.${mate}_fastqc.html"),
        path("${meta.id}.${mate}_fastqc.zip"), emit: reports
    tuple val(meta), val(mate), path("${meta.id}.${mate}.fastqc.version.txt"), emit: version
    script:
    """
    # Fontconfig must not attempt a cache under the container's unwritable default HOME.
    export XDG_CACHE_HOME="\$PWD/.cache"
    mkdir -p "\$XDG_CACHE_HOME"
    fastqc --version > '${meta.id}.${mate}.fastqc.version.txt'
    fastqc --threads ${task.cpus} --memory 512 --format fastq --noextract --outdir . '${read}'
    # Require usable report artifacts, not PASS labels: assay-specific criteria are unspecified.
    test -s '${meta.id}.${mate}_fastqc.html'
    test -s '${meta.id}.${mate}_fastqc.zip'
    """
    stub:
    """
    printf '<html>stub</html>\\n' > '${meta.id}.${mate}_fastqc.html'
    printf 'stub\\n' > '${meta.id}.${mate}_fastqc.zip'
    printf 'stub\\n' > '${meta.id}.${mate}.fastqc.version.txt'
    """
}
