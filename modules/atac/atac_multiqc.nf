include { dotenv } from 'plugin/nf-dotenv'

process ATAC_MULTIQC {
    cpus 1
    memory '2 GB'
    time '1h'
    conda "${moduleDir}/../../environments/DAP_SEQ_MULTIQC.yaml"
    container "${dotenv('DAP_SEQ_MULTIQC_IMAGE')}"
    publishDir "${params.run_dir}/output/multiqc", mode: 'copy', pattern: 'multiqc*'
    input:
    path reports, stageAs: 'inputs/*'
    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_data', emit: data
    script:
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    printf 'data_dir_name: multiqc_data\\nno_version_check: true\\nno_ai: true\\n' > multiqc_config.yaml
    multiqc inputs --config multiqc_config.yaml --filename multiqc_report.html \\
        --data-dir --data-format json
    multiqc --version > multiqc_data/multiqc_version.txt
    """
    stub:
    """
    mkdir multiqc_data
    printf '<html>MultiQC stub</html>\\n' > multiqc_report.html
    """
}
