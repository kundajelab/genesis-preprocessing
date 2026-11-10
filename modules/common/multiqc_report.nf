include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomErrorStrategy } from './retry'

process MULTIQC_REPORT {
    // As DAP-seq's MultiQC, which read 1,422 reports in 35 seconds at 91% CPU, peaking at 0.47 GB.
    // Allow 1 GB plus 1 MB per report, and 15 minutes plus 0.2 seconds per report.
    cpus 1
    memory { oomMemoryOf(task, 1.GB + 1.MB * [reports].flatten().size()) }
    time { 15.min + 200.ms * [reports].flatten().size() }
    // Retry with twice the memory only when an attempt ran out of memory (modules/common/retry.nf).
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_MULTIQC.yaml"
    container "${dotenv('GENESIS_MULTIQC_IMAGE')}"
    publishDir "${params.publish_base}/output/multiqc", mode: 'copy', pattern: 'multiqc*'
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
