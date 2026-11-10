include { dotenv } from 'plugin/nf-dotenv'
include { oomMemoryOf ; oomErrorStrategy } from '../common/retry'

process MULTIQC {
    // On Perlmutter, a report of 142 Sorghum samples took 32-35 seconds at 91% CPU, peaking at
    // 0.47 GB. Allow 1 GB plus 10 MB per sample, and 15 minutes plus 2 seconds per sample.
    cpus 1
    memory { oomMemoryOf(task, 1.GB + 10.MB * samples.size()) }
    time { 15.min + 2.s * samples.size() }
    // Retry with twice the memory only when an attempt ran out of memory (modules/common/retry.nf).
    errorStrategy { oomErrorStrategy(task) }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_MULTIQC.yaml"
    container "${dotenv('GENESIS_MULTIQC_IMAGE')}"
    tag "${provenance.name}"
    publishDir "${params.publish_base}/output/multiqc", mode: 'copy', pattern: 'multiqc*'
    input:
    path reports, stageAs: 'inputs/*'
    val samples
    val provenance
    path config
    path driver
    output:
    path 'multiqc_report.html', emit: html
    path 'multiqc_data', emit: data
    script:
    def manifestJson = groovy.json.JsonOutput.toJson([run: provenance, samples: samples])
    """
    cat > genesis_manifest.json << 'GENESIS_MANIFEST'
${manifestJson}
GENESIS_MANIFEST
    python '${driver}' --manifest genesis_manifest.json --inputs inputs \\
        --config '${config}' --output .
    multiqc --version > multiqc_data/genesis_multiqc_version.txt
    """
    stub:
    def manifestJson = groovy.json.JsonOutput.toJson([run: provenance, samples: samples])
    """
    mkdir multiqc_data
    printf '<html>MultiQC stub</html>\\n' > multiqc_report.html
    cat > multiqc_data/genesis_provenance.json << 'GENESIS_MANIFEST'
${manifestJson}
GENESIS_MANIFEST
    """
}
