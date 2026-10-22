include { dotenv } from 'plugin/nf-dotenv'

process MULTIQC {
    cpus 1
    memory '2 GB'
    time '1h'
    conda "${moduleDir}/../environments/DAP_SEQ_MULTIQC.yaml"
    container "${dotenv('DAP_SEQ_MULTIQC_IMAGE')}"
    tag "${provenance.name}"
    publishDir "${params.run_dir}/output/multiqc", mode: 'copy', pattern: 'multiqc*'
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
