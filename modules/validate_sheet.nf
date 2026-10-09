include { dotenv } from 'plugin/nf-dotenv'

process VALIDATE_SHEET {
    cpus 1
    memory '1 GB'
    time '1h'
    conda "environments/GENESIS_TOOLS.yaml"
    container "${dotenv('GENESIS_TOOLS_IMAGE')}"
    input:
    path sheet
    path references
    path toolsSource
    output:
    path 'validated.tsv', emit: sheet
    script:
    """
    export PYTHONPATH='${toolsSource}'
    export PYTHONDONTWRITEBYTECODE=1
    genesis-tools validate-sheet --sheet '${sheet}' \
        --references '${references}' --output validated.tsv
    """
    stub:
    """
    export PYTHONPATH='${toolsSource}'
    export PYTHONDONTWRITEBYTECODE=1
    genesis-tools validate-sheet --sheet '${sheet}' \
        --references '${references}' --output validated.tsv
    """
}
