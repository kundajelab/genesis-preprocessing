nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { GET_REFERENCES } from '../modules/common/get_references'
include { MULTIQC_REPORT } from '../modules/common/multiqc_report'
include { SCATAC_SRA } from '../modules/scatac-seq/scatac_sra'
include { SCATAC_CHROMAP } from '../modules/scatac-seq/scatac_chromap'
include { SCATAC_CELLS } from '../modules/scatac-seq/scatac_cells'
include { SCATAC_PSEUDOBULKS } from '../modules/scatac-seq/scatac_pseudobulks'

// The barcode read's length for each chemistry: the 16-base cell barcode alone, or after an
// 8-base spacer in multiome ATAC.
def barcodeLength(chemistry) {
    return ['10x-atac-v1': 16, '10x-atac-v2': 16, '10x-multiome-atac': 24][chemistry]
}

// Sheet rules a JSON schema cannot express: the exact columns, existing reference and whitelist
// files, distinct reference IDs, one gene annotation and organelle list per reference, and runs
// used by one library only.
def checkSamples(sheet, rows, references) {
    def header = ['sample_id', 'species', 'sra_runs', 'chemistry', 'cells_url', 'reference_fasta',
                  'reference_genes', 'organelle_contigs', 'barcode_whitelist']
    def columns = sheet.withReader { reader -> reader.readLine() }?.split('\t')?.toList()
    if (columns != header) {
        error("${sheet}: expected header ${header.join(' ')}")
    }
    def referenceIds = [:]
    def annotations = [:]
    def runs = [:]
    rows.each { row ->
        if (!references.resolve(row.reference_fasta).isFile()) {
            error("${sheet}: reference missing: ${row.reference_fasta}")
        }
        if (row.reference_genes != '-' && !references.resolve(row.reference_genes).isFile()) {
            error("${sheet}: gene annotation missing: ${row.reference_genes}")
        }
        if (row.barcode_whitelist != '-' && !references.resolve(row.barcode_whitelist).isFile()) {
            error("${sheet}: barcode whitelist missing: ${row.barcode_whitelist}")
        }
        row.sra_runs.tokenize(',').each { run ->
            def owner = runs.putIfAbsent(run, row.id)
            if (owner != null) {
                error("${sheet}: run ${run} is listed for both ${owner} and ${row.id}")
            }
        }
        def referenceId = row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, '')
        def previous = referenceIds.putIfAbsent(referenceId, row.reference_fasta)
        if (previous != null && previous != row.reference_fasta) {
            error("${sheet}: reference ID collision: ${previous} and ${row.reference_fasta} " +
                  "both use ${referenceId}; use distinct reference basenames")
        }
        def annotation = [row.reference_genes, row.organelle_contigs]
        def recorded = annotations.putIfAbsent(row.reference_fasta, annotation)
        if (recorded != null && recorded != annotation) {
            error("${sheet}: rows using ${row.reference_fasta} disagree on reference_genes or " +
                  "organelle_contigs")
        }
    }
}

workflow {
    // This script lives in workflows/; everything else it uses is relative to the repository.
    def repoDir = projectDir.parent
    // The same end-of-log banner as workflows/dap-seq.nf; see the notes there.
    workflow.onComplete = {
        def rule = '=' * 78
        def report = [
            rule,
            "PIPELINE ${workflow.success ? 'SUCCEEDED' : 'FAILED'}  (${workflow.runName})",
            "  duration: ${workflow.duration}",
            "  tasks:    ${workflow.stats.succeedCount} succeeded, ${workflow.stats.cachedCount} cached, " +
                "${workflow.stats.failedCount} failed, ${workflow.stats.ignoredCount} ignored",
            "  exit:     ${workflow.exitStatus}",
            "  output:   ${params.run_dir}/output",
            rule,
        ]
        log.info(report.join('\n'))
    }

    // assets/scatac-reads.params.schema.json checks the sheet, paths and resources.
    validateParameters()
    log.info(paramsSummaryLog(workflow))
    def sheet = file(params.input)
    def referencesDir = file(params.references).toAbsolutePath()
    def resources = [cpus: params.cpus, memory: params.memory, time: params.time]
    // Validate the whole sheet before any task starts, so a bad sheet downloads nothing.
    def sheetRows = samplesheetToList(sheet, "${repoDir}/assets/scatac-reads.samples.schema.json")
        .collect { row -> row[0] }
    checkSamples(sheet, sheetRows, referencesDir)
    samples = channel.fromList(sheetRows.collect { row ->
        row + [ref_id: row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, ''),
               organelles: row.organelle_contigs == '-' ? [] : row.organelle_contigs.tokenize(','),
               runs: row.sra_runs.tokenize(','), barcode_length: barcodeLength(row.chemistry)]
    })
    // chromap indexes the FASTA itself; the shared reference step gives chromosome sizes, reusing
    // existing products, and leaves the bwa-mem2 index for other modalities.
    references = samples.unique { sample -> sample.ref_id }
        .map { meta ->
            def ref = [id: meta.ref_id, filename: meta.reference_fasta]
            tuple(ref, file("${referencesDir}/${ref.filename}", checkIfExists: true))
        }
    GET_REFERENCES(references, referencesDir)
    sizes_ch = GET_REFERENCES.out.sizes.map { refId, _ref, path -> tuple(refId, path) }
    // One task per run (lane); a library's runs are mapped together.
    SCATAC_SRA(samples.flatMap { meta -> meta.runs.collect { run -> tuple(meta, run) } })
    libraries = SCATAC_SRA.out.reads
        .map { meta, run, r1, r2, barcodes ->
            tuple(groupKey(meta.id, meta.runs.size()), meta, run, r1, r2, barcodes)
        }
        .groupTuple()
        .map { _id, metas, runs, r1, r2, barcodes ->
            def meta = metas[0]
            // Keep the sheet's run order, so the inputs and so the outputs are reproducible.
            def order = meta.runs.collect { run -> runs.indexOf(run) }
            tuple(meta.ref_id, meta, order.collect { i -> r1[i] }, order.collect { i -> r2[i] },
                  order.collect { i -> barcodes[i] })
        }
        .combine(sizes_ch, by: 0)
        .map { _refId, meta, r1, r2, barcodes, sizes ->
            def length = sizes.readLines().sum(0L) { line -> line.tokenize('\t')[1].toLong() }
            def whitelist = meta.barcode_whitelist == '-' ? [] :
                file("${referencesDir}/${meta.barcode_whitelist}")
            tuple(meta + [genome_length: length], r1, r2, barcodes,
                  file("${referencesDir}/${meta.reference_fasta}"), whitelist)
        }
    SCATAC_CHROMAP(libraries)
    SCATAC_CELLS(SCATAC_CHROMAP.out.fragments)
    SCATAC_PSEUDOBULKS(SCATAC_CELLS.out.cells, sizes_ch, referencesDir, resources)
    // Collect channels, not publishDir contents: wait for every library's reports.
    multiqcFiles = SCATAC_PSEUDOBULKS.out.multiqc
        .flatten().collect().map { paths -> paths.sort { path -> path.name } }
    // Record what produced this run beside its outputs, as workflows/dap-seq.nf does.
    def git = { List args ->
        try {
            def process = (['git', '-C', repoDir.toString()] + args).execute()
            // Trim only the end: porcelain status lines start with positional status columns.
            def text = process.text.stripTrailing()
            return process.waitFor() == 0 ? text : null
        } catch (IOException _ignored) {
            return null
        }
    }
    def gitSha = git.call(['rev-parse', 'HEAD']) ?: workflow.commitId ?: 'UNKNOWN'
    def branch = git.call(['rev-parse', '--abbrev-ref', 'HEAD'])
    def changes = branch == null ? null : git.call(['status', '--porcelain'])
    if (changes == null) {
        log.warn('Git unavailable: git-status.txt is UNKNOWN')
    }
    def engine = "container engine: ${workflow.containerEngine ?: 'none'}\n"
    def gitReport = changes == null ?
        "git status unavailable: ${repoDir} is not a Git checkout, or git is not installed\n" +
            engine :
        "branch: ${branch}\ncommit: ${gitSha}\n" + engine +
            (changes ? "dirty files:\n${changes}\n" : "dirty files: none\n")
    def runOutput = file("${params.run_dir}/output")
    runOutput.mkdirs()
    file("${repoDir}/conf/genesis-images.lock").copyTo(runOutput.resolve('genesis-images.lock'))
    runOutput.resolve('git-status.txt').text = gitReport
    MULTIQC_REPORT(multiqcFiles)
}
