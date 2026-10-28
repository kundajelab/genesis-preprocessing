nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { GET_REFERENCES } from '../modules/common/get_references'
include { MULTIQC_REPORT } from '../modules/common/multiqc_report'
include { SCATAC_DOWNLOAD } from '../modules/scatac-seq/scatac_download'
include { SCATAC_FRAGMENTS } from '../modules/scatac-seq/scatac_fragments'
include { SCATAC_PEAKS } from '../modules/scatac-seq/scatac_peaks'
include { SCATAC_TRACKS } from '../modules/scatac-seq/scatac_tracks'
include { SCATAC_QC } from '../modules/scatac-seq/scatac_qc'

// Sheet rules a JSON schema cannot express: the exact columns, existing reference files, distinct
// reference IDs, and one gene annotation and organelle list per reference.
def checkSamples(sheet, rows, references) {
    def header = ['sample_id', 'species', 'fragments_url', 'cells_url', 'reference_fasta',
                  'reference_genes', 'organelle_contigs']
    def columns = sheet.withReader { reader -> reader.readLine() }?.split('\t')?.toList()
    if (columns != header) {
        error("${sheet}: expected header ${header.join(' ')}")
    }
    def referenceIds = [:]
    def annotations = [:]
    rows.each { row ->
        if (!references.resolve(row.reference_fasta).isFile()) {
            error("${sheet}: reference missing: ${row.reference_fasta}")
        }
        if (row.reference_genes != '-' && !references.resolve(row.reference_genes).isFile()) {
            error("${sheet}: gene annotation missing: ${row.reference_genes}")
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

    // assets/scatac-seq.params.schema.json checks the sheet, paths and resources.
    validateParameters()
    log.info(paramsSummaryLog(workflow))
    def sheet = file(params.input)
    def referencesDir = file(params.references).toAbsolutePath()
    def resources = [cpus: params.cpus, memory: params.memory, time: params.time]
    // Validate the whole sheet before any task starts, so a bad sheet downloads nothing.
    def sheetRows = samplesheetToList(sheet, "${repoDir}/assets/scatac-seq.samples.schema.json")
        .collect { row -> row[0] }
    checkSamples(sheet, sheetRows, referencesDir)
    samples = channel.fromList(sheetRows.collect { row ->
        row + [ref_id: row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, ''),
               organelles: row.organelle_contigs == '-' ? [] : row.organelle_contigs.tokenize(',')]
    })
    // Only chromosome sizes are needed, as the fragments are already aligned, but the shared
    // reference step reuses existing products and leaves the bwa-mem2 index for other modalities.
    references = samples.unique { sample -> sample.ref_id }
        .map { meta ->
            def ref = [id: meta.ref_id, filename: meta.reference_fasta]
            tuple(ref, file("${referencesDir}/${ref.filename}", checkIfExists: true))
        }
    GET_REFERENCES(references, referencesDir)
    sizes_ch = GET_REFERENCES.out.sizes.map { refId, _ref, path -> tuple(refId, path) }
    SCATAC_DOWNLOAD(samples)
    // MACS3's genome size counts nuclear contigs only, as in bulk ATAC.
    sized = SCATAC_DOWNLOAD.out.downloads
        .map { meta, fragments, cells -> tuple(meta.ref_id, meta, fragments, cells) }
        .combine(sizes_ch, by: 0)
        .map { _refId, meta, fragments, cells, path ->
            def nuclear = path.readLines().collect { line -> line.tokenize('\t') }
                .findAll { fields -> !(fields[0] in meta.organelles) }
            tuple(meta + [nuclear_size: nuclear.sum(0L) { fields -> fields[1].toLong() }],
                  fragments, cells, path)
        }
    SCATAC_FRAGMENTS(sized)
    // One pseudobulk per library and group; libraries, and so biological replicates, stay separate.
    groups = SCATAC_FRAGMENTS.out.groups.flatMap { meta, files ->
        (files instanceof List ? files : [files]).collect { bed ->
            def group = bed.name.substring(meta.id.length() + 1).replaceFirst(/\.fragments\.bed$/, '')
            tuple(meta + [library: meta.id, group: group, id: "${meta.id}.${group}".toString()], bed)
        }
    }
    SCATAC_PEAKS(groups)
    SCATAC_TRACKS(groups.map { meta, bed -> tuple(meta.ref_id, meta, bed) }
        .combine(sizes_ch, by: 0)
        .map { _refId, meta, bed, path -> tuple(meta, bed, path) }, resources)
    qcInputs = groups.map { meta, bed -> tuple(meta.id, meta, bed) }
        .join(SCATAC_TRACKS.out.bam.map { meta, bam, bai -> tuple(meta.id, bam, bai) },
              failOnDuplicate: true, failOnMismatch: true)
        .join(SCATAC_PEAKS.out.peaks.map { meta, peaks -> tuple(meta.id, peaks) },
              failOnDuplicate: true, failOnMismatch: true)
        .join(SCATAC_TRACKS.out.cuts.map { meta, cuts -> tuple(meta.id, cuts) },
              failOnDuplicate: true, failOnMismatch: true)
        .map { _id, meta, bed, bam, bai, peaks, cuts ->
            def genes = meta.reference_genes == '-' ? [] : file("${referencesDir}/${meta.reference_genes}")
            tuple(meta, bed, bam, bai, peaks, cuts, genes)
        }
    SCATAC_QC(qcInputs, resources)
    // Collect channels, not publishDir contents: wait for every library's reports.
    multiqcFiles = SCATAC_FRAGMENTS.out.summary.map { _meta, summary -> summary }
        .mix(SCATAC_QC.out.summary.map { _meta, summary -> summary })
        .mix(SCATAC_QC.out.tss.map { _meta, files -> files })
        .flatten().collect().map { paths -> paths.sort { it.name } }
    // Record what produced this run beside its outputs, as workflows/dap-seq.nf does.
    def git = { List args ->
        try {
            def process = (['git', '-C', repoDir.toString()] + args).execute()
            // Trim only the end: porcelain status lines start with positional status columns.
            def text = process.text.stripTrailing()
            return process.waitFor() == 0 ? text : null
        } catch (IOException ignored) {
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
