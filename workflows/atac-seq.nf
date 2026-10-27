nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { GET_REFERENCES } from '../modules/common/get_references'
include { DOWNLOAD } from '../modules/common/download'
include { FASTQC } from '../modules/common/fastqc'
include { ATAC_TRIM } from '../modules/atac-seq/atac_trim'
include { ATAC_ALIGN } from '../modules/atac-seq/atac_align'
include { ATAC_FILTER } from '../modules/atac-seq/atac_filter'
include { ATAC_PEAKS } from '../modules/atac-seq/atac_peaks'
include { ATAC_TRACKS } from '../modules/atac-seq/atac_tracks'
include { ATAC_QC } from '../modules/atac-seq/atac_qc'
include { MULTIQC_REPORT } from '../modules/common/multiqc_report'

// Sheet rules a JSON schema cannot express: the exact columns, existing reference files, distinct
// reference IDs, and one gene annotation and organelle list per reference.
def checkSamples(sheet, rows, references) {
    def header = ['sample_id', 'species', 'read1_url', 'read2_url', 'reference_fasta',
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

    // assets/atac-seq.params.schema.json checks the sheet, paths, resources and analysis options.
    validateParameters()
    log.info(paramsSummaryLog(workflow))
    def sheet = file(params.input)
    def referencesDir = file(params.references).toAbsolutePath()
    def resources = [cpus: params.cpus, memory: params.memory, time: params.time]
    def bwaOptions = [
        batch_size: params.bwa_batch_size, seed_length: params.bwa_seed_length,
        max_seed_occurrences: params.bwa_max_seed_occurrences,
        score_threshold: params.bwa_score_threshold
    ]
    // Validate the whole sheet before any task starts, so a bad sheet downloads nothing.
    def sheetRows = samplesheetToList(sheet, "${repoDir}/assets/atac-seq.samples.schema.json")
        .collect { row -> row[0] }
    checkSamples(sheet, sheetRows, referencesDir)
    samples = channel.fromList(sheetRows.collect { row ->
        row + [ref_id: row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, ''), layout: 'PE',
               organelles: row.organelle_contigs == '-' ? [] : row.organelle_contigs.tokenize(',')]
    })
    // Check each reference once; downstream steps rejoin samples to it by ref_id.
    references = samples.unique { sample -> sample.ref_id }
        .map { meta ->
            def ref = [id: meta.ref_id, filename: meta.reference_fasta]
            tuple(ref, file("${referencesDir}/${ref.filename}", checkIfExists: true))
        }
    GET_REFERENCES(references, referencesDir)
    sizes_ch = GET_REFERENCES.out.sizes
    indexes = GET_REFERENCES.out.indexes
    DOWNLOAD(samples.collate(params.download_batch_size))
    // Split each batch back into one tuple per sample, with both mates.
    downloaded = DOWNLOAD.out.downloads.flatMap { batch, reads, _counts ->
        batch.collect { meta ->
            def names = ['read1', 'read2'].collect { mate -> "${meta.id}.${mate}.fastq.gz".toString() }
            tuple(meta, reads.findAll { read -> read.name in names }.sort { read -> read.name })
        }
    }
    // Observe each raw mate independently; trimming continues from the original reads.
    FASTQC(downloaded.flatMap { meta, reads ->
        reads.collect { read ->
            tuple(meta, read.name.endsWith('.read2.fastq.gz') ? 'read2' : 'read1', read)
        }
    })
    ATAC_TRIM(downloaded, params.atac_adapter)
    // Alignment memory follows the whole genome; MACS3's genome size counts nuclear contigs only.
    sized = ATAC_TRIM.out.reads.map { meta, reads -> tuple(meta.ref_id, meta, reads) }
        .combine(sizes_ch, by: 0)
        .map { refId, meta, reads, _ref, path ->
            def lengths = path.readLines().collect { line -> line.tokenize('\t') }
            def nuclear = lengths.findAll { fields -> !(fields[0] in meta.organelles) }
            tuple(refId, meta + [
                genome_size: lengths.sum { fields -> fields[1].toLong() },
                nuclear_size: nuclear.sum(0L) { fields -> fields[1].toLong() }
            ], reads)
        }
    alignmentInputs = sized.combine(indexes, by: 0)
        .map { _refId, meta, reads, _ref, index -> tuple(meta, reads, index) }
    ATAC_ALIGN(alignmentInputs, bwaOptions)
    filterInputs = ATAC_ALIGN.out.bam.map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
        .combine(sizes_ch, by: 0)
        .map { _refId, meta, bam, bai, _ref, path -> tuple(meta, bam, bai, path) }
    ATAC_FILTER(filterInputs, params.atac_mapq)
    usable = ATAC_FILTER.out.bam
    ATAC_PEAKS(usable)
    ATAC_TRACKS(usable, resources)
    qcInputs = usable.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(ATAC_PEAKS.out.peaks.map { meta, peaks -> tuple(meta.id, peaks) },
              failOnDuplicate: true, failOnMismatch: true)
        .join(ATAC_TRACKS.out.cuts.map { meta, cuts -> tuple(meta.id, cuts) },
              failOnDuplicate: true, failOnMismatch: true)
        .join(ATAC_FILTER.out.complexity.map { meta, complexity -> tuple(meta.id, complexity) },
              failOnDuplicate: true, failOnMismatch: true)
        .map { _id, meta, bam, bai, peaks, cuts, complexity ->
            def genes = meta.reference_genes == '-' ? [] : file("${referencesDir}/${meta.reference_genes}")
            tuple(meta, bam, bai, peaks, cuts, complexity, genes)
        }
    ATAC_QC(qcInputs, resources)
    // Collect channels, not publishDir contents: wait for every library's reports.
    multiqcFiles = FASTQC.out.reports.map { _meta, _mate, _html, zip -> zip }
        .mix(ATAC_TRIM.out.report.map { _meta, report -> report })
        .mix(ATAC_ALIGN.out.reports.map { _meta, reports -> reports })
        .mix(ATAC_FILTER.out.reports.map { _meta, report -> report })
        .mix(ATAC_QC.out.summary.map { _meta, summary -> summary })
        .mix(ATAC_QC.out.tss.map { _meta, files -> files })
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
    file("${repoDir}/conf/dap-images.lock").copyTo(runOutput.resolve('dap-images.lock'))
    runOutput.resolve('git-status.txt').text = gitReport
    MULTIQC_REPORT(multiqcFiles)
}
