nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { GET_REFERENCES } from '../modules/common/get_references'
include { DOWNLOAD } from '../modules/common/download'
include { FASTQC } from '../modules/common/fastqc'
include { METADATA } from '../modules/dap-seq/metadata'
include { BWA_MEM2_ALIGN } from '../modules/dap-seq/bwa_mem2_align'
include { QC } from '../modules/dap-seq/qc'
include { MULTIQC } from '../modules/dap-seq/multiqc'
include { TRACKS } from '../modules/dap-seq/tracks'
include { CALL_PEAKS } from '../modules/dap-seq/call_peaks'
include { QUANTIFY } from '../modules/dap-seq/quantify'

// Sheet rules a JSON schema cannot express: the exact columns, existing reference FASTAs,
// distinct reference IDs and valid control assignments. Rows come from samplesheetToList.
def checkSamples(sheet, rows, references) {
    def header = ['sample_id', 'species', 'read1_url', 'read2_url', 'control_sample', 'reference_fasta']
    def columns = sheet.withReader { reader -> reader.readLine() }?.split('\t')?.toList()
    if (columns != header) {
        error("${sheet}: expected header ${header.join(' ')}")
    }
    def referenceIds = [:]
    rows.each { row ->
        if (!references.resolve(row.reference_fasta).isFile()) {
            error("${sheet}: reference missing: ${row.reference_fasta}")
        }
        def referenceId = row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, '')
        def previous = referenceIds.putIfAbsent(referenceId, row.reference_fasta)
        if (previous != null && previous != row.reference_fasta) {
            error("${sheet}: reference ID collision: ${previous} and ${row.reference_fasta} " +
                  "both use ${referenceId}; use distinct reference basenames")
        }
    }
    def samples = rows.collectEntries { row -> [(row.id): row] }
    rows.findAll { row -> row.control != '-' }.each { row ->
        def control = samples[row.control]
        if (control == null || row.control == row.id || control.control != '-' ||
                control.species != row.species || control.reference_fasta != row.reference_fasta ||
                (control.read2_url == '-') != (row.read2_url == '-')) {
            error("${sheet}: invalid control ${row.control} for ${row.id}")
        }
    }
}

// Bowtie 1 modules deliberately have no imports or calls here.
workflow {
    // This script lives in workflows/; everything else it uses is relative to the repository.
    def repoDir = projectDir.parent
    // Print an unmistakable banner as the last thing in the log. The submission script should set
    // NXF_ANSI_LOG=false so the log stays readable in a file, which leaves it as a flat stream of
    // process notices: without this it is hard to tell a run that finished from one that died, or
    // one that is still going.
    // NOTE: assigned here inside the entry workflow rather than declared in nextflow.config, where
    // workflow handlers are deprecated, or as a top-level "workflow.onComplete {}" statement, which
    // nextflow lint rejects as a statement mixed with script declarations. The "onComplete:" entry
    // workflow section would be tidier still, but it only compiles under NXF_SYNTAX_PARSER=v2, and
    // v2 cannot resolve the nf-dotenv "dotenv()" calls that every module uses to pick its container.
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
        // NOTE: deliberately no workflow.errorMessage here. process.shell sets -x, so for a failed
        // task that field holds the whole shell trace and would bury this summary. Nextflow already
        // prints the failing process, its exit status and its work dir immediately above.
        log.info(report.join('\n'))
    }

    // assets/dap-seq.params.schema.json checks the sheet, paths, resources and BWA options.
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
    def sheetRows = samplesheetToList(sheet, "${repoDir}/assets/dap-seq.samples.schema.json")
        .collect { row -> row[0] }
    checkSamples(sheet, sheetRows, referencesDir)
    samples = channel.fromList(sheetRows.collect { row ->
        row + [ref_id: row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, ''),
               layout: row.read2_url == '-' ? 'SE' : 'PE']
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
    // Split each batch back into one tuple per sample: a single read path for SE, a list for PE.
    downloaded = DOWNLOAD.out.downloads.flatMap { batch, reads, _counts ->
        def files = reads instanceof List ? reads : [reads]
        batch.collect { meta ->
            def names = ['read1', 'read2'].collect { mate -> "${meta.id}.${mate}.fastq.gz".toString() }
            def mates = files.findAll { read -> read.name in names }.sort { read -> read.name }
            tuple(meta.id, meta, mates.size() == 1 ? mates[0] : mates)
        }
    }
    // Observe each raw mate independently; alignment continues to consume the original reads.
    fastqcInputs = downloaded.flatMap { _id, meta, reads ->
        def files = reads instanceof List ? reads : [reads]
        files.collect { read ->
            def mate = read.name.endsWith('.read2.fastq.gz') ? 'read2' : 'read1'
            tuple(meta, mate, read)
        }
    }
    FASTQC(fastqcInputs)
    metadataInputs = downloaded.map { _id, meta, reads -> tuple(meta.ref_id, meta, reads) }
        .combine(sizes_ch, by: 0)
        .map { _refId, meta, reads, _ref, path -> tuple(meta, reads, path) }
    METADATA(metadataInputs)
    inferred = METADATA.out.metadata.map { meta, path ->
        def lines = path.readLines()
        def values = [lines[0].split('\t').toList(), lines[1].split('\t').toList()]
            .transpose().collectEntries { it }
        tuple(meta.id, meta + [
            genome_size: values.genome_size.toLong(),
            read_length: values.analysis_read_length.toInteger()
        ], path)
    }
    alignmentInputs = downloaded.join(inferred, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _id, _original, reads, meta, metadata -> tuple(meta.ref_id, meta, reads, metadata) }
        .combine(indexes, by: 0)
        .map { _refId, meta, reads, metadata, _ref, index -> tuple(meta, reads, metadata, index) }
    BWA_MEM2_ALIGN(alignmentInputs, bwaOptions, params.spp_read_length)
    aligned = BWA_MEM2_ALIGN.out.aligned
    mainBams = aligned.map { meta, bam, bai, _qcBam, _qcBai -> tuple(meta, bam, bai) }
    sppBams = BWA_MEM2_ALIGN.out.spp.map { meta, sppBam -> tuple(meta.id, sppBam) }
    qcInputs = aligned
        .map { meta, bam, bai, qcBam, qcBai -> tuple(meta.id, meta, bam, bai, qcBam, qcBai) }
        .join(sppBams, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _id, meta, bam, bai, qcBam, qcBai, sppBam ->
            tuple(meta, bam, bai, qcBam, qcBai, sppBam)
        }
    QC(qcInputs, resources, file(params.spp_script))
    // Collect channels, not publishDir contents: wait for every library's QC and metadata.
    multiqcFiles = FASTQC.out.reports.map { _meta, _mate, _html, zip -> zip }
        .mix(QC.out.reports.map { _meta, reports -> reports })
        .mix(METADATA.out.metadata.map { _meta, path -> path })
        .flatten().collect().map { paths -> paths.sort { it.name } }
    multiqcSamples = samples.map { meta ->
        [id: meta.id, species: meta.species, ref_id: meta.ref_id,
         reference_fasta: meta.reference_fasta, layout: meta.layout, control: meta.control]
    }.collect().map { rows -> rows.sort { it.id } }
    // Record what produced this run beside its outputs: the container image pins, and the
    // checkout's branch, commit and uncommitted changes. Local runs can have a null commitId.
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
    // MultiQC's source status ignores untracked files; git-status.txt lists them too.
    def sourceStatus = changes == null ? 'UNKNOWN' :
        (changes.readLines().any { line -> !line.startsWith('??') } ? 'dirty' : 'clean')
    if (changes == null) {
        log.warn('Git unavailable: git-status.txt and MultiQC source provenance are UNKNOWN')
    }
    // Curation treats a run as complete only from a clean checkout run in pinned containers.
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
    def reportProvenance = [name: file(params.run_dir).name, git_sha: gitSha,
                           source_status: sourceStatus, pipeline_version: workflow.manifest.version]
    MULTIQC(multiqcFiles, multiqcSamples, reportProvenance,
            file("${repoDir}/multiqc_config.yaml"),
            file("${repoDir}/genesis_tools/src/genesis_tools/multiqc_reporting.py"))
    TRACKS(aligned, resources)
    treatments = mainBams.filter { meta, _bam, _bai -> meta.control != '-' }
        .map { meta, bam, bai -> tuple(meta.control, meta, bam, bai) }
    controls = mainBams.filter { meta, _bam, _bai -> meta.control == '-' }
        .map { meta, bam, bai -> tuple(meta.id, bam, bai) }
    peakInputs = treatments.combine(controls, by: 0)
        .map { _controlId, meta, bam, bai, controlBam, controlBai ->
            tuple(meta, bam, bai, controlBam, controlBai)
        }
    CALL_PEAKS(peakInputs)
    peaksBySample = CALL_PEAKS.out.peaks.map { meta, peaks -> tuple(meta.id, meta, peaks) }
    bamsBySample = mainBams.map { meta, bam, bai -> tuple(meta.id, bam, bai) }
    coverageBySample = TRACKS.out.rpkm.map { meta, coverage -> tuple(meta.id, coverage) }
    quantificationInputs = peaksBySample.join(bamsBySample, by: 0, failOnDuplicate: true)
        .join(coverageBySample, by: 0, failOnDuplicate: true)
        .map { _id, meta, peaks, bam, _bai, coverage ->
            tuple(meta.ref_id, meta, peaks, bam, coverage)
        }
        .combine(sizes_ch, by: 0)
        .map { _refId, meta, peaks, bam, coverage, _ref, path ->
            tuple(meta, peaks, bam, coverage, path)
        }
    QUANTIFY(quantificationInputs)
}
