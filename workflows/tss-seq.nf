nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { DOWNLOAD } from '../modules/common/download'
include { FASTQC } from '../modules/common/fastqc'
include { GET_REFERENCES } from '../modules/common/get_references'
include { MULTIQC_REPORT } from '../modules/common/multiqc_report'
include { RNA_ANNOTATION } from '../modules/rna-seq/rna_annotation'
include { GET_STAR_INDEX } from '../modules/rna-seq/get_star_index'
include { TSS_EXTRACT } from '../modules/tss-seq/tss_extract'
include { TSS_TRIM } from '../modules/tss-seq/tss_trim'
include { TSS_ALIGN } from '../modules/tss-seq/tss_align'
include { TSS_CTSS } from '../modules/tss-seq/tss_ctss'
include { TSS_CLUSTERS } from '../modules/tss-seq/tss_clusters'
include { TSS_QC } from '../modules/tss-seq/tss_qc'

// Settings each method implies. Cap-trapping CAGE adds a non-templated G before the first base
// during reverse transcription; oligo-capping (TSS-seq) and csRNA-seq ligate an adapter to the
// RNA's own 5' end. Template-switching methods (RAMPAGE, Smar2C2) add a run of Gs, which local
// alignment soft-clips. Smar2C2 reads start with a spacer, the Smart-seq2 primer, a 12-base UMI
// and AT, in either mate; RAMPAGE's TSS read is the library's read 2.
def methodDefaults(method) {
    def cage = [g_correction: true, soft_clip_5p: false, umi_length: 0, min_length: 20, mates: 1]
    def ligation = cage + [g_correction: false]
    return [
        'CAGE': cage,
        'nAnT-iCAGE': cage,
        'TSS-seq': ligation,
        'oligo-capping': ligation,
        'csRNA-seq': ligation + [min_length: 18],
        'RAMPAGE': ligation + [soft_clip_5p: true],
        'Smar2C2': ligation + [soft_clip_5p: true, umi_length: 12, mates: 2],
    ][method]
}

// Sheet rules a JSON schema cannot express: the exact columns, existing reference files, distinct
// reference IDs, and one annotation and organelle list per reference.
def checkSamples(sheet, rows, references) {
    def header = ['sample_id', 'species', 'read1_url', 'read2_url', 'method', 'trim_5p', 'adapter_3p',
                  'reference_fasta', 'reference_annotation', 'organelle_contigs']
    def columns = sheet.withReader { reader -> reader.readLine() }?.split('\t')?.toList()
    if (columns != header) {
        error("${sheet}: expected header ${header.join(' ')}")
    }
    def referenceIds = [:]
    def annotations = [:]
    rows.each { row ->
        // Smar2C2 needs both mates (its TSS reads are in either); the other methods read one file.
        if ((methodDefaults(row.method).mates == 2) != (row.read2_url != '-')) {
            error("${sheet}: ${row.id}: ${row.method} takes " +
                  (methodDefaults(row.method).mates == 2 ? 'both read1_url and read2_url' :
                   'one FASTQ in read1_url and read2_url -'))
        }
        if (!references.resolve(row.reference_fasta).isFile()) {
            error("${sheet}: reference missing: ${row.reference_fasta}")
        }
        if (!references.resolve(row.reference_annotation).isFile()) {
            error("${sheet}: annotation missing: ${row.reference_annotation}")
        }
        def referenceId = row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, '')
        def previous = referenceIds.putIfAbsent(referenceId, row.reference_fasta)
        if (previous != null && previous != row.reference_fasta) {
            error("${sheet}: reference ID collision: ${previous} and ${row.reference_fasta} " +
                  "both use ${referenceId}; use distinct reference basenames")
        }
        def annotation = [row.reference_annotation, row.organelle_contigs]
        def recorded = annotations.putIfAbsent(row.reference_fasta, annotation)
        if (recorded != null && recorded != annotation) {
            error("${sheet}: rows using ${row.reference_fasta} disagree on reference_annotation " +
                  "or organelle_contigs")
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


    // assets/tss-seq.params.schema.json checks the sheet, paths and resources.
    validateParameters()
    log.info(paramsSummaryLog(workflow))
    def sheet = file(params.input)
    def referencesDir = file(params.references).toAbsolutePath()
    // Validate the whole sheet before any task starts, so a bad sheet downloads nothing.
    def sheetRows = samplesheetToList(sheet, "${repoDir}/assets/tss-seq.samples.schema.json")
        .collect { row -> row[0] }
    checkSamples(sheet, sheetRows, referencesDir)
    samples = channel.fromList(sheetRows.collect { row ->
        row + methodDefaults(row.method) + [
            ref_id: row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, ''),
            layout: row.read2_url == '-' ? 'SE' : 'PE',
            organelles: row.organelle_contigs == '-' ? [] : row.organelle_contigs.tokenize(',')]
    })
    // One annotation and one STAR index per reference; samples rejoin them by ref_id.
    references = samples.unique { sample -> sample.ref_id }
        .map { meta ->
            def ref = [id: meta.ref_id, filename: meta.reference_fasta, annotation: meta.reference_annotation]
            tuple(ref, file("${referencesDir}/${ref.filename}", checkIfExists: true),
                  file("${referencesDir}/${meta.reference_annotation}", checkIfExists: true))
        }
    // The shared reference step reuses existing products; STAR needs only the sizes, but the
    // bwa-mem2 index will be used by the other sequencing modality pipelines.
    GET_REFERENCES(references.map { ref, fasta, _annotation -> tuple(ref, fasta) }, referencesDir)
    // Keep the ref GET_REFERENCES emits: its FASTA hash keys the STAR index record.
    sized = references.map { ref, fasta, annotation -> tuple(ref.id, fasta, annotation) }
        .join(GET_REFERENCES.out.sizes, failOnDuplicate: true, failOnMismatch: true)
        .map { _id, fasta, annotation, ref, path ->
            def lengths = path.readLines().collect { line -> line.tokenize('\t')[1].toLong() }
            tuple(ref + [length: lengths.sum(), contigs: lengths.size()], fasta, annotation, path)
        }
    RNA_ANNOTATION(sized.map { ref, _fasta, annotation, path -> tuple(ref, annotation, path) })
    annotations = RNA_ANNOTATION.out.annotation.map { ref, gtf, rrna, maxIntron ->
        tuple(ref.id, ref, gtf, rrna, maxIntron.text.trim().toInteger())
    }
    GET_STAR_INDEX(sized.map { ref, fasta, _annotation, _path -> tuple(ref.id, ref, fasta) }
        .join(annotations, failOnDuplicate: true, failOnMismatch: true)
        .map { _id, ref, fasta, _ref, gtf, _rrna, _maxIntron -> tuple(ref, fasta, gtf) }, referencesDir)
    DOWNLOAD(samples.collate(params.download_batch_size))
    // Split each batch back into one tuple per sample: a single read path for SE, a list for PE.
    downloaded = DOWNLOAD.out.downloads.flatMap { batch, reads, _counts ->
        def files = reads instanceof List ? reads : [reads]
        batch.collect { meta ->
            def names = ['read1', 'read2'].collect { mate -> "${meta.id}.${mate}.fastq.gz".toString() }
            def mates = files.findAll { read -> read.name in names }.sort { read -> read.name }
            tuple(meta, mates.size() == 1 ? mates[0] : mates)
        }
    }
    // Observe each raw mate independently; STAR aligns the original reads.
    FASTQC(downloaded.flatMap { meta, reads ->
        [reads].flatten().collect { read ->
            tuple(meta, read.name.endsWith('.read2.fastq.gz') ? 'read2' : 'read1', read)
        }
    })
    // Smar2C2's TSS reads are picked out of both mates first; the other methods' reads start at
    // the TSS already.
    byMethod = downloaded.branch { meta, _reads ->
        extract: meta.mates == 2
        direct: true
    }
    TSS_EXTRACT(byMethod.extract)
    TSS_TRIM(byMethod.direct.mix(TSS_EXTRACT.out.reads))
    alignmentInputs = TSS_TRIM.out.reads.map { meta, reads -> tuple(meta.ref_id, meta, reads) }
        .combine(annotations, by: 0)
        .combine(GET_STAR_INDEX.out.index.map { refId, _ref, index -> tuple(refId, index) }, by: 0)
        .map { _refId, meta, reads, ref, _gtf, _rrna, maxIntron, index ->
            tuple(meta + [genome_size: ref.length, max_intron: maxIntron], reads, index)
        }
    TSS_ALIGN(alignmentInputs)
    sizes = GET_REFERENCES.out.sizes.map { refId, _ref, path -> tuple(refId, path) }
    TSS_CTSS(TSS_ALIGN.out.bam.map { meta, bam, bai -> tuple(meta.ref_id, meta, bam, bai) }
        .combine(sizes, by: 0)
        .map { _refId, meta, bam, bai, path -> tuple(meta, bam, bai, path) })
    TSS_CLUSTERS(TSS_CTSS.out.ctss)
    qcInputs = TSS_CTSS.out.ctss.map { meta, ctss -> tuple(meta.id, meta, ctss) }
        .join(TSS_CLUSTERS.out.clusters.map { meta, clusters -> tuple(meta.id, clusters) },
              failOnDuplicate: true, failOnMismatch: true)
        .join(TSS_TRIM.out.log.map { meta, log -> tuple(meta.id, log) },
              failOnDuplicate: true, failOnMismatch: true)
        .join(TSS_ALIGN.out.log.map { meta, log -> tuple(meta.id, log) },
              failOnDuplicate: true, failOnMismatch: true)
        .map { _id, meta, ctss, clusters, trimLog, starLog ->
            tuple(meta.ref_id, meta, ctss, clusters, trimLog, starLog)
        }
        .combine(annotations, by: 0)
        .combine(sizes, by: 0)
        .map { _refId, meta, ctss, clusters, trimLog, starLog, _ref, gtf, rrna, _maxIntron, path ->
            tuple(meta, ctss, clusters, gtf, rrna, trimLog, starLog, path)
        }
    TSS_QC(qcInputs)
    // Collect channels, not publishDir contents: wait for every library's reports.
    multiqcFiles = FASTQC.out.reports.map { _meta, _mate, _html, zip -> zip }
        .mix(TSS_EXTRACT.out.log.map { _meta, log -> log })
        .mix(TSS_TRIM.out.log.map { _meta, log -> log })
        .mix(TSS_ALIGN.out.log.map { _meta, log -> log })
        .mix(TSS_QC.out.summary.map { _meta, summary -> summary })
        .flatten().collect().map { paths -> paths.sort { it -> it.name } }
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
