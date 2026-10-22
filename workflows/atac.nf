nextflow.enable.dsl = 2

include { ATAC_ACQUIRE; ATAC_REFERENCE; ATAC_FASTQC; ATAC_ADAPTERS; ATAC_ALIGN; ATAC_FRAGMENTS; ATAC_PEAKS; ATAC_ENRICHMENT; ATAC_REPORT; ATAC_TRACKS } from '../modules/atac/processes'

workflow {
    def plan = new groovy.json.JsonSlurperClassic().parseText(file(params.manifest).text)
    def code = file(params.source)
    def libraries = Channel.fromList(plan.libraries)
    def manifests = libraries.map { meta -> tuple(meta, file(meta.manifest), meta.lanes.collectMany { lane -> [lane.read1, lane.read2] }.findAll { !it.startsWith('https://') }.collect { file(it) }) }
    def referenceInputs = libraries.map { meta -> meta.reference }.unique { ref -> ref.reference_id }
        .map { ref -> tuple(ref, file(ref.fasta)) }
    ATAC_REFERENCE(referenceInputs)
    ATAC_ACQUIRE(manifests, code, params.source_digest)
    ATAC_FASTQC(ATAC_ACQUIRE.out.reads.flatMap { meta, r1, r2 -> [tuple(meta, 'R1', r1), tuple(meta, 'R2', r2)] })
    ATAC_ADAPTERS(ATAC_ACQUIRE.out.reads)
    def alignmentInputs = ATAC_ADAPTERS.out.reads.map { meta, r1, r2 -> tuple(meta.reference_id, meta, r1, r2) }
        .combine(ATAC_REFERENCE.out.index, by: 0).map { ref, meta, r1, r2, index -> tuple(meta, r1, r2, index) }
    ATAC_ALIGN(alignmentInputs)
    ATAC_FRAGMENTS(ATAC_ALIGN.out.bam.map { meta, bam -> tuple(meta, bam, file(meta.manifest)) }, code, params.source_digest)
    ATAC_PEAKS(ATAC_FRAGMENTS.out.bam)
    ATAC_TRACKS(ATAC_FRAGMENTS.out.cuts)
    def enrichmentInputs = ATAC_FRAGMENTS.out.bam.map { meta, bam -> tuple(meta.library_id, meta, bam) }
        .join(ATAC_FRAGMENTS.out.fragments.map { meta, bed, tbi, sizes -> tuple(meta.library_id, bed, sizes) })
        .join(ATAC_PEAKS.out.peaks.map { meta, peaks -> tuple(meta.library_id, peaks) })
        .map { id, meta, bam, bed, sizes, peaks -> tuple(meta, bam, bed, sizes, peaks, file(meta.reference.tss)) }
    ATAC_ENRICHMENT(enrichmentInputs, code, params.source_digest)
    def reports = ATAC_FASTQC.out.reports.map { meta, zip, html -> tuple(meta.library_id, zip) }
        .mix(ATAC_ALIGN.out.qc.flatMap { meta, files -> files.collect { f -> tuple(meta.library_id, f) } })
        .mix(ATAC_FRAGMENTS.out.qc.flatMap { meta, files -> files.collect { f -> tuple(meta.library_id, f) } })
        .mix(ATAC_ENRICHMENT.out.metrics.map { meta, file -> tuple(meta.library_id, file) })
        .groupTuple().join(libraries.map { meta -> tuple(meta.library_id, meta) })
        .map { id, files, meta -> tuple(meta, files.sort { it.name }) }
    ATAC_REPORT(reports, code, params.source_digest)
}
