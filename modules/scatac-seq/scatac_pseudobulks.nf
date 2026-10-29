include { SCATAC_FRAGMENTS } from './scatac_fragments'
include { SCATAC_PEAKS } from './scatac_peaks'
include { SCATAC_TRACKS } from './scatac_tracks'
include { SCATAC_QC } from './scatac_qc'

// Everything after a library's fragments and cell list exist, shared by workflows/scatac-seq.nf
// (published fragments) and workflows/scatac-reads.nf (fragments made from raw 10x reads): keep
// listed cells on nuclear contigs, then peaks, tracks and QC for each cell group.
workflow SCATAC_PSEUDOBULKS {
    take:
    libraries      // tuple(meta, fragments, cells): gzipped 10x-style fragments and barcode<TAB>group
    sizes          // tuple(ref_id, chrom.sizes)
    referencesDir  // folder holding the gene annotations
    resources      // [cpus, memory, time] limits for the largest tasks

    main:
    // MACS3's genome size counts nuclear contigs only, as in bulk ATAC.
    sized = libraries
        .map { meta, fragments, cells -> tuple(meta.ref_id, meta, fragments, cells) }
        .combine(sizes, by: 0)
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
        .combine(sizes, by: 0)
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

    emit:
    // MultiQC inputs: each library's summary and each pseudobulk's QC and TSS profile.
    multiqc = SCATAC_FRAGMENTS.out.summary.map { _meta, summary -> summary }
        .mix(SCATAC_QC.out.summary.map { _meta, summary -> summary })
        .mix(SCATAC_QC.out.tss.map { _meta, files -> files })
}
