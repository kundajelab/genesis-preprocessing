nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { VALIDATE_SHEET } from './modules/validate_sheet'
include { CHROM_SIZES } from './modules/chrom_sizes'
include { BWA_MEM2_INDEX } from './modules/bwa_mem2_index'
include { DOWNLOAD } from './modules/download'
include { METADATA } from './modules/metadata'
include { BWA_MEM2_ALIGN } from './modules/bwa_mem2_align'
include { QC } from './modules/qc'
include { TRACKS } from './modules/tracks'
include { CALL_PEAKS } from './modules/call_peaks'
include { QUANTIFY } from './modules/quantify'

// Bowtie 1 modules deliberately have no imports or calls here.
workflow {
    // nextflow_schema.json checks the input sheet, paths, resources, and BWA options.
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
    VALIDATE_SHEET(Channel.value(sheet), Channel.value(referencesDir),
                   file("${projectDir}/genesis_tools/src"))
    samples = VALIDATE_SHEET.out.sheet.splitCsv(header: true, sep: '\t')
        .map { row ->
            [id: row.sample_id, species: row.species, read1_url: row.read1_url,
             read2_url: row.read2_url, control: row.control_sample,
             reference_fasta: row.reference_fasta,
             ref_id: row.reference_fasta.replaceFirst(/\.(fa|fasta|fna)\.gz$/, ''),
             layout: row.read2_url == '-' ? 'SE' : 'PE']
        }
    referenceInputs = samples.unique { sample -> sample.ref_id }
        .map { meta ->
            def ref = [id: meta.ref_id, filename: meta.reference_fasta]
            tuple(ref, file("${referencesDir}/${ref.filename}", checkIfExists: true))
        }
    sizesBranches = referenceInputs.branch { ref, _fasta ->
        cached: file("${referencesDir}/${ref.id}.chrom.sizes").exists() &&
                file("${referencesDir}/${ref.id}.chrom.sizes").size() > 0
        missing: true
    }
    CHROM_SIZES(sizesBranches.missing)
    sizes_ch = sizesBranches.cached.map { ref, _fasta ->
            tuple(ref.id, ref, file("${referencesDir}/${ref.id}.chrom.sizes"))
        }
        .mix(CHROM_SIZES.out.sizes.map { ref, path -> tuple(ref.id, ref, path) })
    indexBranches = referenceInputs.branch { ref, _fasta ->
        cached: ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac'].every { suffix ->
            def path = file("${referencesDir}/${ref.id}.bwa-mem2/genome${suffix}")
            path.exists() && path.size() > 0
        }
        missing: true
    }
    // Indexing sizes its memory and time from the reference length.
    genomeLengths = sizes_ch.map { id, _ref, path ->
        tuple(id, path.readLines().sum { line -> line.tokenize('\t')[1].toLong() })
    }
    indexInputs = indexBranches.missing.map { ref, fasta -> tuple(ref.id, ref, fasta) }
        .join(genomeLengths, by: 0, failOnDuplicate: true)
        .map { _id, ref, fasta, length -> tuple(ref + [length: length], fasta) }
    BWA_MEM2_INDEX(indexInputs)
    indexes = indexBranches.cached.map { ref, _fasta ->
            tuple(ref.id, ref, file("${referencesDir}/${ref.id}.bwa-mem2"))
        }
        .mix(BWA_MEM2_INDEX.out.index.map { ref, path -> tuple(ref.id, ref, path) })
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
