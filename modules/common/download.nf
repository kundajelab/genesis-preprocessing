include { DOWNLOAD_SRA } from './download_sra'
include { DOWNLOAD_URLS } from './download_urls'
include { sraRun } from './download_support'

// One sample's downloaded mates, from a batch's FASTQs, in mate order; each sample must have
// exactly its layout's mates.
def matesOf(meta, reads) {
    def names = ['read1', 'read2'].collect { mate -> "${meta.id}.${mate}.fastq.gz".toString() }
    def mates = reads.findAll { read -> read.name in names }.sort { read -> read.name }
    if (mates.size() != (meta.layout == 'PE' ? 2 : 1)) {
        error("${meta.id}: downloaded ${mates*.name}, expected its ${meta.layout} mates")
    }
    return mates
}

// Each sample's raw reads, in batches of --download_batch_size samples per task: from its SRA run
// when its URLs are ENA's FASTQs of one whole run (DOWNLOAD_SRA), otherwise, or when that run is
// unavailable, from its URLs (DOWNLOAD_URLS). Each fallback is logged with its reason.
workflow DOWNLOAD {
    take:
    samples  // one sample's meta map per item

    main:
    def batchSize = params.download_batch_size.toInteger()
    def sources = samples.branch { meta ->
        sra: sraRun(meta) != null
        urls: true
    }
    DOWNLOAD_SRA(sources.sra.collate(batchSize))
    def outcomes = DOWNLOAD_SRA.out.downloads
        .flatMap { batch, status, reads ->
            def recorded = status.readLines().collectEntries { line ->
                def fields = line.split('\t', -1) as List
                [(fields[0]): fields.drop(1)]
            }
            batch.collect { meta ->
                def outcome = recorded[meta.id]
                if (outcome == null || !(outcome[0] in ['ok', 'fallback'])) {
                    error("${meta.id}: DOWNLOAD_SRA recorded no outcome")
                }
                outcome[0] == 'ok' ? tuple(meta, matesOf(meta, reads), null) : tuple(meta, null, outcome[1])
            }
        }
        .branch { _meta, _mates, reason ->
            converted: reason == null
            fallback: true
        }
    def fallbacks = outcomes.fallback.map { meta, _mates, reason ->
        log.warn("${meta.id}: SRA run ${sraRun(meta)} unavailable (${reason}); downloading its FASTQs from their URLs")
        meta
    }
    // Fallbacks are batched apart, so the sheet's own URL samples batch the same on every run.
    DOWNLOAD_URLS(sources.urls.collate(batchSize).mix(fallbacks.collate(batchSize)))
    def fetched = DOWNLOAD_URLS.out.downloads.flatMap { batch, reads ->
        def files = reads instanceof List ? reads : [reads]
        batch.collect { meta -> tuple(meta, matesOf(meta, files)) }
    }

    emit:
    // tuple(meta, mates): the sample's FASTQs, read1 then read2 for paired-end.
    downloads = outcomes.converted.map { meta, mates, _reason -> tuple(meta, mates) }.mix(fetched)
}
