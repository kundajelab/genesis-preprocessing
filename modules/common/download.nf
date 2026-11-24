include { DOWNLOAD_SRA } from './download_sra'
include { DOWNLOAD_URLS } from './download_urls'
include { sraRun } from './download_support'

// A sample's downloaded mates, read1 then read2; exactly its layout's mates.
def checkedMates(meta, mates) {
    if (mates.size() != (meta.layout == 'PE' ? 2 : 1)) {
        error("${meta.id}: downloaded ${mates*.name}, expected its ${meta.layout} mates")
    }
    return mates.sort { read -> read.name }
}

// Each sample's raw reads, one sample per task: from its SRA run when its URLs are ENA's FASTQs
// of one whole run (DOWNLOAD_SRA), otherwise, or when that run is unavailable, from its URLs
// (DOWNLOAD_URLS). Each fallback is logged with its reason.
workflow DOWNLOAD {
    take:
    samples  // one sample's meta map per item

    main:
    def sources = samples.branch { meta ->
        sra: sraRun(meta) != null
        urls: true
    }
    DOWNLOAD_SRA(sources.sra)
    def outcomes = DOWNLOAD_SRA.out.downloads.branch { meta, status, _reason, _mates ->
        converted: status == 'ok'
        fallback: status == 'fallback'
        other: true
    }
    outcomes.other.subscribe { meta, status, _reason, _mates ->
        error("${meta.id}: DOWNLOAD_SRA recorded an unknown status '${status}'")
    }
    def fallbacks = outcomes.fallback.map { meta, _status, reason, _mates ->
        log.warn("${meta.id}: SRA run ${sraRun(meta)} unavailable (${reason}); downloading its FASTQs from their URLs")
        meta
    }
    DOWNLOAD_URLS(sources.urls.mix(fallbacks))

    emit:
    // tuple(meta, mates): the sample's FASTQs, read1 then read2 for paired-end.
    downloads = outcomes.converted.map { meta, _status, _reason, mates -> tuple(meta, checkedMates(meta, mates)) }
        .mix(DOWNLOAD_URLS.out.downloads.map { meta, mates -> tuple(meta, checkedMates(meta, mates)) })
}
