include { CHROMAP_INDEX } from './chromap_index'
// The same record check as the STAR index: the header lines, then exactly the index's files.
include { current } from '../rna-seq/get_star_index'

// The chromap release environments/GENESIS_CHROMAP.yaml pins; its image is built from the same file.
def pinnedChromapVersion() {
    def pin = file("${moduleDir}/../../environments/GENESIS_CHROMAP.yaml").readLines()
        .collect { line -> (line =~ /chromap=([^\s,<>]+)/) }
        .find { match -> match.find() }
    if (pin == null) {
        error('environments/GENESIS_CHROMAP.yaml does not pin a chromap version')
    }
    return pin.group(1)
}

// The lines a chromap index record starts with, before the index files' sizes. CHROMAP_INDEX
// writes the same lines (see chromap_index.nf).
def recordHeader(ref, version) {
    return [[ref.sha256, ref.filename], [version, 'chromap']]
}

// A chromap index for each reference, reused from the references folder only while
// <ref_id>.chromap.sha256 matches the compressed FASTA, the pinned chromap release and every
// index file's size; otherwise it is rebuilt there.
workflow GET_CHROMAP_INDEX {
    take:
    references     // tuple(ref, fasta): ref carries id, filename, sha256 and length
    referencesDir  // folder holding the FASTAs and their derived products

    main:
    def version = pinnedChromapVersion()
    checked = references.map { ref, fasta ->
        def name = "${ref.id}.chromap".toString()
        def reusable = current(referencesDir.resolve("${name}.sha256"), recordHeader(ref, version),
                               referencesDir.resolve(name))
        tuple(ref, fasta, reusable)
    }
    branches = checked.branch { _ref, _fasta, reusable ->
        cached: reusable
        missing: true
    }
    CHROMAP_INDEX(branches.missing.map { ref, fasta, _reusable -> tuple(ref, fasta) })

    emit:
    index = branches.cached
        .map { ref, _fasta, _reusable -> tuple(ref.id, referencesDir.resolve("${ref.id}.chromap")) }
        .mix(CHROMAP_INDEX.out.index.map { ref, path -> tuple(ref.id, path) })
}
