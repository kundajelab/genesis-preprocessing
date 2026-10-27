include { RNA_INDEX } from './rna_index'
include { sha256 } from '../common/get_references'

def log2(value) {
    return Math.log(value as double) / Math.log(2)
}

// STAR index parameters for this genome, from the STAR manual: small genomes and assemblies with
// many scaffolds need smaller suffix-array and bin settings. Computed once here so RNA_INDEX and the
// reuse check agree. ref.length and ref.contigs come from the chromosome sizes.
def starParameters(ref) {
    return [
        sjdbOverhang: 100,
        genomeSAindexNbases: Math.min(14, Math.floor(log2(ref.length) / 2 - 1) as int),
        genomeChrBinNbits: Math.min(18,
            Math.floor(log2(Math.max(ref.length / ref.contigs as double, 100d))) as int),
    ]
}

// The STAR version environments/GENESIS_STAR.yaml pins; its image is built from the same file.
def pinnedStarVersion() {
    def pin = file("${moduleDir}/../../environments/GENESIS_STAR.yaml").readLines()
        .collect { line -> (line =~ /star=([^\s,<>]+)/) }
        .find { match -> match.find() }
    if (pin == null) {
        error('environments/GENESIS_STAR.yaml does not pin a STAR version')
    }
    return pin.group(1)
}

// The lines a STAR index record starts with, before the index files' sizes. RNA_INDEX writes the
// same lines (see rna_index.nf); "STAR" is filled in from STAR --version there.
def recordHeader(ref, version) {
    return [
        [ref.sha256, ref.filename],
        [ref.gtf_sha256, "${ref.id}.genes.gtf".toString()],
        [version, 'STAR'],
    ] + ref.star.collect { name, value -> [value.toString(), name] }
}

// Whether <index>.sha256 has this header and records exactly the files in the index folder, each
// at its recorded size. Anything else, including a malformed record, means rebuild.
def current(saved, header, index) {
    if (!saved.exists() || !index.isDirectory()) {
        return false
    }
    def lines = saved.readLines().findAll { line -> line.trim() }.collect { line -> line.tokenize() }
    def recorded = lines.drop(header.size())
    if (lines.take(header.size()) != header ||
            recorded.any { fields -> fields.size() != 2 || !fields[0].isLong() }) {
        return false
    }
    def sizes = recorded.collectEntries { fields -> [(fields[1]): fields[0].toLong()] }
    def files = index.listFiles().findAll { path -> path.isFile() }
        .collectEntries { path -> [("${index.name}/${path.name}".toString()): path.size()] }
    return sizes.size() == recorded.size() && files == sizes
}

// A STAR index for each reference and annotation, reused from the references folder only while
// <ref_id>.<annotation>.star.sha256 matches the FASTA, the generated GTF, the pinned STAR version,
// the index parameters and every index file's size; otherwise it is rebuilt there.
workflow GET_STAR_INDEX {
    take:
    references     // tuple(ref, fasta, gtf): ref carries sha256, length, contigs and annotation
    referencesDir  // folder holding the FASTAs and their derived products

    main:
    def version = pinnedStarVersion()
    // Rehashed on every run, including -resume, so a changed annotation is always noticed.
    checked = references.map { ref, fasta, gtf ->
        // The annotation's own name, without .gz, keeps indexes for different annotations apart.
        def name = "${ref.id}.${ref.annotation.replaceFirst(/\.gz$/, '')}.star".toString()
        def indexed = ref + [gtf_sha256: sha256(gtf), star: starParameters(ref), star_index: name]
        def reusable = current(referencesDir.resolve("${name}.sha256"), recordHeader(indexed, version),
                               referencesDir.resolve(name))
        tuple(indexed, fasta, gtf, reusable)
    }
    branches = checked.branch { _ref, _fasta, _gtf, current ->
        cached: current
        missing: true
    }
    RNA_INDEX(branches.missing.map { ref, fasta, gtf, _current -> tuple(ref, fasta, gtf) })

    emit:
    index = branches.cached
        .map { ref, _fasta, _gtf, _current -> tuple(ref.id, ref, referencesDir.resolve(ref.star_index)) }
        .mix(RNA_INDEX.out.index.map { ref, path -> tuple(ref.id, ref, path) })
}
