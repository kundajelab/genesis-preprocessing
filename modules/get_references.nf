include { CHROM_SIZES } from './chrom_sizes'
include { BWA_MEM2_INDEX } from './bwa_mem2_index'

// SHA-256 of a file's bytes, streamed rather than read into memory.
def sha256(path) {
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    new java.security.DigestInputStream(path.newInputStream(), digest).withCloseable { stream ->
        stream.transferTo(java.io.OutputStream.nullOutputStream())
    }
    return digest.digest().encodeHex().toString()
}

// Whether <fasta>.sha256 matches this FASTA hash and records exactly these products, each
// present with its recorded size. Anything else, including a malformed record, means rebuild.
def current(saved, hash, referencesDir, products) {
    if (!saved.exists()) {
        return false
    }
    def lines = saved.readLines().findAll { line -> line.trim() }
    def recorded = lines.drop(1).collect { line -> line.tokenize() }
    if (!lines || lines[0].tokenize()[0] != hash ||
            recorded.any { fields -> fields.size() != 2 || !fields[0].isLong() }) {
        return false
    }
    def sizes = recorded.collectEntries { fields -> [(fields[1]): fields[0].toLong()] }
    return sizes.size() == recorded.size() && sizes.keySet() == products.toSet() &&
        products.every { relative ->
            def path = referencesDir.resolve(relative)
            path.exists() && path.size() > 0 && path.size() == sizes[relative]
        }
}

// Chromosome sizes and a bwa-mem2 index for each reference, reused from the references
// folder only while <fasta>.sha256 matches the compressed FASTA and every product's recorded
// size; otherwise both are rebuilt there and a new record is saved.
workflow GET_REFERENCES {
    take:
    references     // tuple(ref, fasta): one per reference, ref = [id, filename]
    referencesDir  // folder holding the FASTAs and their derived products

    main:
    // Rehashed on every run, including -resume, so a replaced FASTA is always noticed.
    checked = references.map { ref, fasta ->
        def hash = sha256(fasta)
        def products = ["${ref.id}.chrom.sizes".toString()] +
            ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac'].collect { suffix ->
                "${ref.id}.bwa-mem2/genome${suffix}".toString()
            }
        def reusable = current(referencesDir.resolve("${ref.filename}.sha256"), hash, referencesDir, products)
        tuple(ref + [sha256: hash], fasta, reusable)
    }
    branches = checked.branch { _ref, _fasta, current ->
        cached: current
        missing: true
    }
    missing = branches.missing.map { ref, fasta, _current -> tuple(ref, fasta) }

    CHROM_SIZES(missing)
    built_sizes = CHROM_SIZES.out.sizes.map { ref, path -> tuple(ref.id, ref, path) }
    // Indexing sizes its memory and time from the reference length.
    index_inputs = missing.map { ref, fasta -> tuple(ref.id, ref, fasta) }
        .join(built_sizes, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _id, ref, fasta, _sizesRef, path ->
            def length = path.readLines().sum { line -> line.tokenize('\t')[1].toLong() }
            tuple(ref + [length: length], fasta, path)
        }
    BWA_MEM2_INDEX(index_inputs)

    cached = branches.cached.map { ref, _fasta, _current -> ref }

    emit:
    sizes = cached.map { ref -> tuple(ref.id, ref, referencesDir.resolve("${ref.id}.chrom.sizes")) }
        .mix(built_sizes)
    indexes = cached.map { ref -> tuple(ref.id, ref, referencesDir.resolve("${ref.id}.bwa-mem2")) }
        .mix(BWA_MEM2_INDEX.out.index.map { ref, path -> tuple(ref.id, ref, path) })
}
