// Runs the raw-read steps of workflows/scatac-reads.nf after the SRA download, on local FASTQs, so
// tests/verify_scatac_reads.py can check chromap and cell calling with real tools offline.
nextflow.enable.dsl = 2

include { GET_CHROMAP_INDEX } from '../modules/scatac-seq/get_chromap_index'
include { SCATAC_CHROMAP } from '../modules/scatac-seq/scatac_chromap'
include { SCATAC_CELLS } from '../modules/scatac-seq/scatac_cells'
include { sha256 } from '../modules/common/get_references'

workflow {
    def inputs = file(params.inputs)
    def fasta = file(params.fasta)
    // The index is built in, or reused from, params.references, as in workflows/scatac-reads.nf.
    def ref = [id: fasta.name.replaceFirst(/\.fa\.gz$/, ''), filename: fasta.name,
               sha256: sha256(fasta), length: 1020000L]
    GET_CHROMAP_INDEX(channel.of(tuple(ref, fasta)), file(params.references))
    def libraries = params.libraries.tokenize(',').collect { spec ->
        // id:barcode_length:whitelist (whitelist '-' for none)
        def (id, length, whitelist) = spec.tokenize(':')
        def meta = [id: id, ref_id: ref.id, species: 'Plant_synthetic',
                    barcode_length: length.toInteger(), barcode_whitelist: whitelist, cells_url: '-',
                    organelles: ['ChrM'], genome_length: 1020000L]
        tuple(ref.id, meta, inputs.resolve("${id}.R1.fastq.gz"), inputs.resolve("${id}.R2.fastq.gz"),
              inputs.resolve("${id}.barcode.fastq.gz"),
              whitelist == '-' ? [] : inputs.resolve(whitelist))
    }
    SCATAC_CHROMAP(channel.fromList(libraries)
        .combine(GET_CHROMAP_INDEX.out.index, by: 0)
        .map { _refId, meta, r1, r2, barcodes, whitelist, index ->
            tuple(meta, r1, r2, barcodes, fasta, index, whitelist)
        })
    SCATAC_CELLS(SCATAC_CHROMAP.out.fragments)
}
