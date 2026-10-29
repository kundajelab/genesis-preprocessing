// Runs the raw-read steps of workflows/scatac-reads.nf after the SRA download, on local FASTQs, so
// tests/verify_scatac_reads.py can check chromap and cell calling with real tools offline.
nextflow.enable.dsl = 2

include { SCATAC_CHROMAP } from '../modules/scatac-seq/scatac_chromap'
include { SCATAC_CELLS } from '../modules/scatac-seq/scatac_cells'

workflow {
    def inputs = file(params.inputs)
    def libraries = params.libraries.tokenize(',').collect { spec ->
        // id:barcode_length:whitelist (whitelist '-' for none)
        def (id, length, whitelist) = spec.tokenize(':')
        def meta = [id: id, species: 'Plant_synthetic', barcode_length: length.toInteger(),
                    barcode_whitelist: whitelist, cells_url: '-', organelles: ['ChrM'],
                    genome_length: 1020000L]
        tuple(meta, inputs.resolve("${id}.R1.fastq.gz"), inputs.resolve("${id}.R2.fastq.gz"),
              inputs.resolve("${id}.barcode.fastq.gz"), file(params.fasta),
              whitelist == '-' ? [] : inputs.resolve(whitelist))
    }
    SCATAC_CHROMAP(channel.fromList(libraries))
    SCATAC_CELLS(SCATAC_CHROMAP.out.fragments)
}
