// Runs SRA_FASTQ's read selection on one run, so tests/verify_scatac_reads.py can check it offline
// with stand-in fasterq-dump and pigz on PATH.
nextflow.enable.dsl = 2

include { SRA_FASTQ } from '../modules/scatac-seq/sra_fastq'

workflow {
    def meta = [id: 'lib', barcode_length: params.barcode_length as int]
    SRA_FASTQ(channel.of(tuple(meta, params.run, file(params.sra))))
}
