// Runs DOWNLOAD_SCATAC_READS on one run, so tests/verify_scatac_reads.py can check its fetch and
// read selection offline: a local server stands in for NCBI's locator and AWS, and stand-in
// fastq-dump, vdb-dump and pigz go first on PATH.
nextflow.enable.dsl = 2

include { DOWNLOAD_SCATAC_READS } from '../modules/scatac-seq/download_scatac_reads'

workflow {
    def meta = [id: 'lib', barcode_length: params.barcode_length as int]
    DOWNLOAD_SCATAC_READS(channel.of(tuple(meta, params.run)))
}
