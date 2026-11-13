// Runs DOWNLOAD on the samples in --input (sample_id, layout, read1_url, read2_url), all in one
// batch, so tests/verify_download.py can check its SRA and URL paths offline with stand-in
// fasterq-dump and pigz on PATH.
nextflow.enable.dsl = 2

include { DOWNLOAD } from '../modules/common/download'

workflow {
    def samples = channel.fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [id: row.sample_id, species: 'Plant', layout: row.layout, read1_url: row.read1_url, read2_url: row.read2_url] }
        .toList()
    DOWNLOAD(samples)
}
