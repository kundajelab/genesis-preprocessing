// Runs DOWNLOAD on the samples in --input (sample_id, layout, read1_url, read2_url), so
// tests/verify_download.py can check its SRA and URL paths offline with stand-in fasterq-dump and
// pigz on PATH. Each sample's downloaded mates are listed in --run_dir/downloads.tsv: ID, then
// the mates' file names.
nextflow.enable.dsl = 2

include { DOWNLOAD } from '../modules/common/download'

workflow {
    def samples = channel.fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [id: row.sample_id, species: 'Plant', layout: row.layout, read1_url: row.read1_url, read2_url: row.read2_url] }
    DOWNLOAD(samples)
    DOWNLOAD.out.downloads
        .map { meta, mates -> ([meta.id] + mates*.name).join('\t') }
        .collectFile(name: 'downloads.tsv', storeDir: params.run_dir, newLine: true, sort: true)
}
