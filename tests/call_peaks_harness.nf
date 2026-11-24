// Runs CALL_PEAKS on the samples in --input (sample_id, layout, treatment, control: paths), so
// tests/verify_call_peaks.py can check it offline with a stand-in macs3 on PATH. The samples that
// get peaks, the input QUANTIFY takes in workflows/dap-seq.nf, are listed in --run_dir/peaks.txt.
nextflow.enable.dsl = 2

include { CALL_PEAKS } from '../modules/dap-seq/call_peaks'

workflow {
    def samples = channel.fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [id: row.sample_id, species: 'Plant', layout: row.layout, genome_size: 1000]
            tuple(meta, file(row.treatment), file(row.treatment), file(row.control),
                  file(row.control))
        }
    CALL_PEAKS(samples)
    CALL_PEAKS.out.peaks
        .map { meta, _peaks -> meta.id }
        .collectFile(name: 'peaks.txt', storeDir: params.run_dir, newLine: true, sort: true)
}
