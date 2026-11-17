// The run's sample sheet, params.input: a path, or an http(s) URL that Nextflow reads directly.
// Its contents are first saved beside the run's trace files, as
// PUBLISH_BASE/trace/samplesheet_TIMESTAMP.tsv with the trace files' timestamp, and that copy is
// returned: the workflow checks and runs exactly the rows the run records, and a remote sheet is
// fetched once, even if it later changes.
def savedSampleSheet() {
    def copy = file("${params.publish_base}/trace/samplesheet_${params.trace_timestamp}.tsv")
    copy.parent.mkdirs()
    copy.bytes = file(params.input).bytes
    return copy
}
