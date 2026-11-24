// Functions shared by DOWNLOAD (download.nf), its processes DOWNLOAD_SRA and DOWNLOAD_URLS, and
// DOWNLOAD_SCATAC_READS (modules/scatac-seq/download_scatac_reads.nf).

// The size in bytes the server reports for a URL (HTTP HEAD), or -1 if it reports none.
def contentLength(url) {
    try {
        def connection = URI.create(url).toURL().openConnection() as java.net.HttpURLConnection
        connection.requestMethod = 'HEAD'
        connection.connectTimeout = 10000
        connection.readTimeout = 10000
        def length = connection.responseCode == 200 ? connection.contentLengthLong : -1L
        connection.disconnect()
        return length
    } catch (Exception _error) {
        return -1L
    }
}

// A sample's read URLs: one for single-end, two for paired-end.
def readUrls(meta) {
    meta.layout == 'PE' ? [meta.read1_url, meta.read2_url] : [meta.read1_url]
}

// The summed size of a sample's FASTQs as their servers report them, or -1 if any reports none.
// Looked up here, on the manager, rather than passed in, so sizes change no task's cache key.
def fastqBytes(meta) {
    def sizes = readUrls(meta).collect { url -> contentLength(url) }
    return sizes.any { size -> size < 0 } ? -1L : sizes.sum() as long
}

// Time to fetch and convert a sample's SRA run: 15 minutes plus 3 per GB of the run's FASTQs. It
// took 3.5 per GB with 4 CPUs and pigz's default level 6, mostly compressing; 8 CPUs at level 4
// should take about 1. Where a server reports no size, 2 hours. Doubled on each retry after an
// attempt ran out of time (retryTimeOf).
def sraTime(meta) {
    def bytes = fastqBytes(meta)
    def seconds = bytes < 0 ? 2L * 3600 : 900L + ((bytes / 1e9 * 180) as long)
    return (seconds * 1000L).ms
}

// Time to fetch a sample's FASTQs from their URLs: 30 minutes plus the FASTQs at 0.5 MB/s, below
// the slowest rate seen from ENA (0.67 MB/s for a whole file), times the attempt. Where a server
// reports no size, a generous 15 hours.
def urlTime(meta, attempt) {
    def bytes = fastqBytes(meta)
    def seconds = bytes < 0 ? 15L * 3600 : 1800L + ((bytes / 0.5e6) as long)
    return "${seconds * attempt}s"
}

// The SRA run whose reads a sample's URLs are, when they are ENA's FASTQs of one whole run:
// .../RUN/RUN.fastq.gz for single-end, .../RUN/RUN_1.fastq.gz and .../RUN/RUN_2.fastq.gz for
// paired-end. Otherwise null, and the sample is fetched from its URLs as given.
def sraRun(meta) {
    def first = meta.read1_url =~ /\/([SED]RR[0-9]+)\/([SED]RR[0-9]+)(_1)?\.fastq\.gz$/
    if (!first.find() || first.group(1) != first.group(2)) {
        return null
    }
    def run = first.group(1)
    if (meta.layout == 'PE') {
        def second = meta.read2_url.endsWith("/${run}/${run}_2.fastq.gz")
        return first.group(3) && second ? run : null
    }
    return first.group(3) ? null : run
}

// The size in bytes of an SRA run's .sra file, as NCBI's locator (params.sra_locate_url) lists it,
// or -1 if it lists none. Looked up here, on the manager, so sizes change no task's cache key. A
// stub run fetches nothing, so it does not ask.
def sraBytes(run) {
    if (workflow.stubRun) {
        return -1L
    }
    try {
        def connection = URI.create("${params.sra_locate_url}${run}").toURL().openConnection()
        connection.connectTimeout = 10000
        connection.readTimeout = 10000
        def located = connection.inputStream.withCloseable { stream -> stream.getText('UTF-8') }
        def size = located =~ /"type": *"sra"[^}]*?"size": *([0-9]+)/
        return size.find() ? size.group(1) as long : -1L
    } catch (Exception _error) {
        return -1L
    }
}

// Time to fetch and convert an SRA run with its technical reads (DOWNLOAD_SCATAC_READS): 30
// minutes plus 4 per GB of its .sra file. On Perlmutter, 8 CPUs converted a 3.4 GB .sra file in
// 3 minutes; fetching it took 1.5 minutes in the RNA-seq run and 16 at a slow time (3.5 MB/s), so
// most of the time is for a slow transfer. Where the locator gives no size, 4 hours. Doubled on each
// retry after an attempt ran out of time (retryTimeOf).
def sraReadsTime(run) {
    def bytes = sraBytes(run)
    def seconds = bytes < 0 ? 4L * 3600 : 1800L + ((bytes / 1e9 * 240) as long)
    return (seconds * 1000L).ms
}
