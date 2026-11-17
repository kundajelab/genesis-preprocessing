# DAP-seq reference guide

Genesis is a Nextflow DSL2 pipeline for single-end and paired-end plant DAP-seq.
It downloads libraries, derives reference indexes and chromosome sizes, aligns
with bwa-mem3, inspects raw reads with FastQC, runs SPP and alignment QC, generates deepTools coverage tracks,
calls MACS3 peaks against each treatment's assigned control, and quantifies peaks.

![Genesis DAP-seq workflow metro map](../../docs/images/genesis_metro_map.svg)

[Diagram source and interpretation](../../docs/diagrams/README.md)

The curated DAP sheets, one per species (`curated-sample-sheets/<Species>/dap-seq-curated.tsv`),
are in the
[Stanford-Genesis Hugging Face bucket](https://huggingface.co/buckets/Stanford-Genesis/genesis-curation).
[`samplesheets/dap-seq-test.tsv`](../../samplesheets/dap-seq-test.tsv) has 24 of their libraries
for testing: 12 Sorghum and A. thaliana treatments from 10 studies, each with its control, SE
and PE. Use one sheet per analysis; [samplesheets/README.md](../../samplesheets/README.md)
describes the references each needs.

## Install and run

Install [Pixi](https://pixi.sh/), then install the project and its Python tools:

```bash
pixi run install-all
```

Place each sheet's gzipped reference FASTA under `references/`, using the exact
basename in `reference_fasta`. Reference genomes are not downloaded automatically.
The `GET_REFERENCES` subworkflow (`modules/common/get_references.nf`) builds each
reference's chromosome sizes and bwa-mem3 index (in `<id>.bwa-mem2/`) and publishes them
beside the FASTA, with `<fasta>.sha256`. Its first line is the SHA-256 of the compressed FASTA (in
`sha256sum` format); each further line records a product's size in bytes. Every run,
including resume, rehashes each FASTA. It reuses `<id>.chrom.sizes` and
`<id>.bwa-mem2/` only when the hash matches and every product is present at its
recorded size; otherwise it rebuilds both and saves a new record. Replacing or
re-compressing a FASTA, or a replaced or partly copied product, therefore rebuilds
automatically. Stub runs (`-stub-run`) never publish into the references folder. Distinct FASTA basenames
must not collapse to the same reference ID after removing their suffix. The first
resume after generating references can rerun downstream tasks as inputs move to
published paths; subsequent resumes reuse computation.

Run exactly one sheet per invocation. On Perlmutter or Sherlock, submit the run as a Slurm job
([running on Perlmutter and Sherlock](clusters.md)); locally, `pixi run pipeline` runs it with
Docker, which suits debugging on a small sheet. Pixi always runs from the checkout folder, so give
sheets elsewhere as absolute paths:

```bash
# Show wrapper options.
pixi run pipeline --help

# Submit a run on Perlmutter or Sherlock; the cluster profile is chosen automatically.
pixi run submit-pipeline --workflow dap-seq -R /path/references /path/dap.tsv

# Continue a previous run using Nextflow's task cache.
pixi run submit-pipeline --workflow dap-seq -R /path/references /path/dap.tsv -resume

# Debug locally with Docker on the test sheet (24 libraries, 70 GB of FASTQ).
pixi run pipeline --workflow dap-seq samplesheets/dap-seq-test.tsv
```

You can also call `scripts/run-pipeline.sh` directly, from any folder; relative paths
then resolve from your current folder. Either way, Nextflow runs in the checkout's
Pixi environment. Choose the profile with `-p`/`--profile`
rather than Nextflow's `-profile`, which the script already sets. Without `-p`,
`scripts/get-default-profile.sh` chooses `nersc` on Perlmutter, `sherlock` on Sherlock and
`local,docker` elsewhere. A profile names an executor and a software environment together:
`nersc` is Slurm with podman-hpc, and `sherlock` is Slurm with Apptainer, on the `normal`
queue (select another with `--queue`, which is a QOS with `nersc`). Locally built images must be published with
approval before a cluster can use them. Cluster work, cache, reference, and output paths must
be visible to compute nodes and writable by the user; see
[running on Perlmutter and Sherlock](clusters.md).

The sheet can be a file or an `http(s)` URL ending in `.tsv`, which Nextflow reads
directly. Each run lives in `WORKSPACE/RUN_NAME/`, where `RUN_NAME` defaults to the sheet's
basename without `.tsv`, for a URL too (override with `-n`/`--run-name`). The default workspace,
from `scripts/get-default-workspace.sh`, is `$SCRATCH/workspace` when `$SCRATCH`
is set and `workspace/` in this repository otherwise; override it with
`-w`/`--workspace`. A run folder contains:

- `work/`: the Nextflow work directory.
- `output/`: published results.
- `trace/`: timestamped execution report, timeline, trace, and DAG files, and
  `samplesheet_TIMESTAMP.tsv`, the sheet the run read.
- `.nextflow.log` and `.nextflow/`: the wrapper launches Nextflow from the run
  folder, so `-resume` uses that run's history.

With `-P`/`--persist-dir DIR`, `output/` and `trace/` go to `DIR/RUN_NAME/` instead, for storage
that is kept when the workspace is purged; the rest stays in the run folder. It defaults to
project storage on Perlmutter and Sherlock, and to none elsewhere.

Use `-R`/`--references /path` to choose the reference location; the wrapper creates it if
missing. It defaults to project storage on Perlmutter and Sherlock (see
[running on Perlmutter and Sherlock](clusters.md)) and to `references/` in this repository
elsewhere. Like the sheet and `-w`, a relative `-R` resolves from your current folder, whereas
relative paths in Nextflow arguments resolve against the run folder. Nextflow retains intermediate FASTQs and BAMs under its
work directory for caching; retain it for `-resume` and clean it only when those
intermediates are no longer needed. Task scripts run with `bash -euxo pipefail`,
so each task's `.command.log` records the commands it ran.

```bash
# Path to the most recent execution report (optional argument: WORKSPACE).
pixi run get-latest-report
# Follow a task's log by the hash Nextflow prints, hiding `set -x` lines.
pixi run watch-log ab/123456
# Totals for a run (hours, disk, utilization) in PERSIST_DIR/RUN_NAME/resources.tsv, with any
# executor, from its trace files and work folder. Only scripts/submit-pipeline.sh runs it
# automatically, when a run finishes; run it yourself after run-pipeline.sh runs.
pixi run summarize-resources PERSIST_DIR/RUN_NAME
```

## Sample sheet

The header and column order are fixed:

```tsv
sample_id	species	read1_url	read2_url	control_sample	reference_fasta
treatment	Arabidopsis_thaliana	https://example.org/treatment.fastq.gz	-	control	TAIR10.fa.gz
control	Arabidopsis_thaliana	https://example.org/control.fastq.gz	-	-	TAIR10.fa.gz
```

Use `-` for an absent mate or control. A mate URL selects paired-end processing.
Sample IDs must be unique within the sheet. Every assigned control must appear
in the same sheet, have `control_sample = -`, and match the treatment's species,
reference, and layout. Controls receive alignment, QC, and tracks; treatments
also receive peaks and quantification. Multiple treatments can share a control.
The whole sheet is validated before any task starts: nf-schema checks it against
`assets/dap-seq.samples.schema.json` (required values, patterns, unique sample IDs), then
`workflows/dap-seq.nf` checks the exact header, the reference FASTAs, reference IDs and control
assignments.

For interspecies DAP, distinguish the TF-source organism from the organism whose DNA
was assayed. The alignment reference must describe that DNA. Preserve the original
sheet labels, and record the separate biological identities in the curation metadata
([genesis-curation](https://github.com/kundajelab/genesis-curation));
do not infer both from a single species label.

Reference names are gzipped FASTA basenames ending in `.fa.gz`, `.fasta.gz`, or
`.fna.gz`. Identifiers and filenames may contain letters, digits, underscores,
periods, and hyphens, with a letter, digit, or underscore first. Read URLs must
be HTTP or HTTPS. Reads are downloaded in batches of `--download_batch_size`
samples (default 1), up to eight batches at once. When a sample's URLs are ENA's
FASTQs of one whole SRA run (`.../RUN/RUN.fastq.gz`, or `RUN_1` and `RUN_2`), the
run's `.sra` file is fetched from NCBI's AWS copy instead, checked against NCBI's
MD5, and converted with `fasterq-dump` into FASTQs whose decompressed bytes are
ENA's, read names included. If that fails, the sample is fetched from its URLs,
with a warning in the task's log. Other URLs are fetched as given; mates with
matching basenames are fine. `--sra_locate_url` sets NCBI's locator.

## Analysis choices

The main workflow uses bwa-mem3 for alignment. With these options its records match bwa-mem2's,
except for a few repeat reads where bwa-mem3 fixes bwa-mem2's inconsistent scores, and for
the `MQ` and `HN` tags it adds; it is about four times as fast.

Main and first-read QC alignments use **full reads**; only SPP's input is cut
(see below). Read length is inferred
from up to the first 100 reads per mate and recorded as metadata; it is not a
trimming instruction. Inspected reads must be nonempty, structurally valid, and
uniform in length across mates. Reads fetched from their URLs are fully
decompressed, which checks gzip integrity; first-read
line counts must be divisible by four, and paired FASTQs must have equal line
counts. Genome size is the sum of chromosome sizes derived from the FASTA.

The main BAM retains primary mapped alignments, including **MAPQ 0**.
`samtools view -F 2308` excludes unmapped, secondary, and supplementary records.
There is no MAPQ threshold, NH weighting, or separate deduplication step in the
alignment. QC aligns full first reads as an unpaired library, used for
first-read statistics and tracks; for a single-end library that is the main
alignment, reused rather than repeated. A third alignment, for SPP only, uses first
reads cut to their first 50 bases (`--spp_read_length`). Each sample's alignment
provenance records the bwa-mem3 version, layout, read length, flag exclusion,
options, and the QC and SPP read policies.

The currently pinned alignment image contains bwa-mem3 **0.16.0**, the newest bioconda build
for all four platforms the environment locks.

BWA defaults are `-K 10000000 -k 19 -c 500 -T 30`. Override them with
`--bwa_batch_size`, `--bwa_seed_length`, `--bwa_max_seed_occurrences`, and
`--bwa_score_threshold`; each must be a positive integer.
Parameters are validated against `assets/dap-seq.params.schema.json` (nf-schema) before any
task runs, and parameters that differ from the defaults are logged at startup.

SPP uses the unmodified, vendored upstream `run_spp.R` revision
[`6984a713aba0218b76bacc63f2fb5087425fd6a3`](https://github.com/kundajelab/phantompeakqualtools/blob/6984a713aba0218b76bacc63f2fb5087425fd6a3/run_spp.R),
with shift range `-s=-0:2:400`, GNU awk, and first reads cut to 50 bases
before alignment, as in ENCODE's cross-correlation QC. With full-length reads as
long as the fragments (e.g. 151-base reads of ~140-base DAP-seq fragments), the
read-length phantom peak hides the fragment peak and SPP cannot estimate a
fragment length. Cutting before alignment, not after, puts the phantom peak at
the cut length, so NSC and RSC keep their usual meaning. If SPP still finds no
fragment peak, QC records a row with `NA` estimates instead of failing.
The provenance and differences from the historical Sherlock script are explained
in [vendor/phantompeakqualtools/README.md](../../vendor/phantompeakqualtools/README.md).
SPP scores can differ from the old workflow because of the aligner, read lengths,
and upstream correlation-baseline changes.

MACS3 uses `BAM` for SE and `BAMPE` for PE, the assigned control, inferred genome
size, and its remaining defaults. In particular, MACS3's default `--keep-dup=1`
limits identical-position/strand SE tags or identical PE fragments during peak
calling. It still reads alignments carrying the BAM duplicate flag; marking a
BAM is not the same as excluding those records. Thus retaining duplicates in the
alignment BAM does not mean peak calling retains every duplicate observation.

Tracks and primary-weighted quantification use the retained BAM, without a
separate duplicate filter; paired read ends count independently in quantification.
Consequently, excluding duplicate records can change coverage and quantification
even when peak calls remain unchanged. These behaviors describe the current
implementation, not a recommendation to remove duplicates from DAP-seq libraries.

deepTools produces total and strand-specific
CPM coverage and strand-specific 5-prime counts at one-base resolution. PE
libraries additionally receive first-read tracks. Tracks have no MAPQ filter.
Strand CPM uses the total alignment count as its shared scaling denominator.

Peak RPM counts overlapping retained primary **read ends**, dividing by all
retained mapped read ends and multiplying by one million; paired ends count
independently. Mean RPKM is the length-weighted average of deepTools' RPKM
bedGraph over each peak, including uncovered bases as zero. Both score files
preserve the narrowPeak rows and append one value. Mean coverage RPKM is a
coverage statistic and differs from read-end RPM divided by peak length in kb.
A single `QUANTIFY` task reads the BAM with pysam and the bedGraph directly,
with no SAM or BED intermediates; overlaps follow `bedtools intersect` (at least
one shared base). The standalone `genesis-tools peak-rpm` command also supports
explicit NH fractional weighting when input alignments have valid NH tags.

## Outputs and resources

Published results live at `WORKSPACE/RUN_NAME/output/<species>/<sample_id>/`:

- FASTQ line counts, metadata, and alignment provenance TSVs.
- Raw-read FastQC HTML/ZIP pairs and version files in `qc/fastqc/`, named
  `<sample_id>.read1_fastqc.*` (SE and PE) and `<sample_id>.read2_fastqc.*` (PE).
- SPP TSV/PDF, samtools flagstat, idxstats, and alignment statistics; PE insert-length distributions.
- CPM and 5-prime count bigWigs, including PE first-read tracks.
- Treatment MACS3 outputs, including gzipped narrowPeak files.
- Treatment `.peaks.RPM.tsv` and `.peaks.mean_RPKM.tsv` scores.

A single run-level aggregate is published at
`WORKSPACE/RUN_NAME/output/multiqc/multiqc_report.html`, with the complete
`multiqc_data/` directory beside it. MultiQC 1.35 consumes FastQC ZIPs, samtools
flagstat, both main/read1 samtools stats, idxstats, and Genesis metadata directly
from completed workflow channels. Full report identifiers are preserved: main and
read1 statistics, and R1/R2 FastQC, never share an entry. The aggregation fails if
any required parser or expected sample is absent, even if MultiQC exits zero.

The Genesis table and header include sample/control associations, species,
reference ID, run name, Git SHA, source-tree status and pipeline version. Read SPP
results in each library's TSV/PDF files. They are not aggregated into MultiQC
because its native phantompeakqualtools parser cannot handle Genesis's valid `NA`
fallback rows. Missing optional SPP/insert-size companion files do not prevent
aggregation; PE insert-size metrics are included in samtools stats.

`multiqc_data/` retains native JSON/parquet, source records, logs and configuration,
plus `genesis_{fastqc,stats,flagstat,idxstats}.json` and `genesis_provenance.json`.
MultiQC's idxstats parser exports mapped counts and reference lengths; original
per-library TSVs also retain unmapped counts. No biological thresholds are imposed.
The committed `multiqc_config.yaml` disables broad name cleaning and keeps
positive-count contigs visible. Three renderer-only fields (gzip timestamp, UUID,
and embedded creation date) are canonicalized for reproducible HTML; native JSON
retains generation metadata. A stable resume caches the complete report/data set.

### Run records

Each run writes two records to its output folder when it starts:

- `genesis-images.lock`: a copy of `conf/genesis-images.lock`, the container image pins
  (by digest) used by container profiles.
- `git-status.txt`: the pipeline checkout's branch and commit, the container
  engine (`docker`, `apptainer`, … or `none` for stub runs without one), and every uncommitted or
  untracked file reported by `git status --porcelain` (`dirty files: none` for a
  clean checkout).

Together with the Nextflow trace and report, these files record how the run was made.

Per-sample `*.alignment.tsv` files record the aligner version and options, and
the MultiQC report's provenance records the run name, Git commit and clean/dirty
status. Use the run's Nextflow traces for measured task runtime and memory.

Full byte/format scans, scientific QC review and eligibility review are separate from
successful execution.

**FASTQs, BAMs, BAM indexes, and bedGraphs are excluded from published results.** Derived reference files
are published separately under the reference directory.

FastQC 0.12.1 reads each raw downloaded mate once on an independent branch,
with one CPU, 1 GB task memory and a 512 MB Java heap. It
neither trims nor changes reads and does not gate alignment on its quality labels.
FastQC WARN/FAIL labels remain in the reports and do **not** reject a dataset:
plant/assay-specific acceptance thresholds are UNSPECIFIED. Execution errors or
missing reports still fail the workflow. Empty/corrupt gzip libraries are rejected
by DOWNLOAD; FastQC itself is not a comprehensive FASTQ validator (for example,
v0.12.1 accepts a short quality string). A successful report is not evidence of
biological suitability.

bwa-mem3 indexing and alignment size their own requests. Indexing reads the
gzipped FASTA directly on four CPUs (more barely help), with 1 GB plus 16 bytes per
reference base (it peaks near 14.5) and one hour plus an hour per Gbp. It writes no
`.0123`, so a record from bwa-mem2, which listed it, makes the index rebuild once.
Alignment uses 16 CPUs, memory from the genome size and CPUs, and time from the size
of the compressed reads. Like every DAP-seq task, both retry only when they run out of
memory or of time, with twice what ran out, or when the system cancels them, with the same
resources. QC uses two CPUs, and memory and time from the size of SPP's BAM. Tracks use
four CPUs (1.7 times as fast as two), 2.0 GB, and time from the genome size and the size
of the BAMs they cover. bamCoverage writes its temporary files in the task's folder,
since a compute node's `/tmp` may be held in memory. It waits for
workers killed for memory until the time limit, so a task that times out after Slurm
logged an out-of-memory kill retries with twice the memory, not the time. Every module
declares its own resources. The `local` profile caps requests at the machine's CPUs and
memory.

The `test` profile caps every task at one CPU, 2 GB, and five minutes; it does
not select a sample subset or replace tools. On clusters, run computation on
allocated compute nodes within your site's resource limits.

## Environments and Docker builds

Pure tool environments are defined in `environments/*.yaml`; their matching
Pixi manifests and locks support container builds. `genesis_tools/` is a typed
Python 3.14.5 uv project that reads BAMs with pysam, built with the dedicated uv Dockerfile.
Each Nextflow module reads its image from `conf/genesis-images.lock` using nf-dotenv.

Build tool images through Pixi:

```bash
pixi run build-dockers --no-push
pixi run build-dockers --no-push -- genesis_alignment genesis_tools
```

Every build targets `linux/amd64` and `linux/arm64`. The build scripts default
to pushing: always pass `--no-push` for local builds. Tags default to the latest
Git tag (`0.1.0`); use `--tag` to override. The scripts update `conf/genesis-images.lock` with image
digests and preserve unrelated entries. Publishing images requires approval.
The Dockerfiles specify pinned base images.

FastQC and MultiQC use separate environments. Docker defaults pin publicly
available Biocontainers images by registry digest: FastQC 0.12.1 and MultiQC 1.35.
They can be pulled on a fresh Linux server without rebuilding or publishing an
image. Use the committed environment locks for local builds.
The FastQC/MultiQC integration has been validated with Linux amd64 Docker.

## Verification

Run the default checks, then use Docker validation to exercise the scientific tools:

```bash
pixi run checks
# Real-tool Docker validation:
pixi run validate-docker
# Retain fixtures, full logs, traces, and intermediate outputs for inspection:
pixi run validate-docker --keep
```

The default suite runs ShellCheck, Ruff, ty, quantification and mocked Docker
build regressions, sample-sheet and metadata cases, and actual Nextflow SE/PE
stub workflows. It checks control fan-out, per-reference generation and reuse,
`-resume`, invalid sheets and references, and publication boundaries. Mock builds cover
both architectures and `--no-push` without contacting Docker or a registry.

Docker validation runs the real download/staging, indexing, bwa-mem3, samtools,
FastQC, SPP, deepTools, MACS3, and Python tools on deterministic 75-base SE/PE
fixtures with enriched regions and duplicated sequence. It checks nonempty peaks,
finite positive scores, provenance, read lengths (full, and 50 bases for SPP),
outputs, and resume.
FastQC coverage requires seven mate tasks, seven HTML/ZIP pairs, report contents,
unchanged raw-read bytes, publication and cached reports on resume. Additional
real-module cases cover advisory quality FAIL labels, empty input, a known
short-quality parser limitation, malformed records and corrupt gzip. Run these
focused cases independently with:

```bash
pixi run uv run --frozen --project genesis_tools python tests/verify_fastqc.py --keep
```

The five-library fixture expects 41 initial tasks and 37 tasks after reference
reuse, including exactly one MULTIQC aggregation.
Task counts depend on the workflow revision; inspect the trace for completion and
cache status. Tool warnings are retained for inspection. MultiQC regressions require seven
FastQC, ten samtools stats, five flagstat and five idxstats identities, with
numerical read-count checks and report sections. They cover shared controls,
optional missing SPP, duplicate identities, collision-prone names, missing required
modules despite a zero tool exit, fresh deterministic HTML/metrics, and cached
full output directories.

These small synthetic Docker runs test integration, not biological performance
on full public datasets.
Native macOS stub traces do not include task runtime metrics; Docker traces do.

See [docs/nextflow-style.md](../../docs/nextflow-style.md) for module conventions.
