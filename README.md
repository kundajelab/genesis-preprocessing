# Genesis plant DAP-seq pipeline

Genesis is a Nextflow DSL2 pipeline for single-end and paired-end plant DAP-seq.
It downloads libraries, derives reference indexes and chromosome sizes, aligns
with bwa-mem2, runs SPP and alignment QC, generates deepTools coverage tracks,
calls MACS3 peaks against each treatment's assigned control, and quantifies peaks.

| Sample sheet | Layout | Libraries | Treatments |
| --- | --- | ---: | ---: |
| `01-Arabidopsis_thaliana-GSE60141.tsv` | SE | 936 | 934 |
| `15-Arabidopsis_lyrata-PRJNA1177479.tsv` | PE | 405 | 378 |
| `16-Arabidopsis_thaliana-PRJNA1177481.tsv` | PE | 800 | 748 |
| `24-Sorghum_bicolor-PRJNA1177471.tsv` | PE | 142 | 134 |

The six-column sheets preserve the original sample IDs, URLs, and control
assignments. The retired Bash pipeline and metadata script are available in Git
history. Stock Bowtie 1 comparison modules remain available in `modules/`, but
are not imported or run by the main workflow.

## Install and run

Install [Pixi](https://pixi.sh/), then install the project and its Python tools:

```bash
pixi run install-all
```

Place each sheet's gzipped reference FASTA under `references/`, using the exact
basename in `reference_fasta`. Reference genomes are not downloaded automatically.
Nextflow builds missing chromosome-size files and bwa-mem2 indexes once per
reference, and copies them back to that directory. The first resume after generating references can rerun downstream
tasks as inputs move to their published cache paths; subsequent resumes reuse
those task results. Existing nonempty size files
and complete indexes are reused. Remove derived files and use a fresh work
directory when replacing a reference; cached reference products are not checked
against the FASTA's contents.

Run exactly one sheet per invocation:

```bash
# Show wrapper options.
pixi run pipeline --help


# Run test data with default profile: inferred from the environment (see below).
pixi run pipeline test-Sorghum_bicolor-PRJNA1177471.tsv

# Local execution with micromamba/conda environments.
pixi run pipeline -p local,conda 15-Arabidopsis_lyrata-PRJNA1177479.tsv

# Continue a previous run using Nextflow's task cache.
pixi run pipeline 15-Arabidopsis_lyrata-PRJNA1177479.tsv -resume

# Submit tasks from Sherlock using accessible container images.
pixi run pipeline -p sherlock,apptainer 15-Arabidopsis_lyrata-PRJNA1177479.tsv
```

Run through the root Pixi environment so Nextflow can find micromamba. Without
`-p`/`--profile` (or a Nextflow `-profile`), `scripts/get-default-profile.sh`
chooses `sherlock,apptainer` where `sbatch` exists, else `local,docker` where
`docker` exists, else `local,conda`. Choose an executor and a software
environment together. `sherlock` uses SLURM and the
`normal` queue; select another queue with `--slurm_queue`. Locally built images
must be published with approval, or made available on the cluster, before
Sherlock's Apptainer can use them. Cluster work, cache, reference, and output
paths must be visible to compute nodes and writable by the user.

Each run lives in `WORKSPACE/RUN_NAME/`, where `RUN_NAME` defaults to the sheet's
basename without `.tsv` (override with `-n`/`--run-name`). The default workspace,
from `scripts/get-default-workspace.sh`, is `$SCRATCH/workspace` when `$SCRATCH`
is set and `workspace/` in this repository otherwise; override it with
`-w`/`--workspace`. A run folder contains:

- `work/`: the Nextflow work directory.
- `output/`: published results.
- `trace/`: timestamped execution report, timeline, trace, and DAG files.
- `.nextflow.log` and `.nextflow/`: the wrapper launches Nextflow from the run
  folder, so `-resume` uses that run's history.

Use `--references /absolute/path` to choose the reference location (default
`references/` in this repository). Relative paths in Nextflow arguments resolve
against the run folder. Nextflow retains intermediate FASTQs and BAMs under its
work directory for caching; retain it for `-resume` and clean it only when those
intermediates are no longer needed. Task scripts run with `bash -euxo pipefail`,
so each task's `.command.log` records the commands it ran.

```bash
# Path to the most recent execution report (optional argument: WORKSPACE).
pixi run get-latest-report
# Follow a task's log by the hash Nextflow prints, hiding `set -x` lines.
pixi run watch-log ab/123456
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
Validation completes before downloads start.

Reference names are gzipped FASTA basenames ending in `.fa.gz`, `.fasta.gz`, or
`.fna.gz`. Identifiers and filenames may contain letters, digits, underscores,
periods, and hyphens, with a letter, digit, or underscore first. Read URLs must
be HTTP or HTTPS. Reads are downloaded in batches of `--download_batch_size`
samples (default 8), one batch at a time; mates with matching basenames are fine.

## Analysis choices

Main and first-read QC alignments use **full reads**; only SPP's input is cut
(see below). Read length is inferred
from up to the first 100 reads per mate and recorded as metadata; it is not a
trimming instruction. Inspected reads must be nonempty, structurally valid, and
uniform in length across mates. Remote reads are fetched with parallel aria2c
connections and fully decompressed with bgzip, which checks gzip integrity; first-read
line counts must be divisible by four, and paired FASTQs must have equal line
counts. Genome size is the sum of chromosome sizes derived from the FASTA.

The main bwa-mem2 BAM retains primary mapped alignments, including **MAPQ 0**.
`samtools view -F 2308` excludes unmapped, secondary, and supplementary records.
There is no MAPQ threshold, NH weighting, or separate deduplication step in the
bwa-mem2 alignment. QC aligns full first reads as an unpaired library, used for
first-read statistics and tracks. A third alignment, for SPP only, uses first
reads cut to their first 50 bases (`--spp_read_length`). Each sample's alignment
provenance records the bwa-mem2 version, layout, read length, flag exclusion,
options, and the QC and SPP read policies.

BWA defaults are `-K 10000000 -k 19 -c 10000 -T 30`. Override them with
`--bwa_batch_size`, `--bwa_seed_length`, `--bwa_max_seed_occurrences`, and
`--bwa_score_threshold`; each must be a positive integer.
Parameters are validated against `nextflow_schema.json` (nf-schema) before any
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
in [vendor/phantompeakqualtools/README.md](vendor/phantompeakqualtools/README.md).
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
- SPP TSV/PDF, samtools flagstat, idxstats, and alignment statistics; PE insert-length distributions.
- CPM and 5-prime count bigWigs, including PE first-read tracks.
- Treatment MACS3 outputs, including gzipped narrowPeak files.
- Treatment `.peaks.RPM.tsv` and `.peaks.mean_RPKM.tsv` scores.

**FASTQs, BAMs, BAM indexes, and bedGraphs are excluded from published results.** Derived reference files
are published separately under the reference directory.

bwa-mem2 indexing and alignment size their own requests. Indexing reads the
gzipped FASTA directly on one CPU, with 1 GB plus 32 bytes per reference base
(it peaks near 28) and one hour plus an hour per Gbp. Alignment uses eight CPUs,
memory from the genome size, and one hour plus an hour per GB of compressed
reads. Both retry twice with doubled, then tripled, requests when killed for
memory or walltime. QC and tracks default to two CPUs, 8 GB, and four hours,
configurable with `--cpus`, `--memory`, and `--time`. Other modules declare their
own resources. The `local` profile caps requests at the machine's CPUs and memory.
The `test` profile caps every task at two CPUs, 2 GB, and five minutes; it does
not select a sample subset or replace tools. For agent-submitted Sherlock validation, every task must explicitly stay
within five minutes, two CPUs, and 8 GB. Run computation on allocated nodes.

## Environments and Docker builds

Pure tool environments are defined in `environments/*.yaml`; their matching
Pixi manifests and locks support container builds. `genesis_tools/` is a typed
Python 3.14.5 uv project that reads BAMs with pysam. Its conda equivalent installs the local
Python project; containers use the dedicated uv Dockerfile. Each Nextflow module
declares its environment and reads its image from `.env` using nf-dotenv.

Keep the original Docker discovery/build scripts and run them through Pixi:

```bash
pixi run build-dockers --no-push
pixi run build-dockers --no-push -- dap_seq_alignment genesis_tools
```

Every build targets `linux/amd64` and `linux/arm64`. The original scripts default
to pushing: always pass `--no-push` for local builds. Tags default to the latest
Git tag (`0.1.0`); use `--tag` to override. The scripts update `.env` with image
digests and preserve unrelated entries. Publishing images requires approval.
The reference Dockerfiles and pinned base images are retained.

## Verification

```bash
pixi run checks
# Opt-in real-tool validation using existing local Docker images:
pixi run validate-docker
# Retain fixtures, full logs, traces, and intermediate outputs for inspection:
pixi run validate-docker --keep
```

The default suite runs ShellCheck, Ruff, ty, quantification and mocked Docker
build regressions, sample-sheet and metadata cases, and actual Nextflow SE/PE
stub workflows. It checks control fan-out, per-reference generation and cache
reuse, `-resume`, invalid inputs, and publication boundaries. Mock builds cover
both architectures and `--no-push` without contacting Docker or a registry.

Docker validation runs the real download/staging, indexing, bwa-mem2, samtools,
SPP, deepTools, MACS3, and Python tools on deterministic 75-base SE/PE
fixtures with enriched regions and duplicated sequence. It checks nonempty peaks,
finite positive scores, provenance, read lengths (full, and 50 bases for SPP),
outputs, and resume.
Validated on 2026-10-06 with the existing linux/arm64 images: all 45 initial
tasks completed, each of the three treatments produced 200 peaks, and all 41
tasks on the stable resume were cached. BAM inspection confirmed retained MAPQ 0
alignments, full 75-base sequences, and unpaired first-read QC alignments. No
FASTQs or BAMs appeared in published outputs. Tool logs contained no warnings.

These small synthetic runs validate integration; comparisons of biological
outputs on the full public datasets remain a separate validation task. Native
macOS stub traces do not include task runtime metrics; Docker traces do.

See [docs/nextflow-style.md](docs/nextflow-style.md) for module conventions.
