# Paired-end bulk ATAC

Complete the [one-time installation](../../README.md#install-and-check-your-environment)
first. Pixi supplies the scientific runtime internally.

```bash
pixi run pipeline-atac /path/to/samples.tsv --references /path/to/references.json --outdir /path/to/results/atac
pixi run pipeline-atac /path/to/samples.tsv --references /path/to/references.json --outdir /path/to/results/atac --resume
```

Pixi runs from the checkout folder, so give absolute paths.

Each row describes one paired-end lane. Repeated `library_id` values concatenate
technical lanes before alignment; all library metadata and analysis choices must
agree. Separate libraries and biological replicates are never pooled automatically.
IDs use letters, numbers, dots, underscores and dashes.

## Sample sheet

Use these exact tab-separated columns:

| Column | Meaning |
| --- | --- |
| `library_id` | Unique library identifier |
| `sample_id` | Biological context used to group library comparisons |
| `biological_replicate` | Explicit replicate identifier within that context |
| `lane_id` | Unique technical lane within a library |
| `reference_id` | Exact entry in the reference registry |
| `read1`, `read2` | Local gzipped FASTQs or public HTTPS URLs |
| `read1_sha256`, `read2_sha256` | SHA256 of each compressed FASTQ |
| `mapq` | Minimum MAPQ for both primary mates, 0–254 |
| `duplicates` | `retain` or `exclude` marked duplicates |
| `adapter_r1`, `adapter_r2` | Adapter sequences, or `-` in both fields for no trimming |

Authenticated downloads must be staged locally; credentials and signed URLs do
not belong in manifests. FASTQs are checked for pairing and integrity. Raw mates
receive separate FastQC HTML and ZIP reports before any chosen adapter processing.
Warnings in FastQC are diagnostic and do not reject a dataset.

Adapter sequences invoke Cutadapt's paired adapter removal without additional
quality or length filtering. The `-` policy uses full reads. These are explicit
scientific choices, not automatic responses to a QC warning.

## Reference registry

The registry is JSON with `schema_version: 1` and a `references` list. Each entry
requires:

- `reference_id`, `species`, `assembly`, `annotation_release`;
- `fasta` and `fasta_sha256`: exact gzipped FASTA and compressed-file checksum;
- `tss` and `tss_sha256`: annotation-derived positions and checksum;
- `genome_size` and `genome_size_method`: positive MACS3 effective genome size
  and its derivation;
- `mitochondrial_contigs` and `plastid_contigs`: separate lists of exact names.

TSS files have three tab-separated columns: contig, **zero-based position**, and
strand (`+` or `-`). Derive positions from the recorded annotation release; the
workflow does not guess a GTF/GFF or infer organelles from contig names. Empty
organelle lists explicitly declare none. Do not omit known organellar contigs.

All local inputs are hashed before execution. Reference index tasks use deep
content caching. Changes to a run's sheet or registry require a new output
folder, preserving earlier results. Task helpers are retained in immutable
source bundles with an explicit SHA256 cache input; changing helper code
invalidates affected tasks. The two-stage cache used by DAP-seq is
separate; ATAC indexes remain in their Nextflow work directory.

## Analysis and metrics

bwa-mem2 aligns reads with a library read group. samtools fixmate and markdup mark
duplicates in an intermediate BAM without removing them. The explicit duplicate
policy is applied when selecting fragments.

A usable fragment requires exactly two primary, mapped, proper mates on the
same nuclear contig, neither QC-failed, with both MAPQs meeting the selected
threshold. MAPQ 255 means unknown and is excluded. Secondary and supplementary
records are counted separately. Invalid pairs and spans shorter than ten bases
are reported. Mitochondrial and plastid pairs are reported separately before
MAPQ filtering and excluded from nuclear analysis.

Fragment coordinates are the outer aligned span, BED0 half-open. MACS3 calls
BAMPE peaks at q=0.01 with `--keep-dup all`, because the preceding selection has
already applied the declared duplicate policy. Cut sites are separate one-base
intervals at `start + 4` and `end - 6`, corresponding to the usual +4/−5 shift
from the forward/reverse aligned 5′ bases. Cut-site bigWigs contain raw counts,
not normalized fragment coverage.

FRiP counts a fragment once when its span overlaps any called peak; the
denominator is all usable fragments. TSS profiles use strand-aware ±2 kb windows
that fit inside the contig. Enrichment is the cut count at the single center
position divided by the mean count in the outer 100 positions on each side.
The stored profile is unsmoothed; a display that smooths the curve must say so.
The metric definition and eligible TSS count are retained with the result.
Undefined metrics are null, not zero. NRF and PBC use nuclear, MAPQ-passing fragment coordinates before duplicate
exclusion; coordinate duplication is not proof of PCR origin.

Libraries with matching sample and reference IDs receive pairwise peak
base-pair Jaccard comparisons. This is descriptive overlap, not IDR or a biological
replicate acceptance criterion. No peak pooling, universal plant threshold,
blacklist, read-depth normalization or single-cell processing is implicit.

## Outputs

Each library has:

- `qc/`: two FastQC reports, raw samtools metrics, duplicate-marking log,
  adapter report, enrichment metrics and MultiQC HTML plus parsed data;
- `fragments/`: indexed `fragments.bed.gz`, cut-count BEDGraph, usable samtools
  metrics, fragment lengths, MAPQ distributions and complexity metrics;
- `tracks/cuts.counts.bw` and `peaks/atac_peaks.narrowPeak`;
- `provenance/acquisition.json`: source checksums and template count.

Run outputs include `libraries.tsv`, `replicate-comparisons.json`, checksums,
Nextflow trace/report/timeline/DAG and invocation provenance. BAMs and FASTQs stay
in the work directory. Keep that directory for resume. Invocation records retain
commands, Git SHA, source checksums, image digests, exit code and elapsed time.

## Scope and deployment

The entry point is `workflows/atac.nf`; reusable tasks are in `modules/atac/`,
execution profiles in `conf/`, and typed helpers in `genesis_tools/atac/`.
The historical prototype in `experimental/atac/` is separate from `pipeline-atac`.

For clusters and isolated dataset campaigns, see the execution guide in
[genesis-curation](https://github.com/kundajelab/genesis-curation/blob/main/docs/usage/execution.md). Large genomes need substantial index memory and scratch space;
Nextflow requests reference-dependent resources before starting indexing.
Computational validation and biological dataset quality are separate decisions.
