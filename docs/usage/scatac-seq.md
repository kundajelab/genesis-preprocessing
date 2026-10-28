# scATAC pseudobulks from fragments

Complete the [one-time installation](../../README.md#install-and-check-your-environment)
first. scATAC runs like DAP-seq and bulk ATAC: the same launcher, profiles, workspace
layout, reference folder and image pins, with `workflows/scatac-seq.nf` as the workflow
and `conf/scatac-seq.config` for its settings.

It starts from published single-cell fragment files, such as the 10x Cell Ranger ATAC
`fragments.tsv.gz`, not from reads. Each library's cells are grouped into pseudobulks,
for example one per cell type, and each pseudobulk gets peaks, a cut-site track and QC.

```bash
pixi run pipeline-scatac -w /scratch/genesis -n scatac-A /data/scatac.tsv --references /data/references
# Continue the same analysis after an interruption.
pixi run pipeline-scatac -w /scratch/genesis -n scatac-A /data/scatac.tsv --references /data/references -resume
```

`pixi run pipeline-scatac` is `scripts/run-pipeline.sh --workflow scatac-seq`; run either
with `--help` for the options. Pixi runs from the checkout folder, so give absolute paths.

## Sample sheet

One row per library, with this exact tab-separated header:

```tsv
sample_id	species	fragments_url	cells_url	reference_fasta	reference_genes	organelle_contigs
root_rep1	Arabidopsis_thaliana	https://example.org/rep1.fragments.tsv.gz	https://example.org/rep1.cells.tsv	TAIR10.fa.gz	TAIR10.genes.bed	Mt,Pt
```

| Column | Meaning |
| --- | --- |
| `sample_id`, `species` | Library identifier and output folder, as for DAP-seq |
| `fragments_url` | HTTP(S) URL of the library's gzipped fragment file |
| `cells_url` | HTTP(S) URL of the library's cell list (below) |
| `reference_fasta` | Gzipped FASTA in the references folder, cached as for DAP-seq; the fragments must be on this assembly |
| `reference_genes` | BED6 gene or transcript intervals in the references folder, or `-` to skip TSS enrichment |
| `organelle_contigs` | Comma-separated mitochondrial and plastid contig names, or `-` for none |

Rows that share a reference must give the same `reference_genes` and
`organelle_contigs`, as in bulk ATAC. Each row is one library: libraries, and so
biological replicates, are never pooled.

The cell list is tab-separated with a `barcode	group` header and one row per cell; a
barcode listed twice stops the library. Barcodes must match the fragment file exactly,
including any `-1` suffix. Groups are
names such as cell types (letters, digits, `_`, `.` and `-`); use one group, for
example `all`, for a single pseudobulk of every cell.

```tsv
barcode	group
AAACGAAAGAAACGCC-1	epidermis
AAACGAAAGACATTGG-1	cortex
```

## Analysis and metrics

Each pseudobulk is processed like a bulk ATAC library, with the differences from
single-cell tools listed below:

1. **Fragments.** Every row must be on a reference contig and inside its length, so
   fragments from another assembly stop the run. Rows of unlisted barcodes and of
   organelle contigs are counted and dropped. Coordinates are kept exactly as written; 10x files are BED0
   half-open with the Tn5 offsets already applied.
2. **Pseudobulks.** Each group's fragments in each library, one row per fragment
   whatever its read-pair count.
3. **Peaks.** MACS3 BEDPE on the pseudobulk's fragments, `-q 0.01 --keep-dup all`,
   with the summed nuclear contig lengths as the genome size.
4. **Cut sites.** The first and last base of each fragment, the same Tn5 insertion
   bases bulk ATAC counts after its +4/−5 shift; the bigWig holds raw per-base counts.

The library table, `qc/*.scatac_library_mqc.tsv`, counts source fragments, those of
unlisted barcodes and of organelles, the listed cells and their kept fragments.
`qc/*.cells.tsv` has, per cell, its group, nuclear fragments, their read-pair support
and the organelle fraction. The pseudobulk table, also shown in MultiQC, reports:

| Metric | Definition |
| --- | --- |
| `cells`, `fragments` | Cells with a kept fragment, and their fragments |
| `insertions`, `FRiP` | Tn5 insertions (two per fragment), and the fraction inside a peak |
| `TSS_enrichment` | As in bulk ATAC: mean cut count at the TSS over the mean of the outer 100 bases at each end of the ±2 kb profile; `NA` without a gene file |

No cell filtering or plant QC thresholds are applied; the cell list decides which
barcodes count.

### Differences from single-cell tools

These follow bulk ATAC and are not yet compared with single-cell tools on plant data.
ArchR is the example here:

- **Cells.** The cell list is taken as given; ArchR filters cells by TSS enrichment
  and fragment count (by default at least 4 and 1,000).
- **Peaks.** MACS3 `-q 0.01` in BEDPE mode on whole fragments, as in bulk ATAC; ArchR
  runs MACS2 on Tn5 insertions extended around each cut site.
- **Peak sets.** Each pseudobulk keeps its own peaks; there is no merged peak set
  across groups and no pseudo-replicates within a group.
- **TSS enrichment.** One value per pseudobulk, as in bulk ATAC; ArchR scores each cell.

## Outputs

Under `WORKSPACE/RUN_NAME/output/<species>/<sample_id>/`:

- `*.fragments.tsv.gz` and its `.tbi` index: the kept fragments of listed cells;
- `qc/`: the library table and the per-cell table;
- `groups/<group>/`: the pseudobulk's MACS3 outputs, `*.cuts.counts.bw` and its QC
  table and TSS profile.

`output/multiqc/` has the run's MultiQC report. As for DAP-seq, `output/` also holds
`genesis-images.lock` and `git-status.txt`, and `trace/` the Nextflow trace, report,
timeline and DAG.
