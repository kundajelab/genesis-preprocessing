# scATAC pseudobulks from fragments

Complete the [one-time installation](../../README.md#install-and-check-your-environment)
first. scATAC runs like DAP-seq and bulk ATAC: the same launcher, profiles, workspace
layout, reference folder and image pins, with `workflows/scatac-seq.nf` as the workflow
and `conf/scatac-seq.config` for its settings.

![Genesis scATAC pseudobulk workflow metro map](../../docs/images/genesis_scatac_metro_map.svg)

[Diagram source and interpretation](../../docs/diagrams/README.md)

It starts from published single-cell fragment files, such as the 10x Cell Ranger ATAC
`fragments.tsv.gz`; [another workflow](#from-raw-10x-reads) starts from raw 10x reads. Each library's cells are grouped into pseudobulks,
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

## From raw 10x reads

Most plant scATAC studies publish only reads. `workflows/scatac-reads.nf` (`pixi run
pipeline-scatac-reads`, `conf/scatac-reads.config`) makes each library's fragments from its
10x reads on NCBI SRA, calls cells, then runs the steps above unchanged.

![Genesis scATAC from 10x reads workflow metro map](../../docs/images/genesis_scatac_reads_metro_map.svg)

```bash
pixi run pipeline-scatac-reads -w /scratch/genesis -n scatac-B /data/scatac-reads.tsv --references /data/references
```

```tsv
sample_id	species	sra_runs	chemistry	cells_url	reference_fasta	reference_genes	organelle_contigs	barcode_whitelist
root_rep1	Arabidopsis_thaliana	SRR12340627,SRR12340628	10x-atac-v1	-	Athaliana_447_TAIR10.fa.gz	Athaliana_447_Araport11.genes.bed	ChrC,ChrM	737K-cratac-v1.txt
```

| Column | Meaning |
| --- | --- |
| `sra_runs` | The library's SRA runs (its lanes), comma-separated; they are mapped together |
| `chemistry` | `10x-atac-v1` or `10x-atac-v2` (16-base barcode read), or `10x-multiome-atac` (24 bases: an 8-base spacer, then the barcode) |
| `cells_url` | The authors' `barcode<TAB>group` cell list, or `-` to call cells |
| `barcode_whitelist` | The 10x barcode list in the references folder, or `-` for no barcode correction |

The other columns are as above. A run may be listed for one library only.

1. **Download.** `prefetch` and `fasterq-dump --include-technical` from NCBI SRA, which keeps
   the 10x barcode read; ENA's FASTQs drop it. Submitters order the reads differently, so the
   barcode read is the one of the chemistry's length and the genomic reads are the two
   longest. A run without such a barcode read stops the library. Downloading and converting
   are separate tasks, so a failed conversion, or `-resume` after one, reuses the download;
   a failed download is tried again. `fasterq-dump` needs temporary disk of about the FASTQs'
   size, about 100 GB for a 116-million-spot run, and takes about 1.5 hours with 4 threads.
   Intermediate FASTQS (tens of GB) are also stored in work/ until the run is cleaned up.
2. **Fragments.** chromap `--preset atac` trims adapters, maps the pairs, corrects barcodes
   against the whitelist and removes duplicates per cell. It writes 10x fragments: BED0, Tn5
   offsets applied (start + 4, end − 5). Sequencers read the barcode on either strand; the
   strand that matches the whitelist for most of the first 200,000 barcodes is used, and a
   library where neither matches half of them stops. The chromap index, about 5 bytes per
   reference base, is built once per reference and published beside the FASTA as
   `<id>.chromap/`, with `<id>.chromap.sha256` recording the FASTA's SHA-256, the chromap
   release and the index's size. Later runs reuse it only while that record matches the FASTA,
   the release pinned in `environments/GENESIS_CHROMAP.yaml` and the index file; otherwise it
   is rebuilt.
3. **Cells.** With `cells_url`, the authors' list is used, after removing Cell Ranger's `-1`
   suffix, which chromap does not write. Otherwise barcodes with at least
   `--scatac_min_fragments` nuclear fragments (default 1000, as in ArchR) are cells, at most
   `--scatac_max_cells` of them (0, the default, for no limit), all in one group, `all`.
   `qc/<sample>.barcodes.tsv` lists every barcode's nuclear and organelle fragments, and
   `<sample>.all_barcodes.fragments.tsv.gz` keeps every barcode's fragments for calling cells
   again with other tools.

10x does not allow its barcode lists to be redistributed, so they are not in this
repository and are never downloaded by the workflow. Copy the list for the chemistry (for
10x ATAC v1 and v2, `737K-cratac-v1.txt` from Cell Ranger ATAC) into the references folder.

Differences from Cell Ranger ATAC and ArchR:

- **Cell calling.** A fragment-count threshold only; Cell Ranger fits a model to the barcode
  counts and removes gel bead doublets, and ArchR also requires a TSS enrichment, which is
  not tuned for plants. Use the authors' cell list when one is published.
- **Barcodes.** chromap corrects barcodes one mismatch away from the list; Cell Ranger also
  uses base qualities.
- **Plate-based designs** (sci-ATAC-seq) have no 10x barcode read and are not supported.
