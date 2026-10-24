# Paired-end bulk ATAC

Complete the [one-time installation](../../README.md#install-and-check-your-environment)
first. Bulk ATAC runs like DAP-seq: the same launcher, profiles, workspace layout,
reference cache and image pins, with `workflows/atac-seq.nf` as the workflow and `conf/atac-seq.config` for its settings.

```bash
pixi run pipeline-atac -w /scratch/genesis -n atac-A /data/atac.tsv --references /data/references
# Continue the same analysis after an interruption.
pixi run pipeline-atac -w /scratch/genesis -n atac-A /data/atac.tsv --references /data/references -resume
```

`pixi run pipeline-atac` is `scripts/run-pipeline.sh --workflow atac-seq`; run either with
`--help` for the options. Pixi runs from the checkout folder, so give absolute paths.

## Sample sheet

One row per paired-end library, with this exact tab-separated header:

```tsv
sample_id	species	read1_url	read2_url	reference_fasta	reference_genes	organelle_contigs
leaf_rep1	Arabidopsis_thaliana	https://example.org/rep1_1.fastq.gz	https://example.org/rep1_2.fastq.gz	TAIR10.fa.gz	TAIR10.genes.bed	ChrC,ChrM
```

| Column | Meaning |
| --- | --- |
| `sample_id`, `species` | Library identifier and output folder, as for DAP-seq |
| `read1_url`, `read2_url` | HTTP(S) URLs of both gzipped mates; single-end libraries are rejected |
| `reference_fasta` | Gzipped FASTA in the references folder, cached as for DAP-seq |
| `reference_genes` | BED6 gene or transcript intervals in the references folder, or `-` to skip TSS enrichment |
| `organelle_contigs` | Comma-separated mitochondrial and plastid contig names, or `-` for none |

Rows that share a reference must give the same `reference_genes` and
`organelle_contigs`. The whole sheet is checked before any task starts, and a
named organelle contig that is not in the FASTA stops the run. Each row is one
library: biological replicates are never pooled.

The gene file only needs a strand per interval; the TSS is the start of `+`
intervals and the end of `-` intervals. From a GFF3, for example:

```bash
awk -v OFS='\t' '$3 == "gene" {print $1, $4 - 1, $5, ".", 0, $7}' genes.gff3 > genes.bed
```

## Analysis and metrics

The steps follow the ENCODE ATAC pipeline:

1. **Trimming.** cutadapt removes the Nextera adapter (`--atac_adapter`, default
   `CTGTCTCTTATACACATCT`) from both mates and drops pairs where either mate is
   shorter than 5 bases. FastQC runs on the raw mates.
2. **Alignment.** bwa-mem2 with the same options as DAP-seq, then samtools fixmate
   and markdup mark duplicates.
3. **Filtering.** A usable pair has both mates primary, properly paired, not
   QC-failed, on a nuclear contig and at least MAPQ `--atac_mapq` (default 30).
   A mate whose partner fails is removed too. Duplicates are removed last.
4. **Fragments.** Each usable pair's outer span, BED0 half-open, sorted and
   tabix-indexed.
5. **Peaks.** MACS3 BAMPE on the usable fragments, `-q 0.01 --keep-dup all`, with
   the summed nuclear contig lengths as the genome size.
6. **Cut sites.** deepTools `alignmentSieve --ATACshift` shifts reads +4/−5; the
   bigWig holds raw per-base counts of the shifted 5′ ends, that is BED0
   `start + 4` and `end − 6` of each fragment.

The `qc/*.atac_qc_mqc.tsv` table, also shown in MultiQC, reports per library:

| Metric | Definition |
| --- | --- |
| `FRiP` | Usable reads overlapping a peak, over all usable reads |
| `TSS_enrichment` | Mean cut count at the TSS over the mean of the outer 100 bases at each end of the ±2 kb profile; `NA` without a gene file |
| `NRF`, `PBC1`, `PBC2` | ENCODE library complexity over usable fragments before duplicate removal |

Organelle reads are visible in the samtools idxstats of the marked BAM. No plant
QC thresholds are applied.

## Outputs

Under `WORKSPACE/RUN_NAME/output/<species>/<sample_id>/`:

- `*.fragments.bed.gz` and its `.tbi` index;
- `*.macs3_peaks.narrowPeak.gz` and the other MACS3 outputs;
- `*.cuts.counts.bw`: raw Tn5 cut counts;
- `qc/`: FastQC, the cutadapt report, samtools reports for the marked and usable
  alignments, `*.complexity.tsv`, the ATAC QC table and the TSS profile.

`output/multiqc/` has the run's MultiQC report. As for DAP-seq, `output/` also
holds `dap-images.lock` and `git-status.txt`, and `trace/` the Nextflow trace,
report, timeline and DAG. FASTQs and BAMs stay in `work/`; keep it for resume.
