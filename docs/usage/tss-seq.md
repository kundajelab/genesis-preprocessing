# TSS mapping

Complete the [one-time installation](../../README.md#install-and-check-your-environment)
first. TSS mapping runs like the other workflows: the same launcher, profiles, workspace
layout, reference folder and image pins, with `workflows/tss-seq.nf` as the workflow and
`conf/tss-seq.config` for its settings.

It handles single-read methods, where the 5′ end of read 1 is the transcription start
site (TSS): cap-trapping CAGE, oligo-capping (TSS-seq) and csRNA-seq. Each library gets
per-base TSS tag counts (CTSSs), TSS clusters and QC. Paired-end methods (RAMPAGE,
Smar2C2, TIF-seq) and template-switching methods (nanoCAGE, nanoPARE), whose reads can
start with several added Gs, are not supported yet.

```bash
pixi run pipeline-tss -w /scratch/genesis -n tss-A /data/tss.tsv --references /data/references
# Continue the same analysis after an interruption.
pixi run pipeline-tss -w /scratch/genesis -n tss-A /data/tss.tsv --references /data/references -resume
```

`pixi run pipeline-tss` is `scripts/run-pipeline.sh --workflow tss-seq`; run either with
`--help` for the options. Pixi runs from the checkout folder, so give absolute paths.

## Sample sheet

One row per library, with this exact tab-separated header:

```tsv
sample_id	species	read1_url	method	trim_5p	adapter_3p	reference_fasta	reference_annotation	organelle_contigs
wt_rep1	Arabidopsis_thaliana	https://example.org/wt1.fastq.gz	CAGE	9	TCGTATGCCGTCTTCTGCTTG	TAIR10.fa.gz	Araport11.gff3	ChrC,ChrM
```

| Column | Meaning |
| --- | --- |
| `sample_id`, `species` | Library identifier and output folder, as for DAP-seq |
| `read1_url` | HTTP(S) URL of the gzipped reads that start at the TSS (read 1 for paired libraries) |
| `method` | `CAGE` or `nAnT-iCAGE` (cap-trapping), `TSS-seq` or `oligo-capping`, or `csRNA-seq` |
| `trim_5p` | Bases to remove from the start of every read, such as a sample barcode or linker, or `-` |
| `adapter_3p` | The 3′ adapter sequence, or `-` for none |
| `reference_fasta`, `reference_annotation`, `organelle_contigs` | As for RNA-seq |

The fixed 5′ bases differ between studies, not only between methods: check the start of
a few reads before writing the sheet. GSE136356's CAGE reads, for example, start with a
3-base sample barcode and `CAGCAG` (`trim_5p` 9) and end in the Illumina small-RNA
adapter. Rows that share a reference must give the same `reference_annotation` and
`organelle_contigs`, and the whole sheet is checked before any task starts.

## Analysis and metrics

1. **Annotation and index.** As for RNA-seq: `RNA_ANNOTATION` writes the GTF and the
   longest-intron limit, and `GET_STAR_INDEX` builds or reuses the STAR index.
2. **Trimming.** cutadapt removes `trim_5p` bases from the 5′ end and the 3′ adapter,
   quality-trims only the 3′ end (Q20) and keeps reads of at least 20 bases (18 for
   csRNA-seq). The 5′ end is the TSS, so it is never quality-trimmed.
3. **Alignment.** STAR with the annotated junctions and `--alignEndsType
   Extend5pOfRead1`, so the 5′ end is never soft-clipped. Only uniquely aligned reads
   (MAPQ 255) give tags.
4. **TSS tags.** Each read's 5′ end on its strand. Cap-trapping CAGE adds a
   non-templated G before the first base; for `CAGE` and `nAnT-iCAGE`, a read that starts
   with a G the genome lacks (from the MD tag) is moved 1 base downstream, as in FANTOM.
   A G that the genome also has cannot be told apart and stays.
5. **Clusters.** paraclu clusters each strand's tags; FANTOM's paraclu-cut then keeps
   clusters of at most `--tss_cluster_max_length` bases (200) whose density rises at
   least 2-fold, and drops clusters inside larger ones. Clusters need at least
   `--tss_cluster_min_tpm` tags per million (1, and at least 2 tags). Each cluster
   reports its tags, TPM and dominant TSS.

The `qc/*.tss_qc_mqc.tsv` table, also shown in MultiQC, reports per library:

| Metric | Definition |
| --- | --- |
| `input_reads`, `trimmed_reads`, `unique_reads` | Reads before and after trimming, and uniquely aligned |
| `tags` | Tags at TSSs: the unique reads |
| `near_TSS_100`, `near_TSS_500` | Fraction of tags within 100 or 500 bases of an annotated transcript start, on its strand |
| `antisense_500` | Fraction of tags within 500 bases of an annotated start on the opposite strand |
| `organelle` | Fraction of tags on organelle contigs |
| `rRNA` | Fraction of tags in rRNA genes' exons; `NA` when the annotation names none (Phytozome) |
| `clusters`, `in_clusters` | Number of clusters, and the fraction of tags in them |

No plant QC thresholds are applied.

### Differences from FANTOM and CAGEr

- **Aligner.** STAR with the annotation, so reads spanning a splice site near the TSS
  align; FANTOM5 used Delve, and many plant CAGE studies use Bowtie on 25-base tags.
- **G correction.** Applied only to cap-trapping methods, from the MD tag; CAGEr's
  `removeFirstG` does the same for BSgenome-based input.
- **Clusters.** paraclu as in FANTOM5; CAGEr also offers distance-based clustering
  (`distclu`) and consensus clusters across libraries, which are not made here. Each
  library is clustered on its own; replicates are never pooled.

## Outputs

Under `WORKSPACE/RUN_NAME/output/<species>/<sample_id>/`:

- `*.ctss.tsv.gz`: tags per TSS (contig, strand, 1-based position, count);
- `*.ctss.plus.bw`, `*.ctss.minus.bw`: tags per base on each strand, and `*.cpm.bw`
  versions in tags per million;
- `*.tss_clusters.bed`: clusters as BED6 (score = tags, capped at 1000 as BED requires)
  plus TPM, the dominant TSS (0-based) and its tags, and the cluster's exact tags;
- `qc/`: the cutadapt and STAR logs, the TSS QC table and FastQC reports.

`output/multiqc/` has the run's MultiQC report. As for DAP-seq, `output/` also holds
`genesis-images.lock` and `git-status.txt`, and `trace/` the Nextflow trace, report,
timeline and DAG. FASTQs and BAMs stay in `work/`.
