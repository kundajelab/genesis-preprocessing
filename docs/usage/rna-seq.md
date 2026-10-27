# Bulk RNA-seq

Complete the [one-time installation](../../README.md#install-and-check-your-environment)
first. RNA-seq runs like DAP-seq and bulk ATAC: the same launcher, profiles, workspace
layout, reference folder and image pins, with `workflows/rna-seq.nf` as the workflow and
`conf/rna-seq.config` for its settings.

```bash
pixi run pipeline-rna -w /scratch/genesis -n rna-A /data/rna.tsv --references /data/references
# Continue the same analysis after an interruption.
pixi run pipeline-rna -w /scratch/genesis -n rna-A /data/rna.tsv --references /data/references -resume
```

`pixi run pipeline-rna` is `scripts/run-pipeline.sh --workflow rna-seq`; run either with
`--help` for the options. Pixi runs from the checkout folder, so give absolute paths.

## Sample sheet

One row per library, with this exact tab-separated header:

```tsv
sample_id	species	read1_url	read2_url	strandedness	reference_fasta	reference_annotation	organelle_contigs
leaf_rep1	Arabidopsis_thaliana	https://example.org/rep1_1.fastq.gz	https://example.org/rep1_2.fastq.gz	reverse	TAIR10.fa.gz	Araport11.gff3	ChrC,ChrM
```

| Column | Meaning |
| --- | --- |
| `sample_id`, `species` | Library identifier and output folder, as for DAP-seq |
| `read1_url`, `read2_url` | HTTP(S) URLs of the gzipped reads; `read2_url` is `-` for single-end libraries |
| `strandedness` | `unstranded`, `forward` (read 1 is the transcript's strand) or `reverse` (dUTP and most stranded kits) |
| `reference_fasta` | Gzipped FASTA in the references folder, cached as for DAP-seq |
| `reference_annotation` | The annotation as published, GFF3 or GTF, optionally gzipped, in the references folder |
| `organelle_contigs` | Comma-separated mitochondrial and plastid contig names, or `-` for none |

Rows that share a reference must give the same `reference_annotation` and
`organelle_contigs`. The whole sheet is checked before any task starts. Annotation
contigs that are not in the FASTA, as when the annotation says `Chr1` and the FASTA
`1`, stop the run, and so does a named organelle contig that is not in the FASTA.
Each row is one library: biological replicates are never pooled.

## Analysis and metrics

1. **Annotation.** gffread converts the GFF3 or GTF (Araport11, MaizeGDB and
   Phytozome publish GFF3) to one GTF, with each gene's ID as `gene_id`. Lifted
   annotations such as Col-CEN Araport11 (Liftoff) reuse a transcript ID for gene
   copies at other loci, for example across rDNA and tRNA arrays; each copy gets its
   own transcript ID and keeps the gene ID, so STAR does not join the copies into one
   transcript and all copies count toward the gene. The `gene:` prefix of Ensembl
   Plants GFF3 IDs is dropped, so genes keep their usual IDs, such as `AT1G01010`.
2. **Index.** STAR indexes the FASTA with the annotated junctions, once per reference
   in each run. Its memory is sized from the genome: about 35 GB for maize B73.
3. **Alignment.** STAR, with no intron longer than the longest annotated intron
   (`--alignIntronMax`; the mate gap may add 1 kb for the unread middle of the
   insert). STAR's default allows about 590 kb, sized for mammals. Annotated plant
   introns are much shorter: 99.9% are under about 2 kb in Araport11, 4 kb in
   A. lyrata and 20 kb in maize, although each annotation also has a few far longer
   ones (14 kb in A. lyrata, 435 kb in maize), which set the limit. Reads may align to up to 20 loci, as in ENCODE. Read ends that do
   not match, such as adapters, are soft-clipped; there is no separate trimming.
4. **Gene counts.** featureCounts on exons, summed per gene, with the sheet's
   strandedness; a pair counts once. Reads that align to several loci, or overlap
   several genes, are not counted.
5. **Tracks.** A per-base bigWig of uniquely aligned reads in CPM (deepTools
   `bamCoverage`), both strands together; spliced reads cover their exons only.

The `qc/*.rna_qc_mqc.tsv` table, also shown in MultiQC beside the STAR and
featureCounts reports, has per library:

| Metric | Definition |
| --- | --- |
| `forward_fraction` | Of the reads featureCounts assigns in either orientation, the share assigned as `forward`: near 1 forward, near 0 reverse, near 0.5 unstranded. Compare it with `strandedness`. |
| `input_reads`, `unique`, `multimapped` | From the STAR log: reads (pairs) in, and the fractions aligned to one locus and to several |
| `assigned` | Fraction of reads (pairs) counted to a gene, from the featureCounts summary |
| `rRNA` | Fraction of assigned reads in rRNA genes: `rRNA` features, or an `rRNA` biotype in Ensembl GTFs; `NA` when the annotation has none, as for MaizeGDB and Phytozome gene sets |
| `organelle` | Fraction of uniquely aligned reads on the organelle contigs |
| `max_intron` | The intron limit used, from the annotation |

Organelle and rRNA genes stay in the gene counts; their fractions are reported so
that libraries dominated by chloroplast or rRNA reads stand out. No plant QC
thresholds are applied.

### Differences from ENCODE and nf-core/rnaseq

These are deliberate but not yet compared with those pipelines on plant data:

- **Intron limits.** From the annotation, as above; ENCODE uses 1,000,000 bases for
  both introns and mate gaps.
- **Quantification.** featureCounts gene counts of unique alignments; ENCODE uses RSEM
  and nf-core/rnaseq Salmon by default, which share multi-mapping reads between
  genes. Identical homeologs in polyploids, and repeated gene families, therefore get
  fewer counts here; `multimapped` shows how many reads that concerns.
- **Strandedness.** Taken from the sheet and checked through `forward_fraction`;
  nf-core/rnaseq can infer it with Salmon.
- **Trimming and duplicates.** No adapter trimming (Trim Galore in nf-core/rnaseq) and
  no duplicate marking.
- **Tracks.** One unstranded CPM track; ENCODE also writes strand-specific signal.
- **Index reuse.** The STAR index is rebuilt in each new run directory; `-resume` reuses
  it within a run.

## Outputs

Under `WORKSPACE/RUN_NAME/output/<species>/<sample_id>/`:

- `*.featureCounts.tsv`: gene counts, with each gene's chromosome, exon starts, ends,
  strands and length;
- `*.cpm.bw`: the CPM track;
- `qc/`: the STAR log and junctions, the featureCounts summary, the RNA-seq QC table
  and FastQC reports.

`output/multiqc/` has the run's MultiQC report. As for DAP-seq, `output/` also holds
`dap-images.lock` and `git-status.txt`, and `trace/` the Nextflow trace, report,
timeline and DAG.
