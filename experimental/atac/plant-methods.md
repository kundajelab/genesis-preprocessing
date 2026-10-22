# Historical prototype: methods and limitations

This page describes the frozen experimental workflow. For current behavior, see
the [bulk ATAC guide](../../docs/usage/bulk-atac.md); for package comparisons, see
the [benchmark guide](https://github.com/kundajelab/genesis-curation/blob/main/benchmarks/README.md).

The integrated prototype preserves the scientific behavior of commit
`0b8850207e94a3d4909bf4ab99d18c03cd010cc2`. Its input is one synthetic,
adapter-free PE library. Computational checks do not establish plant biological
quality. All Genesis plant pass/fail thresholds remain **UNSPECIFIED**.

## Evidence and decisions

| Topic | Primary evidence | Frozen prototype behavior | Remaining work |
| --- | --- | --- | --- |
| Reproducibility and QC | [ENCODE ATAC standards](https://www.encodeproject.org/atac-seq/) assess usable reads, complexity, fragment structure, enrichment and reproducibility; TSS interpretation depends on the reference | Fragment counts, NRF/PBC, FRiP, fragment lengths and a toy TSS profile are diagnostic | Biological replicates, concordance/IDR and reference-aware TSS implementation before production |
| Alignment and filtering | [ENCODE code](https://github.com/ENCODE-DCC/atac-seq-pipeline/tree/47ba8dff9c332e24b48e767303e9fcac98589cf2) separates alignment, MAPQ, duplicate and chromosome filtering | bwa-mem2, samtools duplicate marking; both mates must meet the explicitly supplied MAPQ; marked duplicates excluded from usable fragments | Compare explicit alternatives empirically on exact plant assemblies; no claim of identity with ENCODE Bowtie2/Picard |
| Cut sites versus fragments | [ENCODE bam2ta](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/47ba8dff9c332e24b48e767303e9fcac98589cf2/src/encode_task_bam2ta.py) applies strand-specific +4/−5 shifts | BED0 cut bases at fragment start+4 and end−6; MACS3 BAMPE consumes unshifted fragments | Keep these representations distinct; known-answer tests check endpoints |
| Organellar contamination | [Lu et al., 2017](https://academic.oup.com/nar/article/45/6/e41/2605943) demonstrate nuclei sorting and organellar contamination issues in plant ATAC | Explicit contig names, accounting and exclusion; no guessed chromosome names | Current counter combines mitochondrial/plastid fragments and follows MAPQ filtering. Report each separately and raw denominators in a future tested change |
| Preparation and tissue context | [Maher et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5810565/) examine multiple plant species and cell types using nuclei isolation | Reference and synthetic input provenance retained | Add preparation, tissue, genotype and replicate metadata before real-data comparisons |
| Fragment profiles and controls | [Hsieh et al., 2024](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1370618/full) evaluate maize nuclei preparation, organellar filtering, fragment periodicity and protocol-specific gDNA controls | Fragment distribution retained; no control input in this prototype | Evaluate controls where experimental design supports them; do not adopt one study's peak threshold or fragment profile universally |
| Adapter handling | ENCODE includes explicit adapter handling; the maize study includes trimming | Only explicitly declared synthetic adapter-free inputs accepted; no reads modified | Real plant inputs require adapter evidence and a separately validated trimming policy |
| TSS score | ENCODE uses reference-aware enrichment; plant studies use annotated gene-centered accessibility | Deliberately small ±100 bp toy window, ±10 bp center and 20 bp flanks; undefined background returns null | Not numerically interchangeable with ENCODE or the benchmark's ±2 kb score |

## What this integration does not decide

No human FRiP/TSS cutoff becomes a plant cutoff. MAPQ 30 in the deterministic
validation is an explicit experiment setting, not a newly inferred biological
optimum. Coordinate duplicates do not establish PCR origin. Duplicate marking
retains records in the intermediate BAM; the subsequent analysis exclusion is a
separate, documented policy.

The implementation groups reads in memory and uses a synthetic library identity.
It does not support production-scale multi-library pooling, optical-duplicate
inference, adapter removal, blacklist selection, biological replicate aggregation,
IDR, sc/snATAC barcodes or a universal quality score. Those require separate
designs, fixtures and representative biological data. Local immutable QC image
IDs also need a portable deployment strategy before use on another host.

## Reproducible evaluation

`pixi run checks` runs known-answer arithmetic/filter tests and invalid-policy
rejection. `pixi run validate-atac-legacy --keep` runs the pinned real tools on the
deterministic benchmark fixture, inspects published outputs, verifies raw input
checksums, and requires all eight tasks to cache on resume. The benchmark
manifest and image lock preserve every explicit experimental setting.
