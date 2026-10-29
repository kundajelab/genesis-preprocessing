# Scientific purpose and evidence

The goal is reliable plant sequence-to-function prediction and experimentally
testable sequence design. Genesis supplies traceable preprocessing and curated
measurements. It does not itself establish a plant foundation model, train the
models discussed below, or provide an MPRA processing/experimental platform.

A useful handoff connects **data → declared measurement → model input → downstream
test**. File readability is one step in that chain. A successful pipeline, a
checksum match and a successful model loader answer different questions.
Nextflow executes supported scientific workflows.

## What does a signal measure?

| Readout | Implemented measurement | Current model interface / limitation |
| --- | --- | --- |
| DAP binding | Retained mapped primary read ends; total/strand CPM coverage and strand-specific 5′ counts. Peak RPM counts read ends; mean RPKM is coverage. Paired ends count independently. | Peaks and tracks are available; the scATAC model adapters are not DAP adapters. In-vitro TF binding is not tissue-specific transcription or accessibility. |
| Bulk ATAC accessibility | Usable paired outer fragments in BED0 half-open coordinates; raw one-base cuts at `start + 4` and `end − 6`. FRiP counts usable reads in peaks, as ENCODE does. | Fragment, peak and raw-cut outputs; no automatic conversion through the scATAC adapters. Declare and audit the model convention before reuse. |
| scATAC accessibility | Existing 10x fragments, already offset +4/−5; selected fragment rows count once, independent of supporting read-pair count. Biological replicates stay separate. | The audited `chrombpnet-atac-v1` conversion puts cuts at `s` and `e` for source `[s,e)`. One unstranded raw-count channel, no normalization or further shift; separate pinned ChromBPNet/Cherimoya adapters. |
| TSS / RNA / WGBS / MPRA | Contextual or future prediction targets, not implemented assay workflows here. Annotation-derived TSS QC in bulk ATAC is not TSS-assay processing. | Initiation, expression, methylation and reporter activity require distinct measurement definitions and acceptance criteria. |

Identical bigWig or BED formats do not make targets interchangeable. Preserve count
unit, strand, normalization, coordinates and shift provenance. The
[DAP](usage/dapseq-reference.md#analysis-choices), [bulk
ATAC](usage/bulk-atac.md#analysis-and-metrics) and
[scATAC](https://github.com/kundajelab/genesis-curation/blob/main/docs/usage/scatac.md#process-known-libraries)
guides specify the actual policies.

Sequence-source organism, assayed DNA organism, TF-source organism and experimental
host/context are distinct concepts. This matters for interspecies DAP and for
sequences tested in a heterologous reporter host. Where these distinctions are
absent from an existing manifest, they are future metadata/review requirements;
do not invent values or interpret a single species field as all four.

## Implemented does not mean demonstrated at every scope

| Capability or claim | Evidence scope | Remaining acceptance |
| --- | --- | --- |
| Contracts, fragment ingestion, replicate grouping and inspection | Synthetic software fixtures | Real-data scientific suitability |
| Scientific pipeline integration | Actual DAP/ATAC tools on synthetic inputs | Representative biological cohorts; synthetic peaks are not biological truth |
| Public-data component comparisons | Six Arabidopsis matrices, 1M/5M template caps; seeded nested sampling, seed 60141 | Broader species, full-depth tests and independent accuracy evidence |
| Full real-data analysis | Historical full fragment scans exist for the Arabidopsis scATAC candidate; they are not a full pseudobulk/model release | Re-establish accessible inputs, labels and reference/annotation evidence; reviewed release absent |
| A model loader reads data | Both pinned loaders read synthetic positive/background windows in three folds | Exact-bundle, exact-model real-data loading; no transfer of fixture PASS |
| Model training and independent experimental validation | Not demonstrated by these results | Training evaluation and independent biological/experimental tests |

The six archived `results.json` records reconcile to **64 method/repetition
cases**: SE DAP 12 + 4, PE DAP 12 + 4, bulk ATAC 24 + 8 at the 1M/5M caps.
Their stored result hashes and the 64-row timing table agree. This is inspection of completed records, not a new pipeline rerun.
Technical repetitions and depth/method cases are not extra biological experiments.
The sampling selects lowest seeded hashes of record index and template name from
verified inputs, preserves mates and nests the lower cap inside the higher cap;
actual retained/input denominators are recorded rather than inferred from a cap.

## Three observations supported by completed evidence

| Observation | Representative result and denominator | Interpretation and limits |
| --- | --- | --- |
| No general replacement for current defaults is established | At the 5M SE DAP cap, FUS3 motif hits are 466/998 matched peak windows for baseline, 465/998 for Bowtie2 and 468/998 for MACS2. Two selected FUS3 windows lack matched controls. | Similar independent motif support does not rank biological accuracy. Picard is a performance candidate: 243.2 s versus baseline 289.3 s median for the two-library 1M ATAC cohort, over three technical repetitions. Larger-genome acceptance remains separate. |
| Greater depth improves peak discovery and replicate agreement in this cohort | Two Arabidopsis ATAC biological replicates: peaks rise 4,997→14,805 and 3,043→11,551 from 1M→5M input pairs each. Peak-base-pair Jaccard rises 0.2404→0.4858; cross-replicate FRiP rises 0.143339→0.460099 and 0.163293→0.402303. | Peak coverage also increases; higher held-out FRiP alone is insufficient. TSS ratios are comparatively stable. This is one biological pair, not a universal depth rule. |
| Annotation/reference choices change reported measurements | On the same Col-CEN alignments, all-transcript Araport11/TAIR10 TSS sets contain 46,558/32,986 eligible positions. Ratios are 3.6894/3.8027 and 2.7517/3.2096 for the two libraries. | These are different denominators. For the 18,428 matched eligible single-start genes, positions can still differ between annotation releases. A higher ratio is not proof of better annotation. |

These results come from the 2026-10-14 method comparison. Earlier validation
reports and evidence were removed from this repository; see `validation/` at commit
e2e2aa9.

Cross-assembly agreement applies only after reciprocal coordinate conversion:
for the first ATAC library, 4,992/4,997 TAIR10 and 5,586/5,596 Col-CEN peaks pass;
for the second, 3,040/3,043 and 2,968/2,977 pass. Failed/ambiguous mappings stay
excluded and counted. Matched-region agreement is not ground truth. Keep those
exclusions and mapping assumptions with any comparison.

## What this evidence does not yet establish

These comparisons do not show that foundation-model pretraining improves plant
sequence design. Loader IO is not training, training is not independent biological
validation, and motif plausibility is not experimental confirmation. Neither more
peaks nor higher own-FRiP/mapping rate is a sufficient method-selection criterion.

Chromosome folds do not establish absence of homology leakage. The current model
backgrounds are unmatched genome tiles. ChromBPNet's bias-model plan does not
supply a fitted plant bias model. scATAC TSS enrichment remains unassessed, and
absence of organellar fragments in a filtered input does not mean zero original
organellar contamination. No universal plant QC cutoff is implied.

## The next data-to-model decision

A bounded first acceptance exercise is one reviewed Arabidopsis root cell type,
with separate bundles for its two biological replicates, targeting the audited
ChromBPNet revision. Select the cell type, resolve reference/annotation evidence,
agree on raw-cut and background/bias acceptance, then test each exact bundle with
the model loader. This is a proposed handoff, not a released or approved dataset.

Use [model
inspection](https://github.com/kundajelab/genesis-curation/blob/main/docs/usage/scatac-model-handoff.md#inspect-a-bundle)
to distinguish observed structure, declarations and missing acceptance. The modeling
team still needs to define the held-out prediction task and evidence required before
training.