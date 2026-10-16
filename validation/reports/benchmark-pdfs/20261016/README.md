# Plant preprocessing benchmark reports

These short reports compare processing choices using completed Arabidopsis
DAP-seq and bulk ATAC results. Each defines its measurements and distinguishes
software agreement, biological evidence and computational cost. No pipeline
defaults were changed to produce this collection.

| Report | Pages | Main comparison |
| --- | ---: | --- |
| [Evidence overview](overview.pdf) | 2 | Supported decisions and remaining coverage gaps |
| [Alignment](alignment.pdf) | 2 | BWA-MEM2, BWA and Bowtie2 |
| [Adapter trimming](trimming.pdf) | 2 | Untrimmed input, Cutadapt and fastp |
| [Duplicate handling](duplicates.pdf) | 2 | samtools and Picard |
| [Peak calling](peak-calling.pdf) | 3 | MACS2, MACS3, Genrich and the pinned MACS3 optimization |
| [Interval counting](interval-counting.pdf) | 2 | Genesis/Python and bedtools arithmetic |
| [References and annotations](references.pdf) | 3 | TAIR10/Col-CEN coordinate agreement and TSS annotation sensitivity |

The six public-data subset archives reconcile to **64 method/repetition cases**,
not 64 biological experiments. The optimized-MACS3 experiment is separate:
8 input cases, 3 CPU allocations, 2 builds and 4 invocations per combination
(one warmup and three measured runs), totaling 192 calls. All 24 paired
comparisons preserve the evaluated scientific peak outputs.

The reference report includes a later completed 5-million-pair assessment that
was still pending in the earlier public narrative. Its source receipt and
coordinate exclusions are retained in the [compact local audit](data/local-audit.json).
That extraction is recorded evidence, not an independently rerun analysis.
It does not establish completion of the broader species panel.

## Data and reproducibility

- [Source catalog and revisions](data/sources.json)
- [Sample identities, controls and actual depths](data/samples.tsv)
- [Completed method/repetition cases](data/cases.tsv)
- [Independent counting results](data/counting.tsv)
- [Individual task resource observations](data/resources.tsv)
- [Build manifest and artifact hashes](build-manifest.json)
- [Validation record](VALIDATION.md)
- [Build instructions and shared LaTeX/plot sources](../../../../benchmarks/reports/README.md)

Each PDF has a matching `.tex` file in this directory. Vector figures and their
numerical tables are under `figures/`. The build manifest records hashes of
local compilation logs and page previews; those potentially machine-specific
files are not distributed. The reports contain no private input locations,
credentials or internal meeting material.

A higher own FRiP, mapping rate or peak count is not a tool-accuracy ranking.
The interval-counting report establishes integer agreement, not a measured
standalone speed advantage. QC/track package rankings, model training benefits
and experimental validation remain outside this evidence.
