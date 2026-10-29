# Genesis documentation

If you are new to Genesis, start with [installation and a first analysis](../README.md#install-and-check-your-environment).
Choose the guide for the task you want to complete. For the scientific question,
measurement definitions and evidence boundaries, read the
[scientific overview](scientific-overview.md).

## Process sequencing data

- [DAP-seq](usage/dapseq-reference.md): sample sheets, controls, references, execution and scientific settings.
- [MACS3 builds and resources](usage/macs3.md): select the optimized source build and inspect caller provenance.
- [Bulk ATAC-seq](usage/bulk-atac.md): sample sheets, trimming, filtering, fragments, peaks and QC.
- [scATAC pseudobulks](usage/scatac-seq.md): fragment files, cell lists, pseudobulk peaks, tracks and QC.
- [RNA-seq](usage/rna-seq.md): plant annotations, STAR intron limits, gene counts, tracks and QC.
- [TSS mapping](usage/tss-seq.md): CAGE, TSS-seq and csRNA-seq trimming, cap-G correction, TSS tags, paraclu clusters and QC.

## Related repository

[genesis-curation](https://github.com/kundajelab/genesis-curation) has study
management, metadata curation and review, single-cell ATAC model export and the
method benchmarks.

## Understand the methods

- [Benchmarks](https://github.com/kundajelab/genesis-curation/blob/main/benchmarks/README.md): datasets, method comparisons, metric definitions and measured results.
- [Historical DAP documentation](history/dapseq-readme-fb60096.md): commands and behavior at revision `fb60096`.

## Develop Genesis

- [Contributor guide](development.md): install a development environment, run checks and test changes.
- [Nextflow conventions](nextflow-style.md): module structure and pipeline changes.
- [Diagram maintenance](diagrams/README.md): update the maintained workflow maps.

Current guides describe the installed commands. The historical DAP documentation
describes its own revision; its commands may differ from a current installation.
Earlier validation reports and evidence were removed; see `validation/` at commit
e2e2aa9.
