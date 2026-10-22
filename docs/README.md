# Genesis documentation

If you are new to Genesis, start with [installation and a first analysis](../README.md#install-and-check-your-environment).
Choose the guide for the task you want to complete. For the scientific question,
measurement definitions and evidence boundaries, read the
[scientific overview](scientific-overview.md).

## Process sequencing data

- [DAP-seq](usage/dapseq-reference.md): sample sheets, controls, references, execution and scientific settings.
- [MACS3 builds and resources](usage/macs3.md): select the optimized source build and inspect caller provenance.
- [Bulk ATAC-seq](usage/bulk-atac.md): paired lanes, adapter/duplicate policies, fragments, peaks and QC.

## Related repository

[genesis-curation](https://github.com/kundajelab/genesis-curation) has study
management, metadata curation and review, single-cell ATAC and the method benchmarks.

## Understand the methods

- [Benchmarks](https://github.com/kundajelab/genesis-curation/blob/main/benchmarks/README.md): datasets, method comparisons, metric definitions and measured results.
- [Validation records](../validation/README.md): dated tests and their limitations.
- [Historical DAP documentation](history/dapseq-readme-fb60096.md): commands and behavior at revision `fb60096`.

## Develop Genesis

- [Contributor guide](development.md): install a development environment, run checks and test changes.
- [Nextflow conventions](nextflow-style.md): module structure and pipeline changes.
- [Diagram maintenance](diagrams/README.md): update the maintained workflow maps.

Current guides describe the installed commands. Historical reports describe the
revision and date written in each report; their commands and temporary paths may
differ from a current installation.
