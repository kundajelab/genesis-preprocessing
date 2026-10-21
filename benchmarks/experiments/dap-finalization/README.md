# Pilot-first DAP-seq method comparisons

These tools run DAP-seq method comparisons in small, checked steps before any
larger run. A comparison never changes pipeline defaults or approves a dataset.

## Pilot before expanding

Write down the question, the one changed variable, the fixed settings and the
metrics before running anything. Run every candidate on 50,000 templates sampled
deterministically (seed 60141, mates kept together) with its matched control and
the intended full reference.

`pilot.py` supports that sequence without modifying results:

- `inspect()` re-hashes every recorded artifact, requires exact agreement with
  independent counts, and lists tool warnings that have not been reviewed. A
  zero-peak result is not a malformed output; a caller that cannot fit its model
  has failed the test, which is not a reason to change its settings.
- `forecast()` estimates a larger run as twice the measured setup time plus the
  slowest measured per-template rate. It is an allowance, not a guarantee. A
  200,000-template pilot is required when the estimate exceeds 30 minutes.
- `admit()` starts a larger run only if it still leaves three hours for
  reporting and three for contingencies before the deadline.
- `verify_binding()` rejects a pilot whose specification or dataset changed.

Only pilots with intact evidence and no unreviewed warnings proceed. Record
failed pilots and comparisons excluded before expansion in the report.

The reference-comparison queue accepts `--max-hours` (1–168, default 168). Its
start time is stored, so restarting it cannot extend a shorter deadline. An
interrupted job is recorded as `INTERRUPTED` rather than left running.

## Reusing indexes and alignments

`run_native.py` runs a single-repetition comparison and reuses reference
indexes and alignments from earlier comparisons. Reuse is keyed on the recipe,
this wrapper, the environment, the CPU allocation, the reference record, the
read checksums, the library layout and the alignment settings. Each shared file
is re-hashed before reuse and after the run. Reuse is recorded in `reuse.json`;
elapsed times from such runs are not cold whole-workflow timings. Pipeline
execution and its own caches are unaffected.

`check_reuse.py` is an opt-in real-tool check that needs the benchmark Pixi
environment: a second run must reuse every index and alignment and reproduce
every metric.

## Method options

DAP comparisons can include adapter-only trimming with Cutadapt or fastp, for
single-end as well as paired-end libraries. Declare both adapter fields; the
second applies only to paired-end reads. fastp's quality, length and poly-G
filters are disabled, but its paired-end overlap-based adapter detection
remains part of the method. An exploratory adapter sequence does not establish
a trimming policy.

Bowtie2 methods can set `"bowtie2_mode"` to `"local"` or `"end-to-end"`.
Without it, Bowtie2 keeps `--very-sensitive-local`.

## Recorded results (October 2026)

These runs used a frozen source snapshot without Git metadata, identified by
per-file recipe hashes, and the native benchmark environment with released MACS3
3.0.5. They are not runs of the Nextflow DAP pipeline or its containers, and
they say nothing about MACS3 builds pinned later. Inputs were the nested
seed-60141 samples, identified by SHA-256.

| Comparison | Largest depth | Outcome |
| --- | --- | --- |
| A. thaliana SE adapter trimming (none, Cutadapt, fastp), TAIR10 | 1M | Complete, cold run. MACS3 warned that its 140–141-base fragment model is shorter than twice the read length; reviewed, kept as a limitation, with no settings changed |
| A. thaliana PE adapter trimming | 1M | Complete, cold run |
| PE MAPQ 0, 10 and 30 with duplicates excluded | 1M | Complete for one representative treatment and control, not the full PE panel. All three reused the reference index; MAPQ 10 and 30 also reused both alignments. These timings are not cold |
| A. lyrata adapter trimming | 200k | Complete; not expanded because the conservative estimate exceeded the remaining budget |
| Sorghum Sorbi1 reference | 200k | Complete; not expanded for the same reason |
| PE duplicate markers (samtools, Picard) | 50k | Needs correction: Picard could not parse SRA read names, so optical-duplicate detection was inactive |
| Sorghum BTx623 v5 reference | 50k | Needs correction: a samtools reference-naming warning is unreviewed |
| A. lyrata and sorghum aligners; sorghum trimming | 50k | Stopped at the compute allowance; no result |
| A. lyrata aligners (earlier 1M queue) | — | Failed at the campaign wall-time limit |
| Plate-pooled control designs | — | Not run |

The completed comparisons describe how each choice changes retained reads,
peaks and counts. They do not establish that any alternative is better, and
they do not justify changing alignment, trimming, duplicate, reference or
peak-calling defaults.
