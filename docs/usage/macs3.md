# MACS3 builds and peak-calling resources

DAP-seq peak calling builds MACS3 from upstream source. The DAP environment,
`environments/DAP_SEQ_PEAKS.yaml`, and its image pin revision
`da19e1dcbf3a2ff06ab906abc0658edb0c7e82d1`. It includes the callpeak
optimization, later fixes to worker cleanup, process safety and concurrent
scoring, and a CPU budget (`--threads`) that bounds MACS3's threads and worker
processes. The package still reports `macs3 3.0.5`; the pinned revision in the
environment lock and `dap-images.lock` distinguishes it from the released 3.0.5
package.

## Profiles and images

The DAP profiles (`local,conda`, `local,docker`, `sherlock,apptainer`) use this
environment or its image. Until a rebuilt image is published,
`conf/dap-images.lock` may still name an earlier peak-calling image for container
profiles. Each run copies that lock to its output folder as `dap-images.lock`
(see [DAP run records](dapseq-reference.md#run-records)).

Source builds require a C compiler, zlib headers and the pinned SIMDe Git
submodule; the environment declares these build dependencies. Do not substitute
a GitHub source archive: it omits the submodule. Native macOS and ARM builds
require their own acceptance tests; resolving their lock entries does not
validate them.

The bulk ATAC and scATAC workflows pin the released peak-calling image separately
and are unchanged. The scATAC workflow ingests existing fragments and calls peaks
from BEDPE pseudobulks; it does not process scATAC FASTQs or switch to MACS3's
FRAG mode.

## CPU limits

DAP peak calling runs `macs3 callpeak` directly with `--threads` set to the
task's CPU allocation. The budget includes MACS3's calling thread and bounds its
threads and peak worker processes on any operating system. MACS3 starts peak
worker processes only on Linux; elsewhere those phases run serially. A container
or scheduler CPU limit remains the outer bound.

The DAP peak-calling process requests four CPUs, 1 GB plus a quarter of the
treatment and control BAM sizes, and one hour, retrying twice with more memory and
time if killed. On a Sorghum sample (52M treatment and 12M control read pairs,
8.1 GB of BAMs) with revision `da19e1dcbf`, four CPUs called peaks in 13 s against
73 s for one, eight CPUs were no faster than four, and peaks were byte-identical at
every budget. That sample needed about 1.1 GB at one CPU and ran within 2 GB at
four. Use explicit Nextflow process overrides to change resources; faster
execution can use more memory. The caller's thresholds, controls, genome size,
duplicate policy and input formats are unchanged.

Bulk ATAC and scATAC use the same pinned MACS3 image and call `macs3 callpeak`
directly with `--threads` set to the task's CPU allocation (`${task.cpus}` for
bulk ATAC; two CPUs and 8 GB for scATAC products). On the archived A. thaliana
ATAC inputs (two libraries at 1M and 5M templates), this build gave peak, summit
and XLS outputs identical to released MACS3 3.0.5 at two and four CPUs.

`genesis_tools/macs3_launch.py` remains only for the benchmark harness, which runs
released MACS3 3.0.5 natively; that release has no `--threads` option, so the
launcher bounds it to the benchmark's CPU allocation.

## Earlier measurements

The [measured results](../../validation/reports/macs3-optimization.md) and the
[comparison recipe](../../benchmarks/experiments/macs3-optimization/README.md)
compare the earlier revision `22eb40f7b51f5b20b635614ca0b1759dc5398c5c` with the
released package, using the same inputs and scientific arguments, 1, 2 and 8
CPUs and three measured repetitions. They describe that revision, not the current
DAP pin. Benchmark environments and results remain version-specific; do not
overwrite results from one caller build with those from another.
