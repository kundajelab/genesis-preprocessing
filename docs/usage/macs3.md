# MACS3 builds and peak-calling resources

DAP-seq peak calling builds MACS3 from upstream source. The DAP environment,
`environments/GENESIS_MACS3.yaml`, and its image pin revision
`da19e1dcbf3a2ff06ab906abc0658edb0c7e82d1`. It includes the callpeak
optimization, later fixes to worker cleanup, process safety and concurrent
scoring, and a CPU budget (`--threads`) that bounds MACS3's threads and worker
processes. The package still reports `macs3 3.0.5`; the pinned revision in the
environment lock and `genesis-images.lock` distinguishes it from the released 3.0.5
package.

## Profiles and images

The DAP profiles (`local,conda`, `local,docker`, `sherlock,apptainer`) use this
environment or its image. Until a rebuilt image is published,
`conf/genesis-images.lock` may still name an earlier peak-calling image for container
profiles. Each run copies that lock to its output folder as `genesis-images.lock`
(see [DAP run records](dapseq-reference.md#run-records)).

Source builds require a C compiler, zlib headers and the pinned SIMDe Git
submodule; the environment declares these build dependencies. Do not substitute
a GitHub source archive: it omits the submodule. Native macOS and ARM builds
require their own acceptance tests; resolving their lock entries does not
validate them.

Bulk ATAC uses the same pinned peak-calling image as DAP-seq.

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

Bulk ATAC also calls `macs3 callpeak` directly with `--threads ${task.cpus}`. On
the archived A. thaliana ATAC inputs (two libraries at 1M and 5M templates), this
build gave peak, summit and XLS outputs identical to released MACS3 3.0.5 at two
and four CPUs.

## Earlier measurements

The archived [comparison
recipe](https://github.com/kundajelab/genesis-curation/blob/12bb02cb9ec98fe8c4c0ec660adcfdeb7384f20e/benchmarks/experiments/macs3-optimization/README.md)
compares the earlier revision `22eb40f7b51f5b20b635614ca0b1759dc5398c5c` with the
released package, using the same inputs and scientific arguments, 1, 2 and 8 CPUs and
three measured repetitions. Its results describe that revision, not the current DAP pin.
Benchmark environments and results remain version-specific; do not overwrite results
from one caller build with those from another. The measured results were removed
from this repository; see `validation/reports/macs3-optimization.md` at commit
e2e2aa9.