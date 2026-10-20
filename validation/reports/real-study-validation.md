# Real DAP/ATAC study validation

> Validation record for the revision and date below. Use the
> [current task guides](../../docs/README.md) for a new analysis. Recorded temporary
> paths identify the original test environment.

The complete prepared-study workflow passed real scientific processing, full artifact
validation and the default reviewed export on seven synthetic libraries. Validation ran
on **2026-10-13** from clean revision
`311567b7c327ad1fc9b177f193b4b5e1951e02e5`, extending baseline
`36766e8887c2026030a4b0976c402ac834eded23`. Documentation/evidence follow-ups preserve
that tested implementation. These fixtures establish software integration, not biological
quality or acceptance on another server.

## Delivered behavior

> **Later change:** a subsequent revision removed `dap-run.json` and its per-task
> receipts. DAP runs now record `dap-images.lock` and `git-status.txt` instead; see
> [DAP run records](../../docs/usage/dapseq-reference.md#run-records). This report
> describes the revision it validated.

New DAP runs publish `output/provenance/dap-run.json` through one final Nextflow task.
Its shipped strict version-1 schema binds input checksums, sample/control relationships,
reference/cache identity, contigs, effective scientific parameters, code/runtime identity,
actual task containers, requested resources, script hashes and published artifact hashes.
Analysis-task receipts survive resume. Reported application versions remain in the bound
native reports. Resource observations are explicitly null in receipts; raw Nextflow traces
retain measured task runtime/RSS. No network-byte or whole-host memory estimate is inferred.

Curation derives completeness from evidence. Legacy runs without a sidecar remain readable
and incomplete; malformed present sidecars are rejected. Full validation and human review
remain separate requirements. A reviewer cannot waive corrupt files. No scientific filter,
control assignment, reference choice or plant QC threshold changed.

The public `pixi run validate-study` acceptance runner uses actual scientific containers,
two local worker identities, separately staged inputs and separate worker roots. Only its
newly created synthetic registries can receive fixture approvals. It never accepts an existing
production study as a target. The controller, worker transport and scheduler remain the existing
implementations; there is no second execution framework.

## Commands and results

The recorded host used Linux x86_64, Python 3.14.5, Nextflow 25.10.4, the existing
pinned images, and a separate pyBigWig 0.3.26 interpreter. [Environment/image identities](../evidence/real-study/environment.json)
and the individual JSON receipts record exact paths, hashes and settings. The verified plugin
cache was selected with `NXF_PLUGINS_DIR=/tmp/genesis-nextflow-plugins-v1` and
`NXF_OFFLINE=true`; `UV_LINK_MODE=copy` was explicit. Paths identify this run, not prerequisites
on another host. Some acceptance commands ran concurrently; durations are observations, not
capacity or latency guarantees.

| Check | Result | Observed duration |
| --- | --- | --- |
| `pixi run checks` | **PASS**: lint/format/types, scientific known answers, actual Nextflow stubs, execution, contracts/registry/review, mocks, docs and diagrams | 194.217 s |
| `pixi run validate-study --keep --work-root /var/tmp/genesis-study-MQ6B5glk --bigwig-python /tmp/genesis-bigwig-env/bin/python` | **PASS**: 145 public CLI calls, 5 DAP + 2 ATAC libraries, full validation/default export and recovery | 256.034 s |
| `pixi run validate-docker --keep` | **PASS**: full DAP tools, references, inputs, QC, peaks, quantification, publication, resume and focused FastQC/MultiQC cases | 370.511 s |
| `pixi run validate-atac --keep` | **PASS**: two libraries/two lanes each, repeated resume, helper-code invalidation, unchanged science and adapter case | 93.668 s |
| Build/install wheel outside checkout | **PASS**: all bundled resources, entry points, offline review workflow and strict DAP provenance checks; `PYTHONPATH` unset | 13.711 s |
| Earlier-release scientific comparison | **PASS**: 59 matching scientific products against `36766e8`, using the same DAP input bytes | See [hashes and method](../evidence/real-study/scientific-invariance.json) |
| Final documentation and diagram checks | **PASS**: 197 CLI examples, 241 Pixi recipes, links/anchors and archived README; all seven generated maps reproduce | Recorded with final docs |

Receipts and logs: [full checks](../evidence/real-study/full-checks.json),
[study](../evidence/real-study/real-study.json), [DAP](../evidence/real-study/real-dap.json),
[ATAC](../evidence/real-study/real-atac.json), [installed wheel](../evidence/real-study/installed-wheel.json).
The [evidence index](../evidence/real-study/README.md) distinguishes compact committed records
from full temporary work products. No unavailable test is presented as passed.

## What the study run established

DAP used the existing 75-base SE/PE fixture: three treatments and two controls, including
one control shared by two SE treatments. The two analysis groups completed **25 and 20 tasks**.
ATAC used the existing two-library fixture with two technical lanes per library; each worker
completed **11 tasks**. All four analysis jobs remained at **attempt 1** throughout collection,
revalidation, corruption and repair. [Job counts, export summaries and provenance](../evidence/real-study/study-summary.json)
and [all CLI calls](../evidence/real-study/study-commands.json) record the transitions.

For both assays, the runner checked:

1. Prepared registration, source-backed synthetic input review, immutable planning and actual
   Nextflow execution. A separate controller process reconnected to the existing supervisors.
2. Metadata-only collection excluded all candidates from default export. Full byte/format scans
   passed before explicit QC/eligibility decisions included **5/5 DAP** and **2/2 ATAC** candidates.
3. Repeated unchanged collection imported no additional receipt and preserved reviewed export.
4. A one-byte scan budget produced incomplete validation and excluded all candidates. Revalidation
   made earlier QC/eligibility decisions stale; fresh decisions were required.
5. Corrupting one synthetic bigWig produced an error while the independent analysis group remained
   structurally valid. No damaged candidate was exported. Restoring bytes and recollecting repaired
   validation without launching another scientific attempt.
6. An invalid run sidecar caused a confirmed collection failure while the other group collected.
   Restoring it and retrying collection recovered both final reviewed exports.

The standalone DAP fixture completed **43 initial tasks**. The first resume completed 29 and
cached 10 as reference paths moved to published caches; the stable resume completed only the
uncached sheet/reference validator and cached **38 tasks**, including the provenance task.
[Trace summary](../evidence/real-study/dap-summary.json). All three treatments still yielded
200 peaks. No FASTQs, BAMs, indexes or bedGraphs entered published scientific outputs.

The explicit baseline comparison checked 35 bigWigs, 3 decompressed narrowPeak files,
6 quantification TSVs, 5 metadata TSVs, and 10 flagstat/idxstats files. All 59 hashes match.
Report timestamps, provenance and presentation are outside that scientific comparison.

## Failures found and fixed

The first real run exposed a campaign compiler bug: canonical JSON reordered DAP's fixed
six-column sheet header. The compiler now uses the existing `SAMPLE_HEADER`; a regression
passes its generated sheet through the production sample validator.

Recovery after restoring identical bytes initially hit the registry's historical-import
protection. Each collection generation now adds its exact request to its immutable source
evidence, so a repaired observation becomes a new version while historical-import protection
remains intact. Existing input/review binding checks still govern promotion.

A standalone FastQC regression exposed the first receipt helper's reliance on a caller's
`lib/` directory. The helper is now an explicitly imported Nextflow module resolved relative
to each scientific module. Standalone FastQC and MultiQC cases pass. The generated-map check
also caught a stale process map; source and static/animated SVGs now include DAP_PROVENANCE.
[Initial failure logs](../evidence/real-study/README.md) are retained as failures, not counted
as acceptance passes. The existing Pixi `system-requirements` deprecation warning remains.

## Reproduce and continue

Use a clean checkout, existing approved Docker access, the pinned images and a separately
installed bigWig reader. See [the complete study recipe](../../docs/usage/studies.md#real-tool-acceptance-on-synthetic-inputs)
for input staging, fixture-only review and Docker-visible paths:

```bash
pixi run install-all
pixi run checks
pixi run validate-study --keep --work-root /var/tmp \
  --bigwig-python /path/to/bigwig-environment/bin/python
pixi run validate-docker --keep
pixi run validate-atac --keep
```

New reference recipe hashes require verified fresh derived caches: preserve old caches and
regenerate in a new directory containing the same approved FASTAs. See
[DAP run records](../../docs/usage/dapseq-reference.md#run-records) (current behavior).
Incomplete conda/native or mutable-runtime identities are not upgraded to complete evidence.
Source snapshots and receipts remain bounded to 32 MiB; this tiny acceptance does not establish
large-study throughput or catalog scale.

**NOT RUN:** live SSH, schedulers, API inference, GPU/vLLM deployments and public biological
cohorts. **NOT EVALUATED:** biological metadata accuracy or universal plant QC thresholds.
The next acceptance step is an explicitly configured SSH worker on the destination site,
then a small, reviewed public dataset. No credentials, production policies, infrastructure
or production datasets changed.

Rollback uses reviewable revert commits, preserving public history. Stop affected orchestration
before changing code, retain immutable run/attempt/evidence directories, and restore a consistent
registry backup into a new location when needed. This milestone adds no registry migration;
older code will not recognize the new DAP evidence. Never erase approval history or mutate a
campaign to simulate a rollback.

## Installation from a fresh checkout

A fresh checkout of `f2a1a731f17ba25b0e82dd5b56137beaf29e957a` installed new
environments and passed CLI help/doctor, all 197 documented examples, the offline
tutorial, strict DAP provenance checks and local study acceptance with the
separate binary reader. All **574 tracked files** matched the source revision,
and the checkout remained clean. Elapsed time was **47.989 seconds**, using
existing package/plugin caches on this host. This does not establish compatibility
on other hosts. The harness clears `PYTHONPATH` and the parent `VIRTUAL_ENV` so
the checkout uses its own environment.

The [receipt](../evidence/real-study/fresh-clone.json),
[per-command results](../evidence/real-study/fresh-clone-details.json),
[complete log](../evidence/real-study/fresh-clone.log) and
[documentation checks](../evidence/real-study/final-docs.json) preserve the scope.
See the [contributor guide](../../docs/development.md) for installation instructions.

All 290 study command stdout/stderr files are retained in a
[compact archive](../evidence/real-study/study-cli-logs.tar.gz) as sanitized publication
copies. The [publication manifest](../evidence/real-study/publication-manifest.json)
binds each original per-command hash to its published member hash. Public review-token
values and command arguments are redacted; scientific measurements are unchanged.
Owned temporary checkouts and large scientific work products were
cleaned up after this evidence was saved. Reproduce those products with the documented
fixture commands; the historical temporary paths are not permanent data locations.
