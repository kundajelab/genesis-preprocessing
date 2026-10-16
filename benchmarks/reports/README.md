# Benchmark PDF reports

[Download the seven reports](../../validation/reports/benchmark-pdfs/20261016/README.md)
for package selection, reference choices and the limits of the current evidence.
Each report contains a sample table, metric definitions, measured results and a
scoped recommendation. The reports summarize completed results; building them
does not execute a pipeline or change its defaults.

## Build the reports

Use the dedicated, locked Linux reporting environment. It is independent of the
Genesis processing environments. You also need `pdflatex` with the packages in
[style.tex](style.tex), Ghostscript, and genuine Arial font files. The build rejects
font substitution. Font files are local dependencies and are not committed.

From the repository root:

```bash
mkdir -p workspace/benchmark-reports/tmp
export TMPDIR="$PWD/workspace/benchmark-reports/tmp"
export PIXI_CACHE_DIR="$PWD/workspace/benchmark-reports/pixi-cache"
pixi install --locked --manifest-path benchmarks/reports/pixi.toml
pixi run --manifest-path benchmarks/reports/pixi.toml build \
  --out "$PWD/workspace/benchmark-reports/run-01" \
  --font-dir /path/to/local/Arial-directory \
  --tex-bin /path/to/local/texlive/bin/x86_64-linux
```

Keep those local tool and font directories inside your Genesis workspace. Omit
`--tex-bin` if the prepared `pdflatex` is on PATH. The output directory must not
already exist. The build writes its font cache, TeX cache and temporary files
under that new directory, verifies source hashes, compiles each report twice,
checks page counts and TeX warnings, checks embedded Arial, and renders every
page. Review the rendered PNGs before distributing a new snapshot.

Outputs include seven PDFs, their LaTeX sources, vector figure PDFs, figure data,
sample identities and actual depths, counting results, timing observations,
source hashes, and a build manifest. Logs and rendered previews remain local;
they can contain local tool paths. `build-manifest.json` records their hashes but
does not imply they were published or visually reviewed.

### Tested typesetting dependencies

The October 16 snapshot used TinyTeX-1 v2026.10 (TeX Live 2026), with caption
revision 79618, microtype revision 78228 and TeX Live infrastructure revision
80630. [toolchain.json](toolchain.json) records the release URL, archive digest,
package-file digests and font identities. The Conda `texlive-core` package alone
did not supply a working LaTeX format and is not a dependency of this build.

A local [TinyTeX release](https://github.com/rstudio/tinytex-releases/releases/tag/v2026.10)
can be extracted under your workspace without system installation or PATH
modifications. Additional packages can be installed with that installation's
`tlmgr`; compare their revisions and hashes with the recorded toolchain for an
exact snapshot rebuild. No TeX or font downloads occur during report generation.

The tested Arial files came from the original
[Microsoft Core Fonts archive](https://corefonts.sourceforge.net/), extracted
with `cabextract` without executing the installer. Retain its license locally.
The build checks the font family and embedding flags and records each font hash.
Fonts are embedded as subsets in figures; the original font files are not part
of this repository.

## Evidence and interpretation

[sources.json](sources.json) pins the inputs and their source revisions. Six
completed public-data archives establish **64 method/repetition cases**. The
count is reconstructed from completed, non-resumed processing records; it is
not a count of biological experiments. Archive hashes are checked against the
original publication manifests. The reference assessment discovered in local
storage is a separate compact extraction in [local-audit.json](local-audit.json),
with original result, experiment, mapping and script hashes. Its underlying local
BAMs and full records are not redistributed here.

The reports preserve distinctions between:

- Input depth caps and actual counts, including low-yield treatments and controls.
- Biological replicates, different TF/protocol treatments and technical reruns.
- Exact arithmetic, baseline agreement, independent motifs and biological validation.
- Whole-cohort processing time, individual task observations and indexing.
- Released MACS3 and its pinned optimized source, even though both report 3.0.5.
- Missing measurements (`nd`) and observed zero peaks or zero mismatches.

Whole-cohort timings exclude indexing and common evaluation. The separate
optimized-MACS3 timings use different CPU settings and are not pooled with the
method matrix. Sampled process-group RSS is not whole-host peak memory.

No matched standalone Python-versus-bedtools speed series was established.
No ranking of QC or signal-generation packages is inferred from their presence
in the pipeline. Annotation sensitivity and eligible-coordinate agreement do
not establish reference correctness or model-training benefit.

## Update a snapshot

Review new completed evidence first. Do not change the source catalog merely to
make a failed hash check pass. `freeze_sources.py` deliberately refuses to
replace an existing catalog. For a reviewed new snapshot, update the explicit
source selection, preserve the old snapshot, then freeze new hashes and review
the wording, data and page layouts together. `audit_local.py` can create a new
compact reference audit from an explicitly supplied local result directory;
it never launches missing analyses or follows historical absolute paths.

Run the focused tests with workspace-local caches:

```bash
export MPLCONFIGDIR="$PWD/workspace/benchmark-reports/matplotlib-test-cache"
pixi run --manifest-path benchmarks/reports/pixi.toml test
```

Tests cover provenance failures, duplicate observations, metric incompatibility,
actual depths, missing values, mapping exclusions, arithmetic agreement, absent
Arial and report scope. The published snapshot's validation record states the
exact checks and visual review performed. Scientific pipelines are not rerun by
these report tests.
