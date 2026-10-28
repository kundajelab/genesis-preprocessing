# Contributing to Genesis

This guide covers development setup, tests and contribution conventions.
Clone [Genesis](https://github.com/kundajelab/genesis-preprocessing) to obtain the
source, examples, fixtures and dependency locks. No Git LFS objects or submodules
are required.

## Set up a fresh checkout

Install Git and [Pixi](https://pixi.sh/) for your server, then run:

```bash
git clone --branch main https://github.com/kundajelab/genesis-preprocessing.git
cd genesis-preprocessing
git log -1 --oneline
pixi run install-all
pixi run pipeline --help
```

`install-all` uses the committed `pixi.lock`, `genesis_tools/uv.lock` and Python
selection in `genesis_tools/.python-version` (currently Python 3.14.5 with the GIL).
It creates environments for the checkout. Initial installation needs network
access or populated package caches. Recreate `.pixi/` and `.venv/` when moving a
checkout: environment paths and binaries are machine-specific. The supported Pixi
platforms are Linux x86_64 and macOS x86_64/ARM64.

Docker is needed for the default container profile; the checks run without it.

<a id="command-installation-and-maintenance"></a>
Study management, curation and method benchmarks, including the `genesis` command,
are in [genesis-curation](https://github.com/kundajelab/genesis-curation).

Read [AGENTS.md](../AGENTS.md) before editing. Create a branch for the next change:

```bash
git switch -c work/next-change
```

## Establish your baseline

```bash
pixi run checks
git status --short
```

`checks` includes lint, formatting checks, types, DAP and ATAC Nextflow stubs, the
evidence publication check and README command checks. Its initial Nextflow
run needs access to the pinned plugin registry or a verified cache containing
`nf-schema@2.4.2` and `nf-dotenv@1.0.0`. If plugin resolution fails, restore registry
access or configure an approved cache with `NXF_PLUGINS_DIR`; only enable
`NXF_OFFLINE=true` after populating it. The historical `/tmp/genesis-*` paths in
validation reports identify old test runs and are not paths to reuse on a new server.
The [hybrid validation report](../validation/reports/hybrid-llm-validation.md) records
the verified plugin archive hashes.

### Recovery when plugin lookup fails

If `checks` stops with `Plugin with id nf-schema not found in any repository`,
this recipe downloads the exact two pinned archives from the official registry,
checks their SHA256 values, and installs them in an ignored checkout-local cache.
Run it from the repository root. It refuses to replace an existing plugin directory.

```bash
pixi run uv run --frozen --project genesis_tools python - <<'PY'
import hashlib
import io
import tempfile
import zipfile
from pathlib import Path
from urllib.request import urlopen

cache = Path(".nextflow/genesis-plugins").resolve()
cache.mkdir(parents=True, exist_ok=True)
plugins = [
    ("nf-schema", "2.4.2", "02c2fc38cefd5a38238a208ec6ae52927a3f509f07c97f0e266841924cbbfd9d"),
    ("nf-dotenv", "1.0.0", "bb75972f8557f547b3e1472004e6d71a586b1614431437d3ec64590f550735d4"),
]
for name, version, expected in plugins:
    target = cache / f"{name}-{version}"
    if target.exists():
        raise SystemExit(f"Refusing to overwrite {target}")
    url = f"https://registry.nextflow.io/api/v1/plugins/{name}/{version}/download/{name}-{version}.zip"
    with urlopen(url, timeout=60) as response:
        data = response.read(64 * 1024 * 1024 + 1)
    if len(data) > 64 * 1024 * 1024 or hashlib.sha256(data).hexdigest() != expected:
        raise SystemExit(f"Archive verification failed for {name}")
    with tempfile.TemporaryDirectory(prefix=".stage-", dir=cache) as temporary:
        stage = Path(temporary)
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            for item in archive.infolist():
                if item.filename.startswith("/") or ".." in Path(item.filename).parts:
                    raise SystemExit("Unsafe archive member")
            archive.extractall(stage)
        stage.rename(target)
    print(f"Verified {name}@{version}: {expected}")
PY

NXF_PLUGINS_DIR="$PWD/.nextflow/genesis-plugins" NXF_OFFLINE=true pixi run checks
```

Reuse the final command on subsequent runs. If the server cannot reach the download
URLs either, prepare this verified cache on a connected machine and transfer those
two plugin directories through your approved storage workflow. Nextflow documents
the cache and offline settings in its [plugin guide](https://docs.seqera.io/nextflow/plugins/using-plugins).
Do not change the pinned plugin versions or bypass failed checks to recover access.

Pixi reports a `system-requirements` deprecation warning.
If uv's package cache and the checkout are on different filesystems, its hardlink
attempt may fall back to copying; `UV_LINK_MODE=copy pixi run install-all` selects
copying explicitly. Neither warning is a failed scientific check.

## Change and test a feature

Keep scientific settings explicit. A change to alignment, filtering, peak calling or
count definitions needs a comparison of the affected outputs, not just a successful
exit code. Design that comparison with the [benchmark
guide](https://github.com/kundajelab/genesis-curation/blob/main/benchmarks/README.md).

| Area | Implementation | Focused check |
| --- | --- | --- |
| DAP | `workflows/dap-seq.nf`, `conf/dap-seq.config`, `modules/dap-seq/`, `genesis_tools/samples.py` | `tests/verify_pipeline.py` |
| Bulk ATAC | `workflows/atac-seq.nf`, `conf/atac-seq.config`, `modules/atac-seq/` | `tests/verify_atac.py` |
| scATAC | `workflows/scatac-seq.nf`, `conf/scatac-seq.config`, `modules/scatac-seq/` | `tests/verify_scatac.py` |
| RNA-seq | `workflows/rna-seq.nf`, `conf/rna-seq.config`, `modules/rna-seq/` | `tests/verify_rna.py` |
| Shared steps | `nextflow.config`, `modules/common/` (download, references, FastQC, MultiQC report) | all of the above |
| Documentation | README and task guides | `tests/verify_readme.py` |

Python package paths in this table are relative to `genesis_tools/src/`.
Run a focused check with the project interpreter, then run the full suite before
publishing a code change:

```bash
pixi run uv run --frozen --project genesis_tools python tests/verify_pipeline.py
pixi run checks
git diff --check
git status --short
```

While editing, `pixi run checks-fast` runs everything except the Nextflow pipeline
regressions in about 15 seconds. The test scripts run in parallel, and each runs its
independent phases and Nextflow cases in parallel, each in its own launch folder. At most
`GENESIS_CHECK_JOBS` (default 4) run at once at each level, and a failing script's full
output is printed at the end. Each test run of Nextflow may use half the CPUs and half of
the memory, of Docker's VM for Docker validation and of the machine otherwise. Set
`GENESIS_CHECK_JOBS=1` on a shared machine such as a login node:

```bash
pixi run checks-fast
GENESIS_CHECK_JOBS=1 pixi run checks
```

Use [Nextflow conventions](nextflow-style.md) when editing a workflow. Keep generated
outputs and machine-specific configuration outside Git. Commit a reviewable change
with the commands actually tested and any remaining limitations.

Prepare public validation evidence from a separate private source directory:

```bash
pixi run uv run --frozen --project genesis_tools python -m genesis_tools.publication export /private/evidence /staging/public-evidence
pixi run uv run --frozen --project genesis_tools python -m genesis_tools.publication check /staging/public-evidence
```

The exporter sanitizes review-state hashes in logs and command arguments, along with
recognized credential patterns. It checks archive paths and expansion limits and
records original and published file/member checksums in `publication-manifest.json`.
Keep original receipts and their hashes: sanitized copies do not have the same bytes.
Review the output before committing it. This targeted check is included in
`pixi run checks`; it complements secret scanning rather than certifying that arbitrary
binary data or every credential format is safe to publish. Keep private originals,
environment files and assistant configuration outside Git.

Real-tool acceptance is opt-in: `pixi run validate-docker` runs DAP-seq with its pinned
containers on synthetic inputs. For bulk ATAC, scATAC and RNA-seq, run
`tests/verify_atac.py`, `tests/verify_scatac.py` or `tests/verify_rna.py` with `--docker`,
or `--conda` to build the tools from `environments/*.yaml` where Docker is unavailable.

## Data and configuration

Git transfers code, locks, small fixtures, schemas, templates and validation reports.
It does not transfer references, reads, analysis outputs, caches, live databases,
model weights or credentials. Public container pins for both workflows are in
`conf/genesis-images.lock`. It contains image references only.

Keep `.env`, `.env.*`, `.envrc`, `.claude/` and `.codex/` local. Git and Docker build
contexts ignore them. Genesis reads image pins from `conf/genesis-images.lock`; it
does not load a local `.env`. Check staged files with `git diff --cached --name-only`
before committing.

Keep site-specific configuration outside the checkout. Test schedulers in the
environment where they will run; local checks do not establish cluster compatibility.
[Validation records](../validation/README.md) show which scientific and installation
checks have been performed.