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
pixi run install-all --install-cli
export PATH="$HOME/.local/bin:$PATH"
genesis --help
genesis doctor
```

`install-all` uses the committed `pixi.lock`, `genesis_tools/uv.lock` and Python
selection in `genesis_tools/.python-version` (currently Python 3.14.5 with the GIL).
It creates environments for the checkout. Initial installation needs network
access or populated package caches. Recreate `.pixi/` and `.venv/` when moving a
checkout: environment paths and binaries are machine-specific. The supported Pixi
platforms are Linux x86_64 and macOS x86_64/ARM64; private vLLM deployment requires
an authorized Linux node and its separate serving environment.

`doctor` reports the selected checkout, tools on the host and in Pixi, and Docker access.
Docker, a GPU and API credentials are not prerequisites for the core development checks or mock tutorial;
scientific-container and live-model acceptance have their own prerequisites.

### Command installation and maintenance

`--install-cli` links `~/.local/bin/genesis` to this checkout's
`genesis_tools/.venv/bin/genesis`. Once its directory is on PATH, commands work from
any directory. Add the PATH export above to your shell startup file if necessary;
the installer does not edit it. A custom location uses
`pixi run install-all --install-cli --bin-dir /your/bin`. Check `command -v genesis`
if another installation is selected earlier on PATH. For a temporary session, you
can instead activate `genesis_tools/.venv/bin/activate`; plain `pixi run install-all`
installs environments without creating a user-level command.

After updating this checkout with Git, rerun `pixi run install-all --install-cli`
to synchronize locked dependencies. Reinstalling the same link is safe. Keep the
checkout in place: the link refers to its environment. If moving the checkout,
recreate its environments using the instructions above.

To uninstall just the exposed command, inspect `ls -l "$HOME/.local/bin/genesis"`
and remove that symlink with `rm "$HOME/.local/bin/genesis"` (adjust for a custom
bin directory). This retains the checkout, environments and all analysis data.
To switch checkouts, remove the previous link explicitly and run the installation
command in the new checkout. The installer refuses to replace a different link,
file or directory.

An installed wheel supports core commands independently. To use a wheel's
scientific launch commands, set `GENESIS_REPO=/absolute/path/to/genesis-preprocessing`
to a checkout with installed environments. This overrides editable-source discovery;
an invalid override fails explicitly. Scientific launchers execute the selected
checkout's code with its own Python and locked Pixi environment. They preserve
the caller's working directory until the existing pipeline wrapper selects its
run directory. No workflow is downloaded or chosen from the current directory.
The pipeline commands also need that checkout to show their workflow-specific
`--help`. General `genesis --help` and core curation commands work without it.

If `genesis doctor` reports an unavailable managed runtime, put Pixi on PATH and
run `pixi run install-all` in that checkout. Core curation works without
Pixi. Developer checks and image builds use Pixi tasks.

Read [AGENTS.md](../AGENTS.md) before editing. Create a branch for the next change:

```bash
git switch -c work/next-change
```

## Establish your baseline

```bash
pixi run checks
pixi run uv run --frozen --project genesis_tools python \
  examples/curation/offline_workflow.py /tmp/genesis-development-demo
git status --short
```

The demo requires a new, empty destination. Expected result: `status: PASS`, two
included datasets and six demonstration review decisions. It uses synthetic files
and a mock provider; it does not run sequencing pipelines or call a live model.

`checks` includes lint, formatting checks, types, Nextflow stubs, local worker
tests, prepared-study execution/collection, registry/review tests, fake HTTP,
simulated model services and README
command checks. Its initial Nextflow run needs access to the pinned plugin
registry or a verified cache containing `nf-schema@2.4.2` and `nf-dotenv@1.0.0`.
If plugin resolution fails, restore registry access or configure an approved cache
with `NXF_PLUGINS_DIR`; only enable `NXF_OFFLINE=true` after populating it. The
historical `/tmp/genesis-*` paths in validation reports identify old test runs and
are not paths to reuse on a new server. The [hybrid validation report](../validation/reports/hybrid-llm-validation.md)
records the verified plugin archive hashes.

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

Keep scientific settings explicit. A change to alignment, filtering, peak calling
or count definitions needs a comparison of the affected outputs, not just a
successful exit code. Use the [benchmark guide](../benchmarks/README.md) to design
that comparison. Curation and AI changes should leave ordinary preprocessing usable.

| Area | Implementation | Focused check |
| --- | --- | --- |
| DAP | `workflows/dap-seq.nf`, `modules/`, `genesis_tools/samples.py` | `tests/verify_pipeline.py` |
| Bulk ATAC | `workflows/atac.nf`, `modules/atac/`, `genesis_tools/atac/` | `tests/verify_supported_atac.py` |
| Worker execution | `genesis_tools/execution/` | `tests/verify_execution.py` |
| Study registration and collection | `genesis_tools/study/` | `tests/verify_study.py` |
| Registry and review | `genesis_tools/registry/`, `genesis_tools/curation/` | `tests/verify_curation.py` |
| Optional AI | `genesis_tools/llm/` | `tests/verify_ai.py` |
| Documentation | README and task guides | `tests/verify_readme.py` |

Python package paths in this table are relative to `genesis_tools/src/`.
Run a focused check with the project interpreter, then run the full suite before
publishing a code change:

```bash
pixi run uv run --frozen --project genesis_tools python tests/verify_study.py
pixi run checks
git diff --check
git status --short
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

Real-tool acceptance is opt-in. The [study guide](usage/studies.md#real-tool-acceptance-on-synthetic-inputs)
describes Docker-based scientific acceptance; the [provider](usage/ai-providers.md#opt-in-live-acceptance)
and [model-service](usage/llm-deployment.md#opt-in-real-local-or-ssh-acceptance) guides describe
live API/GPU checks. Mock tests do not establish acceptance on a remote site.

## Data and configuration

Git transfers code, locks, small fixtures, schemas, templates and validation reports.
It does not transfer references, reads, analysis outputs, caches, live databases,
model weights or credentials. Public container pins are in `conf/dap-images.lock`
for DAP and `conf/atac-images.json` for ATAC. These files contain image references only.

Keep `.env`, `.env.*`, `.envrc`, `.claude/` and `.codex/` local. Git and Docker build
contexts ignore them. Genesis reads DAP image pins from `conf/dap-images.lock`;
it does not load a local `.env`. Configure credentials through the environment-variable
or protected-file references described in the [provider guide](usage/ai-providers.md).
Check staged files with `git diff --cached --name-only` before committing.

For a catalog, use [registry backup/restore](usage/curation.md#read-only-status-backup-and-rollback).
With orchestration stopped, preserve study state and worker attempt/output directories
through your normal storage workflow. Keep worker-qualified paths accurate: a new
checkout does not relocate data or reconnect jobs under a different identity.

Keep provider credentials and site-specific configurations outside the checkout.
Test SSH, schedulers and GPUs in the environment where they will run; local checks
do not establish remote compatibility. [Validation records](../validation/README.md)
show which scientific and installation checks have been performed.

## Model-handoff inspection checks

The default suite includes `tests/verify_scatac_inspection.py`: pinned adapters,
bounded scans, invalid region windows/folds, stale evidence, rendering and input
immutability. To include binary readability with an existing track environment:

```bash
pixi run uv run --frozen --project genesis_tools python tests/verify_scatac_inspection.py --bigwig-python /path/to/existing/track-env/bin/python
```

No reader is installed automatically. The [offline inspection example](usage/scatac-model-handoff.md#offline-inspection-example)
is synthetic and cannot establish loader or training acceptance. Keep meeting
notes, local audit receipts and candidate handoff decisions under the ignored
`workspace/meeting-prep/` directory. Public documentation describes scientific
interfaces and evidence, not internal meeting logistics.
