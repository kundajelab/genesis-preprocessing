"""Offline checks of CALL_PEAKS (modules/dap-seq/call_peaks.nf), run through
tests/call_peaks_harness.nf with a stand-in macs3 first on PATH.

A treatment MACS3 cannot build its fragment-size model for must still succeed: no peaks for
QUANTIFY, no MACS3 files published, and SAMPLE.failure-reason.txt published with MACS3's message.
Any other MACS3 failure must stop the run.
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from verify_pipeline import ROOT, executor_limits, failed_task_errors, run_logged

REASON = (
    "MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can only "
    "find 88!"
)
# macs3 callpeak -t TREATMENT ... -n NAME: what it does depends on the treatment file's text.
# "peaks" writes NAME_peaks.narrowPeak and NAME_peaks.xls. "no_model" prints the warnings real
# MACS3 does (its timestamps elided) after writing a partial NAME_peaks.xls, and exits 1, as
# MACS3 does. "crash" fails otherwise.
FAKE_MACS3 = """#!{python}
import sys
from pathlib import Path
args = sys.argv[1:]
option = lambda name: args[args.index(name) + 1]
mode, name = Path(option("-t")).read_text().strip(), option("-n")
if mode == "peaks":
    Path(f"{{name}}_peaks.narrowPeak").write_text("chr1\\t10\\t20\\tp1\\t5\\t.\\t1\\t2\\t3\\t4\\n")
    Path(f"{{name}}_peaks.xls").write_text("# MACS3\\n")
    sys.exit(0)
if mode == "no_model":
    Path(f"{{name}}_peaks.xls").write_text("partial\\n")
    log = "WARNING @ 24 Nov 2026 15:39:06: [47 MB] #2 "
    info = "INFO  @ 24 Nov 2026 15:39:06: [47 MB] #2 "
    print(info + "Total number of paired peaks: 88 ", file=sys.stderr)
    print(log + "{reason} Please make your MFOLD range broader and try again. If MACS3 still "
          "can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed "
          "number instead. ", file=sys.stderr)
    print(log + "Process for pairing-model is terminated! ", file=sys.stderr)
    sys.exit(1)
print("macs3 stand-in: crashed", file=sys.stderr)
sys.exit(1)
"""


def run_harness(folder: Path, tools: Path, samples: dict[str, str]) -> int:
    """Run CALL_PEAKS on SE samples, each a treatment whose text selects the stand-in's mode."""
    folder.mkdir()
    control = folder / "control.bam"
    control.write_text("control\n")
    rows = ["sample_id\tlayout\ttreatment\tcontrol"]
    for sample_id, mode in samples.items():
        treatment = folder / f"{sample_id}.bam"
        treatment.write_text(mode + "\n")
        rows.append(f"{sample_id}\tSE\t{treatment}\t{control}")
    sheet = folder / "samples.tsv"
    sheet.write_text("\n".join(rows) + "\n")
    limits = folder / "limits.config"
    limits.write_text(executor_limits(docker=False))
    command = [
        "nextflow",
        "-log",
        str(folder / "nextflow.log"),
        "run",
        str(ROOT / "tests/call_peaks_harness.nf"),
        "-c",
        str(ROOT / "conf/dap-seq.config"),
        "-profile",
        "local,test",
        "-c",
        str(limits),
        "--input",
        str(sheet),
        "--run_dir",
        str(folder / "published"),
        "-work-dir",
        str(folder / "work"),
        "-ansi-log",
        "false",
    ]
    environment = {
        **os.environ,
        "PATH": ":".join([str(tools), str(Path(sys.executable).parent), os.environ["PATH"]]),
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
    }
    print(f"Running CALL_PEAKS harness in {folder}", flush=True)
    result = run_logged(
        command,
        cwd=folder,
        env=environment,
        log=folder / "output.log",
        timeout=300,
        expect_failure=True,
    )
    return result.returncode


def check_model_failure(work: Path, tools: Path) -> None:
    """A treatment without a model gets only its failure reason; the other gets its peaks."""
    folder = work / "model"
    assert run_harness(folder, tools, {"peaks": "peaks", "weak": "no_model"}) == 0, folder
    assert (folder / "published/peaks.txt").read_text().split() == ["peaks"]
    output = folder / "published/output/Plant"
    assert sorted(path.name for path in (output / "peaks").iterdir()) == [
        "peaks.macs3_peaks.narrowPeak.gz",
        "peaks.macs3_peaks.xls",
    ]
    with gzip.open(output / "peaks/peaks.macs3_peaks.narrowPeak.gz", "rt") as stream:
        assert stream.read().startswith("chr1\t10\t20\tp1")
    assert sorted(path.name for path in (output / "weak").iterdir()) == ["weak.failure-reason.txt"]
    assert (output / "weak/weak.failure-reason.txt").read_text() == REASON + "\n"
    (weak_task,) = [path.parent for path in (folder / "work").glob("*/*/weak.failure-reason.txt")]
    assert not list(weak_task.glob("weak.macs3*")), list(weak_task.iterdir())
    errors = (weak_task / ".command.err").read_text()
    assert f"WARNING: weak: no peaks called: {REASON}" in errors, errors


def check_crash(work: Path, tools: Path) -> None:
    """Any other MACS3 failure stops the run, with no failure reason published."""
    folder = work / "crash"
    assert run_harness(folder, tools, {"broken": "crash"}) != 0, folder
    assert "macs3 stand-in: crashed" in failed_task_errors(folder / "work")
    assert not list(folder.glob("published/**/*.failure-reason.txt"))


def main() -> None:
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="call-peaks-", dir=runs))
    tools = work / "tools"
    tools.mkdir()
    (tools / "macs3").write_text(FAKE_MACS3.format(python=sys.executable, reason=REASON))
    (tools / "macs3").chmod(0o755)
    passed = False
    try:
        check_model_failure(work, tools)
        check_crash(work, tools)
        passed = True
        print("CALL_PEAKS checks passed.")
    finally:
        if passed:
            shutil.rmtree(work)


if __name__ == "__main__":
    try:
        main()
    except subprocess.TimeoutExpired as error:
        raise SystemExit(f"timed out: {error}") from error
