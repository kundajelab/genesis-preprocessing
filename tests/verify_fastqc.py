"""Real-tool FastQC edge cases using the production module and deterministic raw reads."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def check_fastqc_inputs(work: Path) -> None:
    """Quality FAIL labels are advisory; malformed input must not silently lose reports."""
    work.mkdir(parents=True)
    harness = work / "main.nf"
    harness.write_text(
        f"include {{ FASTQC }} from '{ROOT / 'modules/common/fastqc'}'\n"
        "workflow {\n"
        "    FASTQC(Channel.of(tuple([id: 'case', species: 'Synthetic'], 'read1',\n"
        "        file(params.read, checkIfExists: true))))\n"
        "}\n"
    )
    image_line = next(
        line
        for line in (ROOT / "conf/dap-images.lock").read_text().splitlines()
        if line.startswith("DAP_SEQ_FASTQC_IMAGE=")
    )
    (work / "images.lock").write_text(image_line + "\n")
    config = work / "nextflow.config"
    config.write_text(
        "plugins { id 'nf-dotenv@1.0.0' }\n"
        "dotenv { filename = 'images.lock'; relative = '.' }\n"
        "docker.enabled = true\n"
        "docker.runOptions = '-u $(id -u):$(id -g)'\n"
        "params.read = null\nparams.run_dir = null\n"
    )
    valid = ("@read\n" + "A" * 75 + "\n+\n" + "I" * 75 + "\n") * 100
    cases = {
        "quality-fail": gzip.compress(valid.encode(), mtime=0),
        "empty": gzip.compress(b"", mtime=0),
        "short-quality": gzip.compress(b"@read\nACGT\n+\nIII\n", mtime=0),
        "malformed": gzip.compress(b"@read\nACGT\nnot-plus\nIIII\n", mtime=0),
        "corrupt-gzip": b"\x1f\x8b\x08\x00truncated",
    }
    observations = {}
    for name, contents in cases.items():
        folder = work / name
        folder.mkdir()
        read = folder / "case.read1.fastq.gz"
        read.write_bytes(contents)
        digest = hashlib.sha256(contents).hexdigest()
        command = [
            "nextflow",
            "-log",
            str(folder / "nextflow.log"),
            "run",
            str(harness),
            "-c",
            str(config),
            "--read",
            str(read),
            "--run_dir",
            str(folder / "published"),
            "-work-dir",
            str(folder / "work"),
            "-with-trace",
            str(folder / "trace.tsv"),
            "-ansi-log",
            "false",
        ]
        result = subprocess.run(
            command,
            cwd=work,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
            env={**os.environ, "NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true"},
        )
        (folder / "stdout.log").write_text(result.stdout)
        (folder / "stderr.log").write_text(result.stderr)
        observations[name] = {"command": command, "exit": result.returncode, "sha256": digest}
        (work / "observations.json").write_text(json.dumps(observations, indent=2) + "\n")
        assert hashlib.sha256(read.read_bytes()).hexdigest() == digest
        published = folder / "published/output/Synthetic/case/qc/fastqc"
        if name in ("quality-fail", "empty", "short-quality"):
            assert result.returncode == 0, result.stdout + result.stderr
            assert (published / "case.read1_fastqc.html").stat().st_size > 0
            with zipfile.ZipFile(published / "case.read1_fastqc.zip") as archive:
                assert archive.testzip() is None
                data = archive.read("case.read1_fastqc/fastqc_data.txt").decode()
                summary = archive.read("case.read1_fastqc/summary.txt").decode()
                if name == "quality-fail":
                    assert "FAIL\t" in summary
                    assert "Total Sequences\t100\n" in data
                elif name == "short-quality":
                    # Known FastQC 0.12.1 parser limitation, not a validity guarantee.
                    # The pipeline's existing metadata validator rejects this fixture.
                    assert "Total Sequences\t1\n" in data
                else:
                    # The upstream DOWNLOAD process rejects empty libraries; FastQC itself
                    # can report zero reads, which is not a plant-specific biological cutoff.
                    assert "Total Sequences\t0\n" in data
        else:
            assert result.returncode != 0, result.stdout
            assert not list(published.glob("*_fastqc.*"))
            task_logs = "\n".join(p.read_text() for p in (folder / "work").glob("*/*/.command.err"))
            assert "Failed to process" in task_logs, task_logs
    print("FastQC quality flags, empty/malformed/corrupt input and raw-byte checks passed.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--work", type=Path, help="Use and retain this case directory")
    args = parser.parse_args()
    if args.work is not None:
        check_fastqc_inputs(args.work.resolve())
        return
    runs = ROOT / "tests/.runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="fastqc-", dir=runs))
    try:
        check_fastqc_inputs(work / "cases")
    finally:
        if args.keep:
            print(f"Artifacts retained: {work}", flush=True)
        else:
            shutil.rmtree(work)


if __name__ == "__main__":
    main()
