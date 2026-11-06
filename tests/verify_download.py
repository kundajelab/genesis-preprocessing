"""Offline checks of DOWNLOAD (modules/common/download.nf), run through tests/download_harness.nf.

A local server stands in for NCBI's SRA locator, AWS and ENA, and stand-in fasterq-dump and pigz
go first on PATH. SRA runs must be converted with ENA's read names, and a sample must fall back to
its FASTQ URLs when its SRA run is missing, corrupt, or of another layout.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from verify_pipeline import (
    ROOT,
    executor_limits,
    failed_task_errors,
    run_logged,
    serve,
)

# Writes FASTQs as fasterq-dump --split-3 does from the stand-in .sra file (JSON: the run and
# each mate's reads as [name, sequence, quality]), formatting headers with --seq-defline and
# --qual-defline, and records its arguments in fasterq-dump.args.
FAKE_FASTERQ_DUMP = """#!{python}
import json, sys
from pathlib import Path
args = sys.argv[1:]
Path("fasterq-dump.args").write_text(json.dumps(args))
option = lambda name: args[args.index(name) + 1]
outdir = Path(option("--outdir"))
outdir.mkdir()
Path(option("--temp")).mkdir()
run = json.loads(Path(args[-1]).read_text())
accession, mates = run["run"], run["mates"]
for mate, reads in enumerate(mates, 1):
    name = f"{{accession}}.fastq" if len(mates) == 1 else f"{{accession}}_{{mate}}.fastq"
    with (outdir / name).open("w") as stream:
        for spot, (read, sequence, quality) in enumerate(reads, 1):
            fields = {{"$ac": accession, "$si": str(spot), "$sn": read, "$ri": str(mate)}}
            lines = [option("--seq-defline"), sequence, option("--qual-defline"), quality]
            for key, value in fields.items():
                lines[0] = lines[0].replace(key, value)
                lines[2] = lines[2].replace(key, value)
            stream.write("\\n".join(lines) + "\\n")
print(f"spots read      : {{len(mates[0]):,}}")
"""
# pigz -p N -c FILE and pigz -dc FILE, as DOWNLOAD calls it.
FAKE_PIGZ = """#!{python}
import gzip, shutil, sys
args = sys.argv[1:]
source = args[-1]
if "-dc" in args:
    with gzip.open(source, "rb") as stream:
        shutil.copyfileobj(stream, sys.stdout.buffer)
else:
    with open(source, "rb") as stream:
        sys.stdout.buffer.write(gzip.compress(stream.read(), mtime=0))
"""
# md5sum < FILE, for hosts without it (macOS).
FAKE_MD5SUM = """#!{python}
import hashlib, sys
print(hashlib.md5(sys.stdin.buffer.read()).hexdigest() + "  -")
"""

LAYOUTS = {"PE": 2, "SE": 1}


def reads(prefix: str, count: int = 3) -> list[list[str]]:
    """Reads as [name, sequence, quality]."""
    return [
        [f"{prefix}:{index}", "ACGT" * (5 + index), "I" * (20 + 4 * index)]
        for index in range(count)
    ]


def ena_fastq(run: str, mate: int, mate_reads: list[list[str]]) -> str:
    """A mate's reads as ENA writes them."""
    return "".join(
        f"@{run}.{spot} {name}/{mate}\n{sequence}\n+\n{quality}\n"
        for spot, (name, sequence, quality) in enumerate(mate_reads, 1)
    )


def ena_urls(base: str, run: str, layout: str) -> list[str]:
    folder = f"{base}/vol1/fastq/{run[:6]}/00{run[-1]}/{run}"
    if layout == "SE":
        return [f"{folder}/{run}.fastq.gz"]
    return [f"{folder}/{run}_1.fastq.gz", f"{folder}/{run}_2.fastq.gz"]


def publish(root: Path, base: str, url: str, contents: bytes) -> None:
    """Serve contents at url."""
    path = root / url.removeprefix(base + "/")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(contents)


def add_sra_run(
    root: Path, base: str, run: str, mates: list[list[list[str]]], *, md5: str = ""
) -> None:
    """Serve a stand-in .sra file for run, and the locator's answer for it, with md5 if given."""
    contents = json.dumps({"run": run, "mates": mates}).encode()
    link = f"{base}/sra/{run}/{run}"
    publish(root, base, link, contents)
    located = {
        "version": "2",
        "result": [
            {
                "bundle": run,
                "status": 200,
                "msg": "ok",
                "files": [
                    {
                        "accession": run,
                        "type": "sra",
                        "name": run,
                        "size": len(contents),
                        "md5": md5 or hashlib.md5(contents).hexdigest(),
                        "locations": [{"service": "s3", "region": "us-east-1", "link": link}],
                    }
                ],
            }
        ],
    }
    publish(root, base, f"{base}/sdl/{run}", json.dumps(located).encode())


def run_harness(folder: Path, tools: Path, base: str, samples: list[dict[str, str]]) -> int:
    folder.mkdir()
    sheet = folder / "samples.tsv"
    columns = ["sample_id", "layout", "read1_url", "read2_url"]
    sheet.write_text(
        "\t".join(columns)
        + "\n"
        + "".join("\t".join(row[c] for c in columns) + "\n" for row in samples)
    )
    limits = folder / "limits.config"
    limits.write_text(executor_limits(docker=False))
    command = [
        "nextflow",
        "-log",
        str(folder / "nextflow.log"),
        "run",
        str(ROOT / "tests/download_harness.nf"),
        # The profiles are defined with each workflow's settings; DOWNLOAD uses none of the rest.
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
        "--sra_locate_url",
        f"{base}/sdl/",
        "-work-dir",
        str(folder / "work"),
        "-ansi-log",
        "false",
    ]
    environment = {
        **os.environ,
        "PATH": f"{tools}:{Path(sys.executable).parent}:{os.environ['PATH']}",
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
    }
    print(f"Running DOWNLOAD harness in {folder}", flush=True)
    result = run_logged(
        command,
        cwd=folder,
        env=environment,
        log=folder / "output.log",
        timeout=300,
        expect_failure=True,
    )
    return result.returncode


def task_folder(work: Path) -> Path:
    (folder,) = {path.parent for path in work.glob("*/*/.command.sh")}
    return folder


def main() -> None:
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="download-", dir=runs))
    tools = work / "tools"
    tools.mkdir()
    stand_ins = {"fasterq-dump": FAKE_FASTERQ_DUMP, "pigz": FAKE_PIGZ}
    if shutil.which("md5sum") is None:
        stand_ins["md5sum"] = FAKE_MD5SUM
    for name, text in stand_ins.items():
        (tools / name).write_text(text.format(python=sys.executable))
        (tools / name).chmod(0o755)
    root = work / "server"
    root.mkdir()
    passed = False
    try:
        with serve(root, docker=False) as base:
            check_downloads(work, tools, root, base)
            check_failure(work, tools, root, base)
        passed = True
        print("DOWNLOAD SRA and URL checks passed.")
    finally:
        if passed:
            shutil.rmtree(work)


def sample(sample_id: str, layout: str, urls: list[str]) -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "layout": layout,
        "read1_url": urls[0],
        "read2_url": urls[1] if layout == "PE" else "-",
    }


def check_downloads(work: Path, tools: Path, root: Path, base: str) -> None:
    """One batch: SRA runs converted, and samples whose SRA run fails fetched from their URLs."""
    expected: dict[str, list[bytes]] = {}
    samples = []
    # Converted from SRA. Their ENA FASTQs are not served: using them would fail.
    for sample_id, run, layout in (("pe_sra", "SRR0000001", "PE"), ("se_sra", "SRR0000002", "SE")):
        mates = [reads(f"{sample_id}{mate}") for mate in range(LAYOUTS[layout])]
        add_sra_run(root, base, run, mates)
        expected[sample_id] = [
            ena_fastq(run, mate, mate_reads).encode() for mate, mate_reads in enumerate(mates, 1)
        ]
        samples.append(sample(sample_id, layout, ena_urls(base, run, layout)))
    # Fetched from their URLs: a corrupt .sra, a run the locator does not know, a run of the
    # other layout, and a URL that is not ENA's.
    fallbacks = {
        "bad_md5": ("SRR0000003", "PE"),
        "unknown": ("SRR0000004", "SE"),
        "wrong_layout": ("SRR0000005", "PE"),
    }
    for sample_id, (run, layout) in fallbacks.items():
        if sample_id == "bad_md5":
            add_sra_run(root, base, run, [reads("x"), reads("y")], md5="0" * 32)
        elif sample_id == "unknown":
            publish(
                root, base, f"{base}/sdl/{run}", b'{"version": "2","result": [{"status": 404}]}'
            )
        else:
            add_sra_run(root, base, run, [reads("single")])
        urls = ena_urls(base, run, layout)
        expected[sample_id] = []
        for mate, url in enumerate(urls, 1):
            contents = gzip.compress(ena_fastq(run, mate, reads(f"ena{mate}")).encode(), mtime=0)
            publish(root, base, url, contents)
            expected[sample_id].append(contents)
        samples.append(sample(sample_id, layout, urls))
    plain = f"{base}/elsewhere/plain_R1.fastq.gz"
    publish(root, base, plain, gzip.compress(ena_fastq("X", 1, reads("plain")).encode(), mtime=0))
    expected["plain"] = [(root / "elsewhere/plain_R1.fastq.gz").read_bytes()]
    samples.append(sample("plain", "SE", [plain]))

    folder = work / "downloads"
    assert run_harness(folder, tools, base, samples) == 0, folder
    task = task_folder(folder / "work")
    errors = (task / ".command.err").read_text()
    for sample_id, mates in expected.items():
        for mate, contents in enumerate(mates, 1):
            output = (task / f"{sample_id}.read{mate}.fastq.gz").read_bytes()
            if sample_id.endswith("_sra"):
                # Converted and compressed here: the decompressed bytes are ENA's.
                assert gzip.decompress(output) == contents, (sample_id, mate)
            else:
                # Fetched: the server's bytes, unchanged.
                assert output == contents, (sample_id, mate)
        lines = (task / f"{sample_id}.fastq.lines").read_text().strip()
        assert lines == "12", (sample_id, lines)
        fell_back = f"WARNING: {sample_id}: SRA run"
        assert (fell_back in errors) == (sample_id in fallbacks), (sample_id, errors[-3000:])
    assert "MD5 differs from NCBI's" in errors
    assert (
        "fasterq-dump wrote SRR0000005.fastq, not SRR0000005_1.fastq SRR0000005_2.fastq" in errors
    )
    assert not list(task.glob("*.sra")) and not list(task.glob("sra-*")), "conversion leftovers"
    arguments = json.loads((task / "fasterq-dump.args").read_text())
    assert arguments[:2] == ["--split-3", "--threads"], arguments
    assert ["--seq-defline", "@$ac.$si $sn/$ri", "--qual-defline", "+"] == arguments[-5:-1], (
        arguments
    )
    published = sorted(path.name for path in (folder / "published").rglob("*.fastq.lines"))
    assert published == sorted(f"{sample_id}.fastq.lines" for sample_id in expected), published


def check_failure(work: Path, tools: Path, root: Path, base: str) -> None:
    """A sample whose SRA run and URLs both fail stops the task."""
    run = "SRR0000006"
    add_sra_run(root, base, run, [reads("z")], md5="0" * 32)
    (url,) = ena_urls(base, run, "SE")
    publish(root, base, url, b"\x1f\x8b\x08\x00truncated")
    folder = work / "failure"
    assert run_harness(folder, tools, base, [sample("broken", "SE", [url])]) != 0, folder
    errors = failed_task_errors(folder / "work")
    assert "MD5 differs from NCBI's" in errors and "WARNING: broken: SRA run" in errors, errors[
        -3000:
    ]


if __name__ == "__main__":
    try:
        main()
    except subprocess.TimeoutExpired as error:
        raise SystemExit(f"timed out: {error}") from error
