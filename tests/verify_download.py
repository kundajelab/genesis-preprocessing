"""Offline checks of DOWNLOAD (modules/common/download.nf), run through tests/download_harness.nf.

A local server stands in for NCBI's SRA locator, AWS and ENA, and stand-in fastq-dump, vdb-dump
and pigz go first on PATH, then workflows/bin (sra-fetch, sra-split). DOWNLOAD_SRA must convert
SRA runs with ENA's read names, and record a sample for DOWNLOAD_URLS to fetch from its FASTQ URLs
when its SRA run is unavailable, corrupt, or of another layout; a conversion that fails otherwise
stops the run.
"""

from __future__ import annotations

import csv
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

# Streams FASTQ records as fastq-dump -N FIRST -X LAST --split-spot -Z does from the stand-in .sra
# file (JSON: the run and each mate's reads as [name, sequence, quality]): each spot's reads in
# turn, headers formatted with --defline-seq and --defline-qual, and its count on stderr. Appends
# its arguments to fastq-dump.args. Fails, as on a corrupt file, when the file has "crash" set.
FAKE_FASTQ_DUMP = """#!{python}
import json, sys
from pathlib import Path
args = sys.argv[1:]
with open("fastq-dump.args", "a") as log:
    log.write(json.dumps(args) + "\\n")
option = lambda name: args[args.index(name) + 1]
run = json.loads(Path(args[-1]).read_text())
if run.get("crash"):
    sys.exit("fastq-dump stand-in: crashed")
accession, mates = run["run"], run["mates"]
first, last = int(option("-N")), int(option("-X"))
for spot in range(first, last + 1):
    for mate, reads in enumerate(mates, 1):
        read, sequence, quality = reads[spot - 1]
        fields = {{"$ac": accession, "$si": str(spot), "$sn": read, "$ri": str(mate)}}
        lines = [option("--defline-seq"), sequence, option("--defline-qual"), quality]
        for key, value in fields.items():
            lines[0] = lines[0].replace(key, value)
            lines[2] = lines[2].replace(key, value)
        sys.stdout.write("\\n".join(lines) + "\\n")
print(f"Read {{last - first + 1}} spots for {{accession}}", file=sys.stderr)
"""
# vdb-dump --info FILE: the stand-in .sra file's spot count, as SEQ. vdb-dump --table_enum FILE:
# its tables, with PRIMARY_ALIGNMENT when it has "aligned" set.
FAKE_VDB_DUMP = """#!{python}
import json, sys
from pathlib import Path
run = json.loads(Path(sys.argv[-1]).read_text())
if "--table_enum" in sys.argv:
    print("tbl : PRIMARY_ALIGNMENT" if run.get("aligned") else "tbl : SEQUENCE")
else:
    print(f"SEQ    : {{len(run['mates'][0]):,}}")
"""
# pigz -p N -4 (stdin to stdout, as sra-split calls it) and pigz -dc FILE (as DOWNLOAD_URLS
# does); each call's arguments are appended to pigz.args.
FAKE_PIGZ = """#!{python}
import gzip, json, shutil, sys
args = sys.argv[1:]
with open("pigz.args", "a") as log:
    log.write(json.dumps(args) + "\\n")
if "-dc" in args:
    with gzip.open(args[-1], "rb") as stream:
        shutil.copyfileobj(stream, sys.stdout.buffer)
else:
    sys.stdout.buffer.write(gzip.compress(sys.stdin.buffer.read(), mtime=0))
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
    root: Path,
    base: str,
    run: str,
    mates: list[list[list[str]]],
    *,
    md5: str = "",
    crash: bool = False,
    aligned: bool = False,
) -> None:
    """Serve a stand-in .sra file for run, and the locator's answer for it, with md5 if given."""
    contents = json.dumps({"run": run, "mates": mates, "crash": crash, "aligned": aligned}).encode()
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
        "-with-trace",
        str(folder / "trace.tsv"),
        "-work-dir",
        str(folder / "work"),
        "-ansi-log",
        "false",
    ]
    environment = {
        **os.environ,
        # Nextflow adds the main script's bin/ to PATH, and the harness is in tests/.
        "PATH": ":".join(
            [str(tools), str(ROOT / "workflows/bin"), str(Path(sys.executable).parent)]
            + [os.environ["PATH"]]
        ),
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


def task_folders(folder: Path) -> dict[tuple[str, str], list[Path]]:
    """Each attempt's task folder, from the harness's trace, by process and sample: the process
    is DOWNLOAD_SRA or DOWNLOAD_URLS, and the sample its task's tag."""
    with (folder / "trace.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    tasks: dict[tuple[str, str], list[Path]] = {}
    for row in rows:
        process, tag = row["name"].split(" ", 1)
        (path,) = (folder / "work").glob(f"{row['hash']}*")
        key = (process.removeprefix("DOWNLOAD:"), tag.strip("()"))
        tasks.setdefault(key, []).append(path)
    return tasks


def main() -> None:
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="download-", dir=runs))
    tools = work / "tools"
    tools.mkdir()
    stand_ins = {"fastq-dump": FAKE_FASTQ_DUMP, "vdb-dump": FAKE_VDB_DUMP, "pigz": FAKE_PIGZ}
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
            check_crash(work, tools, root, base)
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
    """SRA runs converted, and samples whose SRA run fails fetched from their URLs."""
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
    # Fetched from their URLs, each with the reason recorded: a corrupt .sra, a run the locator
    # does not know, a run of the other layout, a locator that fails (HTTP 404), and a run that
    # stores aligned reads.
    fallbacks = {
        "bad_md5": ("SRR0000003", "PE", "the .sra file's MD5 differs from NCBI's"),
        "unknown": ("SRR0000004", "SE", "NCBI's locator lists no .sra file"),
        "wrong_layout": ("SRR0000005", "PE", "its .sra file is not of the sample's layout"),
        "no_locator": ("SRR0000007", "SE", "NCBI's locator failed"),
        "aligned": ("SRR0000009", "PE", "its .sra file stores aligned reads, too slow to convert"),
    }
    for sample_id, (run, layout, _reason) in fallbacks.items():
        if sample_id == "bad_md5":
            add_sra_run(root, base, run, [reads("x"), reads("y")], md5="0" * 32)
        elif sample_id == "unknown":
            publish(
                root, base, f"{base}/sdl/{run}", b'{"version": "2","result": [{"status": 404}]}'
            )
        elif sample_id == "wrong_layout":
            add_sra_run(root, base, run, [reads("single")])
        elif sample_id == "aligned":
            add_sra_run(root, base, run, [reads("a1"), reads("a2")], aligned=True)
        urls = ena_urls(base, run, layout)
        expected[sample_id] = []
        for mate, url in enumerate(urls, 1):
            contents = gzip.compress(ena_fastq(run, mate, reads(f"ena{mate}")).encode(), mtime=0)
            publish(root, base, url, contents)
            expected[sample_id].append(contents)
        samples.append(sample(sample_id, layout, urls))
    # Not ENA's FASTQs of a run: fetched from its URL, without asking NCBI.
    plain = f"{base}/elsewhere/plain_R1.fastq.gz"
    publish(root, base, plain, gzip.compress(ena_fastq("X", 1, reads("plain")).encode(), mtime=0))
    expected["plain"] = [(root / "elsewhere/plain_R1.fastq.gz").read_bytes()]
    samples.append(sample("plain", "SE", [plain]))

    folder = work / "downloads"
    assert run_harness(folder, tools, base, samples) == 0, folder
    tasks = task_folders(folder)
    # One task per sample in DOWNLOAD_SRA, with the plain URL never asked of NCBI, and one per
    # sample fetched from its URLs: the plain one and each fallback.
    sra_ids = sorted(sample_id for sample_id in expected if sample_id != "plain")
    assert sorted(tasks) == sorted(
        [("DOWNLOAD_SRA", sample_id) for sample_id in sra_ids]
        + [("DOWNLOAD_URLS", sample_id) for sample_id in [*fallbacks, "plain"]]
    ), sorted(tasks)
    assert all(len(attempts) == 1 for attempts in tasks.values()), tasks
    task = {key: attempts[0] for key, attempts in tasks.items()}
    assert (
        "MD5 differs from NCBI's" in (task["DOWNLOAD_SRA", "bad_md5"] / ".command.err").read_text()
    )
    assert (
        "SRR0000005: reads per spot are not this PE sample's 2"
        in (task["DOWNLOAD_SRA", "wrong_layout"] / ".command.err").read_text()
    )
    # A fallback leaves nothing of its run or sample in its SRA task.
    for sample_id in fallbacks:
        leftovers = [
            path.name
            for path in task["DOWNLOAD_SRA", sample_id].iterdir()
            if path.name.endswith(".sra")
            or path.name.startswith(("sra-split.", "counts.tsv", f"{sample_id}."))
        ]
        assert not leftovers, (sample_id, leftovers)
    log = (folder / "nextflow.log").read_text()
    for sample_id, (run, _, reason) in fallbacks.items():
        assert f"{sample_id}: SRA run {run} unavailable ({reason})" in log, sample_id
    assert log.count("unavailable (") == len(fallbacks), log[-3000:]
    downloads = {
        line.split("\t")[0]: line.split("\t")[1:]
        for line in (folder / "published" / "downloads.tsv").read_text().splitlines()
    }
    assert downloads == {
        sample_id: [f"{sample_id}.read{mate}.fastq.gz" for mate in range(1, len(mates) + 1)]
        for sample_id, mates in expected.items()
    }, downloads
    for sample_id, mates in expected.items():
        source = "DOWNLOAD_SRA" if sample_id.endswith("_sra") else "DOWNLOAD_URLS"
        for mate, contents in enumerate(mates, 1):
            output = (task[source, sample_id] / f"{sample_id}.read{mate}.fastq.gz").read_bytes()
            if source == "DOWNLOAD_SRA":
                # Converted and compressed here: the decompressed bytes are ENA's.
                assert gzip.decompress(output) == contents, (sample_id, mate)
            else:
                # Fetched: the server's bytes, unchanged.
                assert output == contents, (sample_id, mate)
    sra_task = task["DOWNLOAD_SRA", "pe_sra"]
    # Converted reads are compressed fast (level 4), one pigz thread per output.
    calls = [json.loads(line) for line in (sra_task / "pigz.args").read_text().splitlines()]
    assert calls == [["-p", "1", "-4"]] * 2, calls
    # One fastq-dump per slice of the spots (the task's 1 CPU in the test profile, so one), with
    # technical reads skipped and ENA's read names.
    (arguments,) = [
        json.loads(line) for line in (sra_task / "fastq-dump.args").read_text().splitlines()
    ]
    assert arguments[:6] == ["--skip-technical", "-N", "1", "-X", "3", "--split-spot"], arguments
    assert ["--defline-seq", "@$ac.$si $sn/$ri", "--defline-qual", "+"] == arguments[-5:-1], (
        arguments
    )
    published = sorted(path.name for path in (folder / "published").rglob("*.fastq.lines"))
    assert published == sorted(f"{sample_id}.fastq.lines" for sample_id in expected), published
    for sample_id in expected:
        (counts,) = (folder / "published").rglob(f"{sample_id}.fastq.lines")
        assert counts.read_text().strip() == "12", (sample_id, counts.read_text())


def check_failure(work: Path, tools: Path, root: Path, base: str) -> None:
    """A sample whose SRA run and URLs both fail stops the run, in DOWNLOAD_URLS."""
    run = "SRR0000006"
    add_sra_run(root, base, run, [reads("z")], md5="0" * 32)
    (url,) = ena_urls(base, run, "SE")
    publish(root, base, url, b"\x1f\x8b\x08\x00truncated")
    folder = work / "failure"
    assert run_harness(folder, tools, base, [sample("broken", "SE", [url])]) != 0, folder
    tasks = task_folders(folder)
    (sra_task,) = tasks["DOWNLOAD_SRA", "broken"]
    assert "MD5 differs from NCBI's" in (sra_task / ".command.err").read_text()
    # The first attempt and its two retries.
    assert len(tasks["DOWNLOAD_URLS", "broken"]) == 3, tasks
    errors = failed_task_errors(folder / "work")
    assert "broken.read1.fastq.gz" in errors, errors[-3000:]


def check_crash(work: Path, tools: Path, root: Path, base: str) -> None:
    """A conversion that fails for another reason than the run stops the run, with no fallback."""
    run = "SRR0000008"
    add_sra_run(root, base, run, [reads("c")], crash=True)
    (url,) = ena_urls(base, run, "SE")
    publish(root, base, url, gzip.compress(ena_fastq(run, 1, reads("c")).encode(), mtime=0))
    folder = work / "crash"
    assert run_harness(folder, tools, base, [sample("crash", "SE", [url])]) != 0, folder
    tasks = task_folders(folder)
    assert sorted(tasks) == [("DOWNLOAD_SRA", "crash")], tasks
    errors = failed_task_errors(folder / "work")
    assert "fastq-dump stand-in: crashed" in errors, errors[-3000:]
    assert "sra-split failed" in errors, errors[-3000:]


if __name__ == "__main__":
    try:
        main()
    except subprocess.TimeoutExpired as error:
        raise SystemExit(f"timed out: {error}") from error
