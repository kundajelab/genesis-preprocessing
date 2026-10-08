"""Nextflow regressions and optional real-tool Docker validation; no test framework."""

from __future__ import annotations

import argparse
import csv
import functools
import gzip
import hashlib
import http.server
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import threading
from collections import Counter
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

from genesis_tools.metadata import chromosome_sizes, infer_metadata, inspect_read_length
from genesis_tools.samples import SAMPLE_HEADER, load_samples

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "01-Arabidopsis_thaliana-GSE60141": (936, 934, "SE"),
    "15-Arabidopsis_lyrata-PRJNA1177479": (405, 378, "PE"),
    "16-Arabidopsis_thaliana-PRJNA1177481": (800, 748, "PE"),
    "24-Sorghum_bicolor-PRJNA1177471": (142, 134, "PE"),
}
SHEET_DIGESTS = {
    "01-Arabidopsis_thaliana-GSE60141": (
        "fcb598098b6943241633f66d1a9882a1329ec9fc738658e3a7aaf3526a245a3c"
    ),
    "15-Arabidopsis_lyrata-PRJNA1177479": (
        "ee3323ecfbf0010dd727291da44699065c17e017b4ad155b6391f9289098ab9f"
    ),
    "16-Arabidopsis_thaliana-PRJNA1177481": (
        "9150a1ec44a26af02f0996d98a07e2f83e9a265a9ed271d77ccb51a794145a62"
    ),
    "24-Sorghum_bicolor-PRJNA1177471": (
        "48a02373ead600ec41af5b230f9e2f2c43415af5c86c090b5de450af4bf2bf07"
    ),
}
SAMPLES = {
    "se_treatment": "SE",
    "se_second": "SE",
    "se_control": "SE",
    "pe_treatment": "PE",
    "pe_control": "PE",
}


def rejected(function: Callable[[], object], message: str) -> None:
    """Require a validation failure with an actionable diagnostic."""
    try:
        function()
    except ValueError as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError(f"Expected rejection: {message}")


def write_sheet(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, SAMPLE_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_reads(path: Path, sequences: list[str]) -> None:
    with gzip.open(path, "wt") as stream:
        for index, sequence in enumerate(sequences):
            stream.write(f"@read{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def check_wrapper(work: Path) -> None:
    """Accept TSV-valued Nextflow options while rejecting multiple input sheets."""
    binaries = work / "bin"
    binaries.mkdir()
    executable = binaries / "nextflow"
    executable.write_text(
        f"#!{sys.executable}\nimport json, os, sys\nfrom pathlib import Path\n"
        "Path(os.environ['WRAPPER_LOG']).write_text(json.dumps(sys.argv[1:]))\n"
    )
    executable.chmod(0o755)
    sheet = work / "wrapper.tsv"
    sheet.touch()
    env = {
        **os.environ,
        "PATH": f"{binaries}:{os.environ['PATH']}",
        "WRAPPER_LOG": str(work / "wrapper.json"),
    }
    wrapper = ROOT / "scripts/run-pipeline.sh"
    default_profile = subprocess.run(
        [str(ROOT / "scripts/get-default-profile.sh")],
        capture_output=True,
        text=True,
        check=True,
        timeout=10,
    ).stdout.strip()
    workspace = work / "wrapper-workspace"
    for options, extra, profile, run_name in [
        ([], [], default_profile, "wrapper"),
        (["-p", "local,conda", "-n", "named"], [], "local,conda", "named"),
        ([], ["-profile", "local,docker", "-with-trace", "trace.tsv"], "local,docker", "wrapper"),
    ]:
        result = subprocess.run(
            [str(wrapper), "-w", str(workspace), *options, str(sheet), *extra],
            env=env,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        assert result.returncode == 0, result.stderr
        arguments = json.loads((work / "wrapper.json").read_text())
        assert arguments[arguments.index("--input") + 1] == str(sheet)
        assert arguments[arguments.index("--workspace") + 1] == str(workspace)
        assert arguments[arguments.index("--run_name") + 1] == run_name
        assert arguments.count("-profile") == 1
        assert arguments[arguments.index("-profile") + 1] == profile
        assert (workspace / run_name).is_dir()
    usage = subprocess.run(
        [str(wrapper), "--help"], capture_output=True, text=True, check=False, timeout=10
    )
    assert usage.returncode == 0 and usage.stdout.startswith("Usage:"), usage.stderr
    before = (work / "wrapper.json").read_bytes()
    for options in (
        [str(sheet)],
        ["--input", str(sheet)],
        ["--input=" + str(sheet)],
        ["--workspace", str(workspace)],
    ):
        result = subprocess.run(
            [str(wrapper), "-w", str(workspace), str(sheet), *options],
            env=env,
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        assert result.returncode == 2 and (
            "exactly one" in result.stderr or "-w/--workspace" in result.stderr
        ), result.stderr
        assert (work / "wrapper.json").read_bytes() == before


def check_sheets(work: Path) -> None:
    """Check all four migrated datasets without requiring real reference downloads."""
    references = work / "sheet-references"
    references.mkdir()
    for dataset, (count, treatments, layout) in EXPECTED.items():
        sheet = ROOT / f"{dataset}.tsv"
        with sheet.open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        for row in rows:
            (references / row["reference_fasta"]).touch()
        samples = load_samples(sheet, references)
        conserved = [
            [
                row[key]
                for key in (
                    "sample_id",
                    "read1_url",
                    "read2_url",
                    "control_sample",
                    "reference_fasta",
                )
            ]
            for row in samples
        ]
        digest = hashlib.sha256(json.dumps(conserved, separators=(",", ":")).encode()).hexdigest()
        assert digest == SHEET_DIGESTS[dataset], f"Changed sample identity or assignment: {dataset}"
        assert len(samples) == count
        assert sum(row["control_sample"] != "-" for row in samples) == treatments
        assert {row["species"] for row in samples} == {dataset.split("-")[1]}
        assert all((row["read2_url"] == "-") == (layout == "SE") for row in samples)
    # Exercise control validation independently of the real sheets.
    good = {
        "sample_id": "sample",
        "species": "Plant",
        "read1_url": "https://example.org/a.gz",
        "read2_url": "-",
        "control_sample": "control",
        "reference_fasta": "ref.fa.gz",
    }
    control = {**good, "sample_id": "control", "control_sample": "-"}
    (references / "ref.fa.gz").touch()
    sheet = work / "invalid.tsv"
    cases = [
        ([good], "invalid control"),
        ([good, good, control], "duplicate sample_id"),
        ([{**good, "control_sample": "sample"}, control], "invalid control"),
        ([good, {**control, "species": "Other"}], "invalid control"),
        (
            [good, {**control, "read2_url": "https://example.org/b.gz"}],
            "invalid control",
        ),
        ([good, {**control, "reference_fasta": "missing.fa.gz"}], "reference missing"),
        ([{**control, "read1_url": "https://example.org/a'bad"}], "invalid read1_url"),
        ([{**control, "read1_url": "file:///tmp/a.gz"}], "invalid read1_url"),
        ([{**control, "sample_id": "bad;id"}], "unsafe sample_id"),
        ([], "no samples"),
    ]
    for rows, diagnostic in cases:
        write_sheet(sheet, rows)
        rejected(lambda: load_samples(sheet, references), diagnostic)
    sheet.write_text("old\theader\n")
    rejected(lambda: load_samples(sheet, references), "expected header")


def check_metadata(work: Path) -> None:
    """Retain bounded sampling, malformed input, and atomic-output regressions."""
    sizes, reads, mate, output = [
        work / name for name in ("sizes.tsv", "reads.gz", "mate.gz", "metadata.tsv")
    ]
    sizes.write_text("chr1\t100\nchr2\t250\n")
    write_reads(reads, ["A" * 75] * 100 + ["A" * 76])
    assert inspect_read_length(reads) == (75, 100)
    write_reads(mate, ["T" * 75] * 3)
    infer_metadata(
        sample_id="sample",
        species="Plant",
        reference_fasta="ref.fa.gz",
        chrom_sizes=sizes,
        read1=reads,
        read2=mate,
        output=output,
    )
    with output.open() as stream:
        metadata = next(csv.DictReader(stream, delimiter="\t"))
    assert (
        metadata["layout"],
        metadata["genome_size"],
        metadata["analysis_read_length"],
        metadata["read2_reads_checked"],
    ) == ("PE", "350", "75", "3")
    before = output.read_bytes()
    write_reads(mate, ["T" * 76])
    rejected(
        lambda: infer_metadata(
            sample_id="sample",
            species="Plant",
            reference_fasta="ref.fa.gz",
            chrom_sizes=sizes,
            read1=reads,
            read2=mate,
            output=output,
        ),
        "mate length",
    )
    assert output.read_bytes() == before and not list(work.glob("*.partial"))
    write_reads(reads, [])
    rejected(lambda: inspect_read_length(reads), "no reads")
    write_reads(reads, ["A" * 75, "A" * 74])
    rejected(lambda: inspect_read_length(reads), "identical lengths")
    with gzip.open(reads, "wt") as stream:
        stream.write("@read\nACGT\n+\nIII\n")
    rejected(lambda: inspect_read_length(reads), "malformed FASTQ")
    for contents, diagnostic in [
        ("", "empty"),
        ("chr1\t0\n", "positive"),
        ("chr1\t1.5\n", "integer size"),
        ("chr1\t100\nchr1\t10\n", "names unique"),
    ]:
        sizes.write_text(contents)
        rejected(lambda: chromosome_sizes(sizes), diagnostic)


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        pass


@contextmanager
def serve(folder: Path, *, docker: bool) -> Iterator[str]:
    """Serve fixture reads over HTTP; yield the base URL that tasks should use."""
    # Docker Desktop forwards host.docker.internal to the host's loopback, but on Linux it maps
    # to the bridge gateway, so the server must listen beyond loopback there.
    host = "0.0.0.0" if docker and sys.platform.startswith("linux") else "127.0.0.1"
    handler = functools.partial(QuietHandler, directory=str(folder))
    server = http.server.ThreadingHTTPServer((host, 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        name = "host.docker.internal" if docker else "127.0.0.1"
        yield f"http://{name}:{server.server_address[1]}"
    finally:
        server.shutdown()
        server.server_close()


def create_fixture(work: Path, base_url: str, *, docker: bool) -> Path:
    """Make deterministic full-length SE/PE libraries with enrichment and multimappers."""
    generator = random.Random(60141)
    sequence = "".join(generator.choices("ACGT", k=2000000))
    # An exact duplicated block ensures genuine MAPQ-zero primary alignments.
    sequence = sequence[:1750000] + sequence[1700000:1700800] + sequence[1750800:]
    references = work / "references"
    references.mkdir()
    for name in ("se.fa.gz", "pe.fa.gz"):
        with gzip.open(references / name, "wt") as stream:
            stream.write(">chr1\n" + sequence + "\n")
    rows = []
    complement = str.maketrans("ACGT", "TGCA")
    for sample, layout in SAMPLES.items():
        folder = work / "inputs" / sample
        folder.mkdir(parents=True)
        control = sample.endswith("control")
        first, second = [], []
        for index in range(6000 if control else 32000):
            if index < 100:
                start = 1700000 + generator.randrange(500)
            elif not control and index < 30100:
                peak = (index - 100) // 150
                start = 20000 + peak * 8000 + generator.randrange(150)
            else:
                start = generator.randrange(500, len(sequence) - 500)
            fragment = 180 + generator.randrange(-20, 21)
            read1 = sequence[start : start + 75]
            read2 = sequence[start + fragment - 75 : start + fragment].translate(complement)[::-1]
            if index % 2:
                read1, read2 = read2, read1
            first.append(read1)
            second.append(read2)
        write_reads(folder / "reads.fastq.gz", first)
        mate_folder = folder / "mate"
        mate_folder.mkdir()
        write_reads(mate_folder / "reads.fastq.gz", second)
        prefix = layout.lower()
        rows.append(
            {
                "sample_id": sample,
                "species": f"Plant_{layout}",
                "read1_url": f"{base_url}/{sample}/reads.fastq.gz",
                "read2_url": f"{base_url}/{sample}/mate/reads.fastq.gz" if layout == "PE" else "-",
                "control_sample": "-" if control else f"{prefix}_control",
                "reference_fasta": f"{prefix}.fa.gz",
            }
        )
    sheet = work / "samples.tsv"
    write_sheet(sheet, rows)
    config = work / "limits.config"
    limits = "executor.cpus = 4\nexecutor.memory = '8 GB'\n"
    if docker:
        # Linux Docker only resolves host.docker.internal with an explicit host-gateway mapping.
        limits += (
            "docker.runOptions = '-u $(id -u):$(id -g) "
            "--add-host=host.docker.internal:host-gateway'\n"
        )
    config.write_text(limits)
    return sheet


def run_pipeline(
    work: Path,
    sheet: Path,
    *,
    docker: bool = False,
    resume: bool = False,
    suffix: str | None = None,
    extra: list[str] | None = None,
    expect_failure: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run the actual workflow; write full logs for actionable failure diagnostics."""
    suffix = suffix or ("resume" if resume else "run")
    command = [
        "nextflow",
        "-log",
        str(work / f"nextflow-{suffix}.log"),
        "run",
        str(ROOT / "main.nf"),
        "-profile",
        "local,test,docker" if docker else "local,test",
        "--input",
        str(sheet),
        "--references",
        str(work / "references"),
        "--workspace",
        str(work / "published"),
        "-work-dir",
        str(work / "work"),
        "-c",
        str(work / "limits.config"),
        "-with-trace",
        str(work / f"trace-{suffix}.tsv"),
        # Small batches exercise several DOWNLOAD tasks, including a partial final batch.
        "--download_batch_size",
        "2",
        "-ansi-log",
        "false",
    ]
    if not docker:
        command.append("-stub-run")
    if resume:
        command.append("-resume")
    command.extend(extra or [])
    print(
        f"Running {'Docker tools' if docker else 'Nextflow stubs'} ({suffix}) in {work}",
        flush=True,
    )
    environment = {
        **os.environ,
        "PATH": f"{Path(sys.executable).parent}:{os.environ['PATH']}",
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
    }
    environment.pop("FTP_PROXY", None)
    environment.pop("ftp_proxy", None)
    result = subprocess.run(
        command,
        cwd=work,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=900,
        check=False,
    )
    (work / f"output-{suffix}.log").write_text(result.stdout)
    if result.returncode and not expect_failure:
        print(result.stdout[-12000:], flush=True)
    return result


def trace_rows(work: Path, suffix: str = "run") -> list[dict[str, str]]:
    with (work / f"trace-{suffix}.tsv").open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check_outputs(work: Path, *, docker: bool) -> None:
    """Assert treatment/control wiring, numerical scores, and safe publication."""
    rows = trace_rows(work)
    counts = Counter(row["name"].split(" ")[0] for row in rows)
    assert all(row["status"] == "COMPLETED" for row in rows), rows
    assert counts == {
        "VALIDATE_SHEET": 1,
        "CHROM_SIZES": 2,
        "BWA_MEM2_INDEX": 2,
        "DOWNLOAD": 3,
        "METADATA": 5,
        "BWA_MEM2_ALIGN": 5,
        "QC": 5,
        "TRACKS": 5,
        "CALL_PEAKS": 3,
        "QUANTIFY": 3,
    }, counts
    reports = list((work / "published" / "samples" / "trace").glob("execution_report_*.html"))
    assert len(reports) == 1, reports
    latest = subprocess.run(
        [str(ROOT / "scripts/get-latest-report.sh"), str(work / "published")],
        capture_output=True,
        text=True,
        check=True,
        timeout=10,
    ).stdout.strip()
    assert latest == str(reports[0]), latest
    published = work / "published" / "samples" / "output"
    for path in published.rglob("*"):
        assert not path.name.endswith(
            (".fastq", ".fastq.gz", ".bam", ".bai", ".sam", ".bedGraph")
        ), path
        assert path.name not in (
            "quantification.json",
            "alignments.bed",
            "indexed_peaks.bed",
        ), path
    for sample, layout in SAMPLES.items():
        folder = published / f"Plant_{layout}" / sample
        assert folder.is_dir(), folder
        with (folder / f"{sample}.metadata.tsv").open() as stream:
            meta = next(csv.DictReader(stream, delimiter="\t"))
        assert (meta["layout"], meta["analysis_read_length"], meta["genome_size"]) == (
            layout,
            "75" if docker else "50",
            "2000000" if docker else "1000",
        )
        assert (folder / f"{sample}.fastq.lines").read_text().strip() == (
            ("24000" if sample.endswith("control") else "128000") if docker else "4"
        )
        assert len(list(folder.glob("*.bw"))) == (10 if docker and layout == "PE" else 5)
        assert (folder / f"{sample}.qc.report.spp.tsv").stat().st_size > 0
        if docker:
            assert (folder / f"{sample}.qc.report.spp.pdf").stat().st_size > 0
            with (folder / f"{sample}.alignment.tsv").open() as stream:
                provenance = next(csv.DictReader(stream, delimiter="\t"))
            assert provenance["read_length"] == "75"
            assert provenance["mapq_filter"] == "none" and provenance["flag_exclude"] == "2308"
            assert provenance["qc_reads"] == "full_read1_unpaired"
            assert provenance["spp_reads"] == "read1_first_50bp_unpaired"
            assert provenance["bwa_options"] == "-K 10000000 -k 19 -c 10000 -T 30"
        if sample.endswith("control"):
            assert not list(folder.glob("*.macs3*")) and not list(folder.glob("*.peaks.*.tsv"))
            continue
        with gzip.open(folder / f"{sample}.macs3_peaks.narrowPeak.gz", "rt") as stream:
            peaks = stream.read().splitlines()
        assert peaks, f"Expected enriched peaks for {sample}"
        for metric in ("RPM", "mean_RPKM"):
            scores = (folder / f"{sample}.peaks.{metric}.tsv").read_text().splitlines()
            assert [line.rsplit("\t", 1)[0] for line in scores] == peaks
            assert all(math.isfinite(float(line.rsplit("\t", 1)[1])) for line in scores)
            assert any(float(line.rsplit("\t", 1)[1]) > 0 for line in scores) if docker else True
    # Inspect executed commands, so changes to workflow wiring cannot bypass assertions.
    for command in (work / "work").glob("*/*/.command.sh"):
        source = command.read_text()
        assert "bowtie " not in source and "--sam-nh" not in source
        assert "--minMappingQuality" not in source
        if "bwa-mem2 mem" in source:
            # Main and QC/track alignments use full reads; only the SPP input is cut, to 50 bp.
            assert source.count("bwa-mem2 mem") == 3
            assert source.count("samtools view -u -F 2308 -") == 3 and " -q " not in source
            assert "trimfastq" not in source and "head " not in source
            assert source.count("substr($0, 1, bases)") == 1 and "-v bases=50 " in source
        if "macs3 callpeak" in source:
            assert "_control.primary.bam" in source
        if "Rscript" in source:
            assert "-s=-0:2:400" in source and "qc-bin/awk" in source
            assert ".spp.bam'" in source


def check_resume(work: Path, sheet: Path, *, docker: bool) -> None:
    """Ensure complete reference caches are reused and task caching is effective."""
    references = work / "references"
    before = {
        path: (path.stat().st_mtime_ns, path.read_bytes())
        for path in references.rglob("*")
        if path.is_file()
    }
    assert run_pipeline(work, sheet, docker=docker, resume=True).returncode == 0
    rows = trace_rows(work, "resume")
    names = {row["name"].split(" ")[0] for row in rows}
    assert not names.intersection({"CHROM_SIZES", "BWA_MEM2_INDEX"}), names
    # Newly published references change staged input paths on the first resume.
    # A second resume must cache all downstream computation at those stable paths.
    downloads = [row for row in rows if row["name"].split(" ")[0] == "DOWNLOAD"]
    assert len(downloads) == 3 and all(row["status"] == "CACHED" for row in downloads)
    assert run_pipeline(work, sheet, docker=docker, resume=True, suffix="stable").returncode == 0
    rows = trace_rows(work, "stable")
    for process in ("DOWNLOAD", "BWA_MEM2_ALIGN", "QUANTIFY"):
        tasks = [row for row in rows if row["name"].split(" ")[0] == process]
        assert tasks and all(row["status"] == "CACHED" for row in tasks), tasks
    assert before == {path: (path.stat().st_mtime_ns, path.read_bytes()) for path in before}


def check_failures(work: Path, sheet: Path) -> None:
    """A rejected sheet must stop before download; invalid parameters fail early."""
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    bad = work / "bad.tsv"
    write_sheet(bad, [{**rows[0], "control_sample": "missing"}])
    assert run_pipeline(work, bad, suffix="invalid-sheet", expect_failure=True).returncode != 0
    tasks = trace_rows(work, "invalid-sheet")
    assert (
        len(tasks) == 1 and tasks[0]["name"] == "VALIDATE_SHEET" and tasks[0]["status"] == "FAILED"
    ), tasks
    assert (
        run_pipeline(
            work, sheet, suffix="invalid-bwa", expect_failure=True, extra=["--bwa_seed_length", "0"]
        ).returncode
        != 0
    )
    assert "must be a positive integer" in (work / "output-invalid-bwa.log").read_text()
    assert (
        run_pipeline(
            work,
            sheet,
            suffix="multiple-sheets",
            expect_failure=True,
            extra=["--input", f"{sheet},{bad}"],
        ).returncode
        != 0
    )
    assert "Supply exactly one" in (work / "output-multiple-sheets.log").read_text()


def check_bams(work: Path) -> None:
    """Read real BAMs to verify MAPQ-zero retention, flags, and untrimmed sequences."""
    images = dict(
        line.split("=", 1) for line in (ROOT / ".env").read_text().splitlines() if "=" in line
    )
    for row in trace_rows(work):
        if not row["name"].startswith("BWA_MEM2_ALIGN "):
            continue
        folder = work / "work" / row["hash"]
        # Trace hashes abbreviate the task directory's second component.
        folders = list(folder.parent.glob(folder.name + "*"))
        assert len(folders) == 1, folders
        for bam in folders[0].glob("*.bam"):
            result = subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "--user",
                    f"{os.getuid()}:{os.getgid()}",
                    "--network",
                    "none",
                    "-v",
                    f"{work}:{work}:ro",
                    images["DAP_SEQ_ALIGNMENT_IMAGE"],
                    "samtools",
                    "view",
                    str(bam),
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=60,
            )
            records = [line.split("\t") for line in result.stdout.splitlines()]
            assert records and any(fields[4] == "0" for fields in records), bam
            assert all(not int(fields[1]) & 2308 for fields in records), bam
            # Only the SPP alignment uses read 1 cut to 50 bp; the others keep all 75 bases.
            length = 50 if bam.name.endswith(".spp.bam") else 75
            assert all(len(fields[9]) == length for fields in records), bam
            if bam.name.endswith((".qc.bam", ".spp.bam")):
                assert all(not int(fields[1]) & 1 for fields in records), bam


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--docker", action="store_true", help="Use real tools in existing Docker images"
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Retain fixtures, logs, and work directories",
    )
    args = parser.parse_args()
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="docker-" if args.docker else "regression-", dir=runs))
    try:
        assert (
            hashlib.sha256(
                (ROOT / "vendor/phantompeakqualtools/run_spp.R").read_bytes()
            ).hexdigest()
            == "778511418f32602383da97526a8f56033fa136577372dff14d014bef57866aaa"
        )
        if not args.docker:
            check_wrapper(work)
            check_sheets(work)
            check_metadata(work)
        (work / "inputs").mkdir()
        with serve(work / "inputs", docker=args.docker) as base_url:
            sheet = create_fixture(work, base_url, docker=args.docker)
            assert run_pipeline(work, sheet, docker=args.docker).returncode == 0
            check_outputs(work, docker=args.docker)
            if args.docker:
                check_bams(work)
            check_resume(work, sheet, docker=args.docker)
            if not args.docker:
                check_failures(work, sheet)
        label = "Real-tool Docker end-to-end" if args.docker else "Nextflow regression"
        print(f"{label} validation passed.", flush=True)
    finally:
        if args.keep:
            print(f"Artifacts retained: {work}", flush=True)
        else:
            shutil.rmtree(work)


if __name__ == "__main__":
    main()
