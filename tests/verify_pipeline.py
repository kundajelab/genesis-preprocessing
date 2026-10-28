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
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import zipfile
from collections import Counter
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

from genesis_tools.metadata import chromosome_sizes, infer_metadata, inspect_read_length
from genesis_tools.samples import SAMPLE_HEADER, load_samples

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "Arabidopsis_thaliana/dap-seq-GSE60141": (936, 934, "SE"),
    "Arabidopsis_lyrata/dap-seq-PRJNA1177479": (405, 378, "PE"),
    "Arabidopsis_thaliana/dap-seq-PRJNA1177481": (800, 748, "PE"),
    "Sorghum_bicolor/dap-seq-PRJNA1177471": (142, 134, "PE"),
}
SHEET_DIGESTS = {
    "Arabidopsis_thaliana/dap-seq-GSE60141": (
        "24b5929ae27606d655cf1f190e6cb8703e1db89cf3442e8b2413e0d615ae20d5"
    ),
    "Arabidopsis_lyrata/dap-seq-PRJNA1177479": (
        "ee3323ecfbf0010dd727291da44699065c17e017b4ad155b6391f9289098ab9f"
    ),
    "Arabidopsis_thaliana/dap-seq-PRJNA1177481": (
        "4ceb6b00bb6393ce8702293a0d38dbc9c64c27b30bce5b842fa0bc9dc850cf6f"
    ),
    "Sorghum_bicolor/dap-seq-PRJNA1177471": (
        "20df482448afde0a126e9f61de367115ff7803c12449fe1d63e732b0883c233d"
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


def check_sheets(work: Path) -> None:
    """Check all four migrated datasets without requiring real reference downloads."""
    references = work / "sheet-references"
    references.mkdir()
    for dataset, (count, treatments, layout) in EXPECTED.items():
        sheet = ROOT / "samplesheets" / f"{dataset}.tsv"
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
        assert {row["species"] for row in samples} == {dataset.split("/")[0]}
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
    (references / "ref.fna.gz").touch()
    (references / "ref.fasta.gz").touch()
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
        (
            [control, {**control, "sample_id": "other", "reference_fasta": "ref.fna.gz"}],
            "reference ID collision",
        ),
        (
            [control, {**control, "sample_id": "other", "reference_fasta": "ref.fasta.gz"}],
            "reference ID collision",
        ),
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
        str(ROOT / "workflows/dap-seq.nf"),
        # The workflow lives in workflows/, so pass its config explicitly, as the launcher does.
        "-c",
        str(ROOT / "conf/dap-seq.config"),
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
    trace = work / f"trace-{suffix}.tsv"
    if not trace.exists():  # A run rejected before any task may not write a trace.
        return []
    with trace.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check_outputs(work: Path, *, docker: bool) -> None:
    """Assert treatment/control wiring, numerical scores, and safe publication."""
    rows = trace_rows(work)
    counts = Counter(row["name"].split(" ")[0] for row in rows)
    assert all(row["status"] == "COMPLETED" for row in rows), rows
    assert counts == {
        "GET_REFERENCES:CHROM_SIZES": 2,
        "GET_REFERENCES:BWA_MEM2_INDEX": 2,
        "DOWNLOAD": 3,
        "FASTQC": 7,
        "MULTIQC": 1,
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
    aggregate = published / "multiqc"
    assert (aggregate / "multiqc_report.html").stat().st_size > 0
    assert len(list(aggregate.glob("*.html"))) == 1
    provenance = json.loads((aggregate / "multiqc_data/genesis_provenance.json").read_text())
    assert provenance["run"]["name"] == "samples"
    assert provenance["run"]["pipeline_version"] == "0.1.0"
    assert len(provenance["run"]["git_sha"]) == 40
    # The run records its image pins and source checkout beside its outputs, not task receipts.
    assert (published / "genesis-images.lock").read_bytes() == (
        ROOT / "conf/genesis-images.lock"
    ).read_bytes()
    status = (published / "git-status.txt").read_text().splitlines()
    assert (
        status[0].startswith("branch: ") and status[1] == "commit: " + provenance["run"]["git_sha"]
    )
    assert status[2] == "container engine: " + ("docker" if docker else "none"), status
    assert status[3] in ("dirty files: none", "dirty files:"), status
    # Porcelain lines keep their two positional status columns, e.g. " M path" or "?? path".
    assert all(re.fullmatch(r"[ MADRCU?!]{2} .+", line) for line in status[4:]), status
    assert (
        not list(published.rglob("*.genesis-task.json")) and not (published / "provenance").exists()
    )
    identities = {sample["id"]: sample for sample in provenance["samples"]}
    assert set(identities) == set(SAMPLES)
    assert (
        identities["se_treatment"]["control"] == identities["se_second"]["control"] == "se_control"
    )
    if docker:
        html = (aggregate / "multiqc_report.html").read_text()
        assert provenance["run"]["git_sha"] in html
        for section in ("fastqc", "samtools-stats", "samtools-flagstat", "samtools-idxstats"):
            assert f'id="{section}"' in html, section
        data = json.loads((aggregate / "multiqc_data/multiqc_data.json").read_text())
        parsed = data["report_saved_raw_data"]
        assert len(parsed["multiqc_fastqc"]) == 7
        assert len(parsed["multiqc_samtools_stats"]) == 10
        assert len(parsed["multiqc_samtools_flagstat"]) == 5
        assert len(parsed["multiqc_samtools_idxstats"]) == 5
        assert set(parsed["multiqc_genesis_metadata"]) == set(SAMPLES)
        for sample, layout in SAMPLES.items():
            reads = 6000 if sample.endswith("control") else 32000
            assert parsed["multiqc_samtools_stats"][f"{sample}.qc.report.main.stats.txt"][
                "raw_total_sequences"
            ] == reads * (2 if layout == "PE" else 1)
            assert (
                parsed["multiqc_samtools_stats"][f"{sample}.qc.report.read1.stats.txt"][
                    "raw_total_sequences"
                ]
                == reads
            )
            assert parsed["multiqc_samtools_flagstat"][f"{sample}.qc.report.flagstat.txt"][
                "flagstat_total"
            ] == reads * (2 if layout == "PE" else 1)
            assert "chr1" in parsed["multiqc_samtools_idxstats"][f"{sample}.qc.report.idxstats.tsv"]
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
        mates = ("read1", "read2") if layout == "PE" else ("read1",)
        fastqc = folder / "qc" / "fastqc"
        assert len(list(fastqc.glob("*_fastqc.html"))) == len(mates)
        assert len(list(fastqc.glob("*_fastqc.zip"))) == len(mates)
        assert len(list(fastqc.iterdir())) == 3 * len(mates)
        for mate in mates:
            stem = f"{sample}.{mate}_fastqc"
            assert (fastqc / f"{stem}.html").stat().st_size > 0
            assert (fastqc / f"{stem}.zip").stat().st_size > 0
            if docker:
                assert (fastqc / f"{sample}.{mate}.fastqc.version.txt").read_text().strip() == (
                    "FastQC v0.12.1"
                )
                with zipfile.ZipFile(fastqc / f"{stem}.zip") as archive:
                    assert archive.testzip() is None
                    data = archive.read(f"{stem}/fastqc_data.txt").decode()
                    assert data.startswith("##FastQC\t0.12.1")
                    assert f"Filename\t{sample}.{mate}.fastq.gz" in data
                    expected_reads = 6000 if sample.endswith("control") else 32000
                    assert f"Total Sequences\t{expected_reads}\n" in data
                    assert "Sequence length\t75\n" in data
                    assert archive.read(f"{stem}/summary.txt")
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
    fastqc_tasks = [row for row in rows if row["name"].startswith("FASTQC ")]
    expected_tags = {
        f"FASTQC ({sample}:{mate})"
        for sample, layout in SAMPLES.items()
        for mate in (("read1", "read2") if layout == "PE" else ("read1",))
    }
    assert {row["name"] for row in fastqc_tasks} == expected_tags
    if docker:
        # Verify compressed bytes through serving, downloading, and FastQC staging.
        for row in fastqc_tasks:
            prefix = work / "work" / row["hash"]
            folders = list(prefix.parent.glob(prefix.name + "*"))
            assert len(folders) == 1
            reads = list(folders[0].glob("*.fastq.gz"))
            assert len(reads) == 1
            read = reads[0]
            sample, mate = read.name.removesuffix(".fastq.gz").rsplit(".", 1)
            original = work / "inputs" / sample
            if mate == "read2":
                original /= "mate"
            original /= "reads.fastq.gz"
            assert (
                hashlib.sha256(read.read_bytes()).digest()
                == hashlib.sha256(original.read_bytes()).digest()
            ), read
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
        if docker and "bamCoverage " in source:
            assert (
                "Matplotlib created a temporary cache directory"
                not in (command.parent / ".command.err").read_text()
            )
        assert "macs3_launch" not in source, command
        if "macs3 callpeak" in source:
            # MACS3 is invoked directly, bounded by the task's CPUs.
            assert "_control.primary.bam" in source
            assert re.search(r"--threads [0-9]+\n", source), source
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
    fastqc_before = {
        path: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (work / "published").glob("samples/output/*/*/qc/fastqc/*")
    }
    assert len(fastqc_before) == 21
    assert run_pipeline(work, sheet, docker=docker, resume=True).returncode == 0
    rows = trace_rows(work, "resume")
    reference_tasks = [
        row
        for row in rows
        if row["name"].split(" ")[0]
        in {"GET_REFERENCES:CHROM_SIZES", "GET_REFERENCES:BWA_MEM2_INDEX"}
    ]
    if docker:
        # Published references are reused from disk without any reference task.
        assert not reference_tasks, reference_tasks
    else:
        # Stub runs publish no references, so they are rebuilt, from Nextflow's task cache.
        assert len(reference_tasks) == 4, reference_tasks
        assert all(row["status"] == "CACHED" for row in reference_tasks), reference_tasks
        assert sorted(p.name for p in references.iterdir()) == ["pe.fa.gz", "se.fa.gz"]
    fastqc_tasks = [row for row in rows if row["name"].split(" ")[0] == "FASTQC"]
    assert len(fastqc_tasks) == 7 and all(row["status"] == "CACHED" for row in fastqc_tasks)
    # Newly published references change staged input paths on the first resume.
    # A second resume must cache all downstream computation at those stable paths.
    downloads = [row for row in rows if row["name"].split(" ")[0] == "DOWNLOAD"]
    assert len(downloads) == 3 and all(row["status"] == "CACHED" for row in downloads)
    multiqc_before = {
        str(path.relative_to(work / "published")): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (work / "published/samples/output/multiqc").rglob("*")
        if path.is_file()
    }
    assert run_pipeline(work, sheet, docker=docker, resume=True, suffix="stable").returncode == 0
    rows = trace_rows(work, "stable")
    for process in (
        "DOWNLOAD",
        "FASTQC",
        "MULTIQC",
        "BWA_MEM2_ALIGN",
        "QUANTIFY",
    ):
        tasks = [row for row in rows if row["name"].split(" ")[0] == process]
        assert tasks and all(row["status"] == "CACHED" for row in tasks), tasks
    assert multiqc_before == {
        str(path.relative_to(work / "published")): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (work / "published/samples/output/multiqc").rglob("*")
        if path.is_file()
    }
    assert before == {path: (path.stat().st_mtime_ns, path.read_bytes()) for path in before}
    assert fastqc_before == {
        path: hashlib.sha256(path.read_bytes()).hexdigest() for path in fastqc_before
    }


def check_failures(work: Path, sheet: Path) -> None:
    """A rejected sheet must stop before download; invalid parameters fail early."""
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    bad = work / "bad.tsv"
    write_sheet(bad, [{**rows[0], "control_sample": "missing"}])
    # The whole sheet is checked before any task starts, so a rejected sheet runs nothing.
    assert run_pipeline(work, bad, suffix="invalid-sheet", expect_failure=True).returncode != 0
    assert trace_rows(work, "invalid-sheet") == []
    assert "invalid control missing" in (work / "output-invalid-sheet.log").read_text()
    schema = work / "schema-invalid.tsv"
    write_sheet(
        schema,
        [
            rows[0],
            {**rows[1], "sample_id": rows[0]["sample_id"]},
            {**rows[2], "read1_url": "ftp://example.org/reads.fastq.gz"},
        ],
    )
    assert run_pipeline(work, schema, suffix="schema-sheet", expect_failure=True).returncode != 0
    assert trace_rows(work, "schema-sheet") == []
    log = (work / "output-schema-sheet.log").read_text()
    assert "Validation of samplesheet failed" in log and "duplicate entries" in log, log[-2000:]
    assert "read1_url must be an HTTP or HTTPS URL" in log, log[-2000:]
    lines = sheet.read_text().splitlines(keepends=True)
    for suffix, text, diagnostic in (
        # nf-schema reports a missing column; the header check rejects what it allows.
        (
            "missing-column",
            sheet.read_text().replace("control_sample", "control", 1),
            "Missing required field(s): control_sample",
        ),
        (
            "extra-column",
            "".join(line.rstrip("\n") + "\tnote\n" for line in lines),
            "expected header",
        ),
        (
            "reordered-columns",
            "".join(
                "\t".join(f[1:2] + f[:1] + f[2:]) + "\n"
                for f in (line.rstrip("\n").split("\t") for line in lines)
            ),
            "expected header",
        ),
    ):
        candidate = work / f"{suffix}.tsv"
        candidate.write_text(text)
        assert run_pipeline(work, candidate, suffix=suffix, expect_failure=True).returncode != 0
        assert trace_rows(work, suffix) == []
        log = (work / f"output-{suffix}.log").read_text()
        assert diagnostic in log, (suffix, log[-2000:])
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


def check_download_failures(work: Path, base_url: str) -> None:
    """Empty or corrupt raw libraries must fail before reaching FastQC or alignment."""
    for name, contents in {
        "empty": gzip.compress(b"", mtime=0),
        "corrupt": b"\x1f\x8b\x08\x00truncated",
    }.items():
        raw = work / "inputs" / f"{name}.fastq.gz"
        raw.write_bytes(contents)
        case = work / f"download-{name}"
        case.mkdir()
        (case / "references").symlink_to(work / "references", target_is_directory=True)
        shutil.copyfile(work / "limits.config", case / "limits.config")
        sheet = case / "samples.tsv"
        write_sheet(
            sheet,
            [
                {
                    "sample_id": "bad_control",
                    "species": "Plant_SE",
                    "read1_url": f"{base_url}/{raw.name}",
                    "read2_url": "-",
                    "control_sample": "-",
                    "reference_fasta": "se.fa.gz",
                }
            ],
        )
        assert run_pipeline(case, sheet, docker=True, expect_failure=True).returncode != 0
        rows = trace_rows(case)
        assert any(
            row["name"].startswith("DOWNLOAD ") and row["status"] == "FAILED" for row in rows
        ), rows
        assert not any(row["name"].split(" ")[0] in {"FASTQC", "BWA_MEM2_ALIGN"} for row in rows), (
            rows
        )
        assert raw.read_bytes() == contents


def check_bams(work: Path) -> None:
    """Read real BAMs to verify MAPQ-zero retention, flags, and untrimmed sequences."""
    images = dict(
        line.split("=", 1)
        for line in (ROOT / "conf/genesis-images.lock").read_text().splitlines()
        if "=" in line
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
                    images["GENESIS_ALIGNMENT_IMAGE"],
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


def check_reference_failures(work: Path, sheet: Path, *, docker: bool) -> None:
    """Reject sheets with unusable references before any task, using copied references."""
    parent = work / "reference-failures"
    parent.mkdir()
    original = {
        p.relative_to(work / "references"): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (work / "references").rglob("*")
        if p.is_file()
    }
    for scenario, diagnostic in (
        ("id-collision", "reference ID collision"),
        ("missing-reference", "reference missing: absent.fa.gz"),
    ):
        case = parent / scenario
        case.mkdir()
        references = case / "references"
        shutil.copytree(work / "references", references)
        shutil.copy2(work / "limits.config", case / "limits.config")
        samples = list(csv.DictReader(sheet.open(), delimiter="\t"))
        control = next(row for row in samples if row["sample_id"] == "se_control")
        if scenario == "id-collision":
            shutil.copy2(references / "se.fa.gz", references / "se.fna.gz")
            samples.append({**control, "sample_id": "collision", "reference_fasta": "se.fna.gz"})
        else:
            samples.append({**control, "sample_id": "orphan", "reference_fasta": "absent.fa.gz"})
        candidate = case / "samples.tsv"
        write_sheet(candidate, samples)
        assert run_pipeline(case, candidate, docker=docker, expect_failure=True).returncode != 0
        assert diagnostic in (case / "output-run.log").read_text()
        assert trace_rows(case) == []
    assert original == {
        p.relative_to(work / "references"): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (work / "references").rglob("*")
        if p.is_file()
    }


INDEX_SUFFIXES = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")


def reference_record(references: Path, ref_id: str, fasta: str) -> str:
    """The <fasta>.sha256 content BWA_MEM2_INDEX writes: FASTA hash, then product sizes."""
    lines = [f"{hashlib.sha256((references / fasta).read_bytes()).hexdigest()}  {fasta}"]
    for relative in [f"{ref_id}.chrom.sizes"] + [
        f"{ref_id}.bwa-mem2/genome{s}" for s in INDEX_SUFFIXES
    ]:
        lines.append(f"{(references / relative).stat().st_size}  {relative}")
    return "\n".join(lines) + "\n"


def synthetic_products(references: Path) -> None:
    """Stand-in reference products with a valid record, as a real run would publish."""
    for ref_id in ("se", "pe"):
        (references / f"{ref_id}.chrom.sizes").write_text("chr1\t2000000\n")
        (references / f"{ref_id}.bwa-mem2").mkdir()
        for suffix in INDEX_SUFFIXES:
            (references / f"{ref_id}.bwa-mem2/genome{suffix}").write_text(f"{ref_id}{suffix}\n")
        (references / f"{ref_id}.fa.gz.sha256").write_text(
            reference_record(references, ref_id, f"{ref_id}.fa.gz")
        )


def check_reference_rebuild(work: Path, sheet: Path, *, docker: bool) -> None:
    """Rebuild a reference only when its FASTA hash or a recorded product size is off."""
    parent = work / "reference-rebuild"
    parent.mkdir()
    built = work / "references"
    if docker:
        # The real run published products and a record the reader must accept.
        assert (built / "se.fa.gz.sha256").read_text() == reference_record(built, "se", "se.fa.gz")

    def changed_fasta(references: Path) -> None:
        fasta = references / "se.fa.gz"
        content = gzip.decompress(fasta.read_bytes())
        fasta.write_bytes(gzip.compress(content.replace(b"A", b"T", 1), mtime=0))

    def resized_index(references: Path) -> None:
        index = references / "se.bwa-mem2/genome.ann"
        index.write_bytes(index.read_bytes() + b"x")  # e.g. replaced or partly copied

    def hash_only(references: Path) -> None:
        record = references / "se.fa.gz.sha256"
        record.write_text(record.read_text().splitlines()[0] + "\n")  # earlier one-line format

    cases: list[tuple[str, Callable[[Path], object] | None, bool]] = [
        ("unchanged", None, False),
        ("changed-fasta", changed_fasta, True),
        ("missing-record", lambda r: (r / "se.fa.gz.sha256").unlink(), True),
        ("empty-index", lambda r: (r / "se.bwa-mem2/genome.ann").write_bytes(b""), True),
        ("resized-index", resized_index, True),
        ("hash-only-record", hash_only, True),
        (
            "stale-hash",
            lambda r: (r / "se.fa.gz.sha256").write_text(
                "0" * 64 + (r / "se.fa.gz.sha256").read_text()[64:]
            ),
            True,
        ),
    ]
    for scenario, change, rebuilt in cases:
        case = parent / scenario
        case.mkdir()
        references = case / "references"
        shutil.copytree(built, references)
        if not docker:
            synthetic_products(references)
        shutil.copy2(work / "limits.config", case / "limits.config")
        if change:
            change(references)
        before = {
            path: path.read_bytes() for path in sorted(references.rglob("*")) if path.is_file()
        }
        assert run_pipeline(case, sheet, docker=docker).returncode == 0, scenario
        names = [row["name"] for row in trace_rows(case)]
        for process in ("GET_REFERENCES:CHROM_SIZES", "GET_REFERENCES:BWA_MEM2_INDEX"):
            # Only the altered SE reference is rebuilt; the PE reference is always reused.
            assert (f"{process} (se)" in names) == rebuilt, (scenario, names)
            assert f"{process} (pe)" not in names, (scenario, names)
        after = {
            path: path.read_bytes() for path in sorted(references.rglob("*")) if path.is_file()
        }
        if docker:
            # A rebuild publishes fresh products and a record matching them.
            record = (references / "se.fa.gz.sha256").read_text()
            assert record == reference_record(references, "se", "se.fa.gz"), scenario
        else:
            # Stub runs never publish into the references folder.
            assert after == before, scenario


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
    passed = False
    try:
        assert (
            hashlib.sha256(
                (ROOT / "vendor/phantompeakqualtools/run_spp.R").read_bytes()
            ).hexdigest()
            == "778511418f32602383da97526a8f56033fa136577372dff14d014bef57866aaa"
        )
        if not args.docker:
            check_sheets(work)
            check_metadata(work)
        (work / "inputs").mkdir()
        with serve(work / "inputs", docker=args.docker) as base_url:
            sheet = create_fixture(work, base_url, docker=args.docker)
            raw_before = {
                str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (work / "inputs").rglob("*.fastq.gz")
            }
            (work / "raw-input-sha256.json").write_text(json.dumps(raw_before, indent=2) + "\n")
            assert run_pipeline(work, sheet, docker=args.docker).returncode == 0
            check_outputs(work, docker=args.docker)
            if args.docker:
                check_bams(work)
            check_resume(work, sheet, docker=args.docker)
            check_reference_failures(work, sheet, docker=args.docker)
            check_reference_rebuild(work, sheet, docker=args.docker)
            assert raw_before == {
                path: hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in raw_before
            }
            if args.docker:
                check_download_failures(work, base_url)
            else:
                check_failures(work, sheet)
        if args.docker:
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tests/verify_fastqc.py"),
                    "--work",
                    str(work / "fastqc-cases"),
                ],
                check=True,
                timeout=900,
            )
        if args.docker:
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "tests/verify_multiqc.py"),
                    "--source",
                    str(work / "published/samples/output"),
                    "--manifest",
                    str(
                        work
                        / "published/samples/output/multiqc/multiqc_data/genesis_provenance.json"
                    ),
                    "--work",
                    str(work / "multiqc-cases"),
                ],
                check=True,
                timeout=900,
            )
        label = "Real-tool Docker end-to-end" if args.docker else "Nextflow regression"
        print(f"{label} validation passed.", flush=True)
        passed = True
    finally:
        if args.keep or not passed:
            print(f"Artifacts retained: {work}", flush=True)
        else:
            shutil.rmtree(work)


if __name__ == "__main__":
    main()
