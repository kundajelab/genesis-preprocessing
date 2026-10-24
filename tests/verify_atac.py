"""Bulk ATAC Nextflow regressions and optional real-tool validation; no test framework."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pysam
from verify_pipeline import serve, trace_rows, write_reads

ROOT = Path(__file__).resolve().parents[1]
HEADER = [
    "sample_id",
    "species",
    "read1_url",
    "read2_url",
    "reference_fasta",
    "reference_genes",
    "organelle_contigs",
]
ADAPTER = "CTGTCTCTTATACACATCT"
READ_LENGTH = 50
# Planted accessible sites; genes start at each, alternating strands.
PEAKS = [50000 + 40000 * index for index in range(10)]
# Reads inside the first copy of this block also match its second copy, so they get MAPQ 0.
REPEAT = (300000, 301000, 600000)


def write_sheet(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def create_fixture(work: Path, base_url: str, *, engine: str | None) -> Path:
    """Paired libraries with planted peaks, organelle reads, duplicates, repeats and adapters."""
    generator = random.Random(4105)
    chromosome = "".join(generator.choices("ACGT", k=1000000))
    start, end, copy = REPEAT
    chromosome = chromosome[:copy] + chromosome[start:end] + chromosome[copy + end - start :]
    contigs = {
        "chr1": chromosome,
        "ChrC": "".join(generator.choices("ACGT", k=20000)),
        "ChrM": "".join(generator.choices("ACGT", k=20000)),
    }
    references = work / "references"
    references.mkdir()
    for name in ("atac.fa.gz", "bare.fa.gz"):
        with gzip.open(references / name, "wt") as stream:
            stream.writelines(f">{contig}\n{sequence}\n" for contig, sequence in contigs.items())
    (references / "genes.bed").write_text(
        "".join(
            f"chr1\t{peak}\t{peak + 2000}\tgene{index}\t0\t+\n"
            if index % 2 == 0
            else f"chr1\t{peak - 2000}\t{peak}\tgene{index}\t0\t-\n"
            for index, peak in enumerate(PEAKS)
        )
    )
    filler = "".join(generator.choices("ACGT", k=READ_LENGTH))
    rows = []
    for sample, reference in (("rep1", "atac"), ("rep2", "atac"), ("bare", "bare")):
        first, second = [], []
        fragments: list[tuple[str, int, int]] = []
        for _ in range(12000):
            roll = generator.random()
            if fragments and roll < 0.05:
                fragments.append(fragments[generator.randrange(len(fragments))])
                continue
            length = generator.randrange(30, 45) if roll < 0.1 else generator.randrange(80, 300)
            if roll < 0.55:
                contig, position = "chr1", generator.choice(PEAKS) + generator.randrange(-150, 150)
            elif roll < 0.65:
                contig, position = generator.choice(["ChrC", "ChrM"]), generator.randrange(19000)
            elif roll < 0.7:
                contig, position = "chr1", start + generator.randrange(end - start - 400)
            else:
                contig, position = "chr1", generator.randrange(1000, 990000)
            fragments.append((contig, position - length // 2 if roll < 0.55 else position, length))
        for index, (contig, position, length) in enumerate(fragments):
            insert = contigs[contig][position : position + length]
            # Short inserts read through into the adapter, as Tn5 libraries do.
            read1 = (insert + ADAPTER + filler)[:READ_LENGTH]
            read2 = (reverse_complement(insert) + ADAPTER + filler)[:READ_LENGTH]
            if index % 2:
                read1, read2 = read2, read1
            first.append(read1)
            second.append(read2)
        folder = work / "inputs" / sample
        (folder / "mate").mkdir(parents=True)
        write_reads(folder / "reads.fastq.gz", first)
        write_reads(folder / "mate" / "reads.fastq.gz", second)
        rows.append(
            {
                "sample_id": sample,
                "species": "Plant_synthetic",
                "read1_url": f"{base_url}/{sample}/reads.fastq.gz",
                "read2_url": f"{base_url}/{sample}/mate/reads.fastq.gz",
                "reference_fasta": f"{reference}.fa.gz",
                "reference_genes": "genes.bed" if reference == "atac" else "-",
                "organelle_contigs": "ChrC,ChrM" if reference == "atac" else "-",
            }
        )
    sheet = work / "samples.tsv"
    write_sheet(sheet, rows)
    limits = "executor.cpus = 4\nexecutor.memory = '8 GB'\n"
    if engine == "docker":
        # Linux Docker only resolves host.docker.internal with an explicit host-gateway mapping.
        limits += (
            "docker.runOptions = '-u $(id -u):$(id -g) "
            "--add-host=host.docker.internal:host-gateway'\n"
        )
    if engine == "conda":
        # Keep the environments built for earlier runs.
        limits += f"conda.cacheDir = '{work.parent / 'conda_cache'}'\n"
    (work / "limits.config").write_text(limits)
    return sheet


def run_atac(
    work: Path,
    sheet: Path,
    *,
    engine: str | None = None,
    resume: bool = False,
    suffix: str | None = None,
    expect_failure: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run the actual workflow; write full logs for actionable failure diagnostics."""
    suffix = suffix or ("resume" if resume else "run")
    command = [
        "nextflow",
        "-log",
        str(work / f"nextflow-{suffix}.log"),
        "run",
        str(ROOT / "workflows/atac-seq.nf"),
        "-c",
        str(ROOT / "conf/atac-seq.config"),
        "-profile",
        f"local,test,{engine}" if engine else "local,test",
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
        "--download_batch_size",
        "2",
        "-ansi-log",
        "false",
    ]
    if not engine:
        command.append("-stub-run")
    if resume:
        command.append("-resume")
    print(f"Running {engine or 'stub'} bulk ATAC ({suffix}) in {work}", flush=True)
    environment = {
        **os.environ,
        "PATH": f"{Path(sys.executable).parent}:{os.environ['PATH']}",
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
    }
    result = subprocess.run(
        command,
        cwd=work,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=3600,
        check=False,
    )
    (work / f"output-{suffix}.log").write_text(result.stdout)
    if result.returncode and not expect_failure:
        print(result.stdout[-12000:], flush=True)
    return result


def task_folder(work: Path, row: dict[str, str]) -> Path:
    folder = work / "work" / row["hash"]
    # Trace hashes abbreviate the task directory's second component.
    folders = list(folder.parent.glob(folder.name + "*"))
    assert len(folders) == 1, folders
    return folders[0]


def check_outputs(work: Path, *, engine: str | None) -> None:
    """Assert the process wiring, publication and recorded provenance."""
    rows = trace_rows(work)
    assert all(row["status"] == "COMPLETED" for row in rows), rows
    counts = Counter(row["name"].split(" ")[0] for row in rows)
    assert counts == {
        "GET_REFERENCES:CHROM_SIZES": 2,
        "GET_REFERENCES:BWA_MEM2_INDEX": 2,
        "DOWNLOAD": 2,
        "FASTQC": 6,
        "ATAC_TRIM": 3,
        "ATAC_ALIGN": 3,
        "ATAC_FILTER": 3,
        "ATAC_PEAKS": 3,
        "ATAC_TRACKS": 3,
        "ATAC_QC": 3,
        "ATAC_MULTIQC": 1,
    }, counts
    published = work / "published" / "samples" / "output"
    for path in published.rglob("*"):
        assert not path.name.endswith((".fastq", ".fastq.gz", ".bam", ".bai")), path
    for sample in ("rep1", "rep2", "bare"):
        library = published / "Plant_synthetic" / sample
        for name in (
            f"{sample}.fragments.bed.gz",
            f"{sample}.fragments.bed.gz.tbi",
            f"{sample}.macs3_peaks.narrowPeak.gz",
            f"{sample}.cuts.counts.bw",
            f"qc/{sample}.cutadapt.log",
            f"qc/{sample}.marked.stats.txt",
            f"qc/{sample}.usable.flagstat.txt",
            f"qc/{sample}.complexity.tsv",
            f"qc/{sample}.atac_qc_mqc.tsv",
            f"qc/fastqc/{sample}.read1_fastqc.zip",
            f"qc/fastqc/{sample}.read2_fastqc.zip",
        ):
            assert (library / name).is_file(), library / name
    assert (published / "multiqc/multiqc_report.html").stat().st_size > 0
    assert (published / "dap-images.lock").read_bytes() == (
        ROOT / "conf/dap-images.lock"
    ).read_bytes()
    status = (published / "git-status.txt").read_text().splitlines()
    assert status[0].startswith("branch: ") and status[1].startswith("commit: "), status
    assert status[2] == f"container engine: {engine if engine == 'docker' else 'none'}", status
    log = (work / "output-run.log").read_text()
    assert "PIPELINE SUCCEEDED" in log, log[-2000:]
    if engine:
        check_results(work, published)


def check_results(work: Path, published: Path) -> None:
    """Known answers from the synthetic libraries, through the real tools."""
    for sample in ("rep1", "rep2", "bare"):
        library = published / "Plant_synthetic" / sample
        lines = [
            line
            for line in (library / f"qc/{sample}.atac_qc_mqc.tsv").read_text().splitlines()
            if not line.startswith("#")
        ]
        summary = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
        # Short inserts were trimmed back to the insert, not dropped.
        log = (library / f"qc/{sample}.cutadapt.log").read_text()
        assert "Read 1 with adapter:" in log and "Pairs written (passing filters):" in log, log
        with gzip.open(library / f"{sample}.fragments.bed.gz", "rt") as stream:
            fragments = [line.split("\t") for line in stream]
        contigs = Counter(fields[0] for fields in fragments)
        if sample == "bare":
            # No organelles declared: their fragments stay.
            assert contigs["ChrC"] and contigs["ChrM"], contigs
            assert summary["TSS_enrichment"] == "NA", summary
        else:
            assert set(contigs) == {"chr1"}, contigs
            assert float(summary["TSS_enrichment"]) > 5, summary
        # Each usable pair gives one fragment, and planted duplicates lower NRF below one.
        assert int(summary["usable_reads"]) == 2 * len(fragments), (summary, len(fragments))
        assert 0.5 < float(summary["NRF"]) < 1, summary
        assert float(summary["FRiP"]) > 0.4, summary
        with gzip.open(library / f"{sample}.macs3_peaks.narrowPeak.gz", "rt") as stream:
            peaks = [line.split("\t") for line in stream]
        found = {
            planted
            for planted in PEAKS
            for fields in peaks
            if fields[0] == "chr1" and int(fields[1]) <= planted < int(fields[2])
        }
        assert len(found) >= 9, (sample, sorted(found))
    for row in trace_rows(work):
        if not row["name"].startswith("ATAC_FILTER "):
            continue
        folder = task_folder(work, row)
        (bam,) = folder.glob("*.usable.bam")
        names: Counter[str] = Counter()
        with pysam.AlignmentFile(str(bam)) as stream:
            for read in stream:
                assert read.is_proper_pair and not read.is_duplicate, read
                assert not (read.is_secondary or read.is_supplementary or read.is_qcfail), read
                assert read.mapping_quality >= 30, read
                names[read.query_name] += 1
        # Both mates of every pair, or neither: fixmate removed mates whose partner failed.
        assert names and set(names.values()) == {2}, Counter(names.values())


def check_resume(work: Path, sheet: Path, *, engine: str | None) -> None:
    """Resumed runs reuse the downloads, and then every task."""
    assert run_atac(work, sheet, engine=engine, resume=True).returncode == 0
    rows = trace_rows(work, "resume")
    for process in ("DOWNLOAD", "FASTQC", "ATAC_TRIM"):
        tasks = [row for row in rows if row["name"].split(" ")[0] == process]
        assert tasks and all(row["status"] == "CACHED" for row in tasks), tasks
    # Newly published references change staged input paths on the first resume.
    # A second resume must cache all downstream computation at those stable paths.
    assert run_atac(work, sheet, engine=engine, resume=True, suffix="stable").returncode == 0
    rows = trace_rows(work, "stable")
    assert rows and all(row["status"] == "CACHED" for row in rows), rows


def check_failures(work: Path, sheet: Path) -> None:
    """A rejected sheet stops before any task."""
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    lines = sheet.read_text().splitlines(keepends=True)
    cases: list[tuple[str, str, str]] = []
    for suffix, changed, diagnostic in (
        ("single-end", [{**rows[0], "read2_url": "-"}], "bulk ATAC libraries are paired-end"),
        ("missing-genes", [{**rows[0], "reference_genes": "absent.bed"}], "annotation missing"),
        # Both rows use atac.fa.gz but name different organelles.
        ("organelles", [rows[0], {**rows[1], "organelle_contigs": "ChrM"}], "disagree on"),
    ):
        candidate = work / f"{suffix}.tsv"
        write_sheet(candidate, changed)
        cases.append((suffix, str(candidate), diagnostic))
    extra = work / "extra-column.tsv"
    extra.write_text("".join(line.rstrip("\n") + "\tnote\n" for line in lines))
    cases.append(("extra-column", str(extra), "expected header"))
    for suffix, candidate, diagnostic in cases:
        assert run_atac(work, Path(candidate), suffix=suffix, expect_failure=True).returncode != 0
        assert trace_rows(work, suffix) == []
        log = (work / f"output-{suffix}.log").read_text()
        assert diagnostic in log, (suffix, log[-2000:])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    engines = parser.add_mutually_exclusive_group()
    engines.add_argument(
        "--docker", action="store_true", help="Use real tools in the pinned Docker images"
    )
    engines.add_argument(
        "--conda", action="store_true", help="Use real tools from the environments/*.yaml files"
    )
    parser.add_argument(
        "--keep", action="store_true", help="Retain fixtures, logs, and work directories"
    )
    args = parser.parse_args()
    engine = "docker" if args.docker else "conda" if args.conda else None
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f"atac-{engine or 'regression'}-", dir=runs))
    passed = False
    try:
        (work / "inputs").mkdir()
        with serve(work / "inputs", docker=engine == "docker") as base_url:
            sheet = create_fixture(work, base_url, engine=engine)
            assert run_atac(work, sheet, engine=engine).returncode == 0
            check_outputs(work, engine=engine)
            check_resume(work, sheet, engine=engine)
            if not engine:
                check_failures(work, sheet)
                # Stub runs never publish into the references folder.
                assert sorted(p.name for p in (work / "references").iterdir()) == [
                    "atac.fa.gz",
                    "bare.fa.gz",
                    "genes.bed",
                ]
        label = f"Real-tool ({engine}) bulk ATAC" if engine else "Bulk ATAC Nextflow regression"
        print(f"{label} validation passed.", flush=True)
        passed = True
    finally:
        if args.keep or not passed:
            print(f"Artifacts retained: {work}", flush=True)
        else:
            shutil.rmtree(work)


if __name__ == "__main__":
    main()
