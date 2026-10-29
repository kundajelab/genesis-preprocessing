"""Regression checks for scATAC from raw 10x reads (workflows/scatac-reads.nf).

Without options: stub runs of the whole workflow and its sheet checks. With --conda: real chromap
and cell calling on simulated 10x reads (tests/scatac_reads_harness.nf); the SRA download needs
NCBI and is checked on real data instead.
"""

from __future__ import annotations

import argparse
import csv
import functools
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from verify_pipeline import executor_limits, in_parallel, launch_folder, trace_rows

ROOT = Path(__file__).resolve().parents[1]
HEADER = [
    "sample_id",
    "species",
    "sra_runs",
    "chemistry",
    "cells_url",
    "reference_fasta",
    "reference_genes",
    "organelle_contigs",
    "barcode_whitelist",
]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def write_sheet(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def nextflow(
    work: Path, command: list[str], suffix: str, *, expect_failure: bool = False
) -> subprocess.CompletedProcess[str]:
    environment = {
        **os.environ,
        "PATH": f"{Path(sys.executable).parent}:{os.environ['PATH']}",
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
    }
    log = str(work / f"nextflow-{suffix}.log")
    full = ["nextflow", "-log", log, *command, "-ansi-log", "false"]
    print(f"Running scATAC reads ({suffix}) in {work}", flush=True)
    result = subprocess.run(
        full,
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


def run_workflow(
    work: Path, sheet: Path, *, expect_failure: bool = False
) -> subprocess.CompletedProcess[str]:
    """A stub run of the whole workflow."""
    limits = work / "limits.config"
    if not limits.exists():
        limits.write_text(executor_limits(docker=False))
    command = [
        "run",
        str(ROOT / "workflows/scatac-reads.nf"),
        "-c",
        str(ROOT / "conf/scatac-reads.config"),
        "-profile",
        "local,test",
        "--input",
        str(sheet),
        "--references",
        str(work / "references"),
        "--workspace",
        str(work / "published"),
        "-work-dir",
        str(work / "work"),
        "-c",
        str(limits),
        "-with-trace",
        str(work / "trace-run.tsv"),
        "-stub-run",
    ]
    return nextflow(work, command, "run", expect_failure=expect_failure)


def create_stub_fixture(work: Path) -> Path:
    references = work / "references"
    references.mkdir()
    with gzip.open(references / "sc.fa.gz", "wt") as stream:
        stream.write(">chr1\n" + "ACGT" * 2500 + "\n>ChrM\n" + "ACGT" * 500 + "\n")
    (references / "genes.bed").write_text("chr1\t100\t2000\tgene0\t0\t+\n")
    (references / "atac.whitelist.txt").write_text("ACGTACGTACGTACGT\n")
    rows = [
        {
            "sample_id": "lanes",
            "species": "Plant_synthetic",
            "sra_runs": "SRR0000001,SRR0000002",
            "chemistry": "10x-atac-v1",
            "cells_url": "-",
            "reference_fasta": "sc.fa.gz",
            "reference_genes": "genes.bed",
            "organelle_contigs": "ChrM",
            "barcode_whitelist": "atac.whitelist.txt",
        },
        {
            "sample_id": "multiome",
            "species": "Plant_synthetic",
            "sra_runs": "SRR0000003",
            "chemistry": "10x-multiome-atac",
            "cells_url": "https://example.org/cells.tsv",
            "reference_fasta": "sc.fa.gz",
            "reference_genes": "genes.bed",
            "organelle_contigs": "ChrM",
            "barcode_whitelist": "-",
        },
    ]
    sheet = work / "samples.tsv"
    write_sheet(sheet, rows)
    return sheet


def check_stub(work: Path) -> None:
    rows = trace_rows(work)
    assert all(row["status"] == "COMPLETED" for row in rows), rows
    counts = Counter(row["name"].split(" ")[0] for row in rows)
    assert counts == {
        "GET_REFERENCES:CHROM_SIZES": 1,
        "GET_REFERENCES:BWA_MEM2_INDEX": 1,
        "SCATAC_SRA:SRA_PREFETCH": 3,
        "SCATAC_SRA:SRA_FASTQ": 3,
        "SCATAC_CHROMAP": 2,
        "SCATAC_CELLS": 2,
        "SCATAC_PSEUDOBULKS:SCATAC_FRAGMENTS": 2,
        "SCATAC_PSEUDOBULKS:SCATAC_PEAKS": 2,
        "SCATAC_PSEUDOBULKS:SCATAC_TRACKS": 2,
        "SCATAC_PSEUDOBULKS:SCATAC_QC": 2,
        "MULTIQC_REPORT": 1,
    }, counts
    # Runs are downloaded and converted one task each, then mapped per library.
    for process in ("SCATAC_SRA:SRA_PREFETCH", "SCATAC_SRA:SRA_FASTQ"):
        names = sorted(row["name"] for row in rows if row["name"].startswith(f"{process} "))
        assert names == [
            f"{process} (lanes:SRR0000001)",
            f"{process} (lanes:SRR0000002)",
            f"{process} (multiome:SRR0000003)",
        ], names
    # Stub runs never publish into the references folder.
    assert sorted(p.name for p in (work / "references").iterdir()) == [
        "atac.whitelist.txt",
        "genes.bed",
        "sc.fa.gz",
    ]


def check_failures(work: Path, sheet: Path) -> None:
    """A rejected sheet stops before any task."""
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    lines = sheet.read_text().splitlines(keepends=True)
    cases: list[tuple[str, str, str]] = []
    for suffix, changed, diagnostic in (
        ("chemistry", [{**rows[0], "chemistry": "10x-rna-v3"}], "chemistry must be"),
        ("run-format", [{**rows[0], "sra_runs": "SRR1;SRR2"}], "sra_runs must be"),
        ("whitelist", [{**rows[0], "barcode_whitelist": "absent.txt"}], "whitelist missing"),
        ("cells", [{**rows[0], "cells_url": "cells.tsv"}], "cells_url must be"),
        ("shared-run", [rows[0], {**rows[1], "sra_runs": "SRR0000002"}], "is listed for both"),
        ("organelles", [rows[0], {**rows[1], "organelle_contigs": "-"}], "disagree on"),
    ):
        candidate = work / f"{suffix}.tsv"
        write_sheet(candidate, changed)
        cases.append((suffix, str(candidate), diagnostic))
    extra = work / "extra-column.tsv"
    extra.write_text("".join(line.rstrip("\n") + "\tnote\n" for line in lines))
    cases.append(("extra-column", str(extra), "expected header"))

    def reject(suffix: str, candidate: str, diagnostic: str) -> None:
        case = launch_folder(work, suffix)
        assert run_workflow(case, Path(candidate), expect_failure=True).returncode != 0
        assert trace_rows(case) == [], suffix
        log = (case / "output-run.log").read_text()
        assert diagnostic in log, (suffix, log[-2000:])

    in_parallel([functools.partial(reject, *case) for case in cases])


def fastq(path: Path, reads: list[str]) -> None:
    with gzip.open(path, "wt") as stream:
        for index, read in enumerate(reads):
            stream.write(f"@r{index}\n{read}\n+\n{'I' * len(read)}\n")


def simulate(work: Path, generator: random.Random) -> dict[str, Counter]:
    """10x reads for two libraries; returns each library's expected chromap fragments."""
    genome = {
        "chr1": "".join(generator.choices("ACGT", k=1000000)),
        "ChrM": "".join(generator.choices("ACGT", k=20000)),
    }
    with gzip.open(work / "sc.fa.gz", "wt") as stream:
        stream.writelines(f">{name}\n{sequence}\n" for name, sequence in genome.items())
    whitelist = ["".join(generator.choices("ACGT", k=16)) for _ in range(60)]
    (work / "inputs" / "whitelist.txt").write_text("".join(f"{b}\n" for b in whitelist))
    expected: dict[str, Counter] = {}
    # 10x ATAC v1: a 16-base barcode read as is. Multiome: an 8-base spacer, then the barcode,
    # read on the other strand.
    for library, multiome in (("v1", False), ("arc", True)):
        r1, r2, codes = [], [], []
        fragments: Counter = Counter()
        cells = whitelist[:6] if not multiome else whitelist[10:16]
        noise = whitelist[30:33] if not multiome else whitelist[40:43]
        for barcode, count in [(b, 40) for b in cells] + [(b, 2) for b in noise]:
            chosen = []
            for _ in range(count):
                contig = "ChrM" if generator.random() < 0.1 else "chr1"
                start = generator.randrange(1000, len(genome[contig]) - 1000)
                chosen.append((contig, start, start + generator.randrange(150, 400)))
            # Some fragments are read twice: PCR duplicates, kept once.
            for contig, start, end in chosen + chosen[:5]:
                code = barcode
                if generator.random() < 0.2:
                    # One sequencing error, corrected against the whitelist.
                    position = generator.randrange(16)
                    base = generator.choice([b for b in "ACGT" if b != code[position]])
                    code = code[:position] + base + code[position + 1 :]
                r1.append(genome[contig][start : start + 50])
                r2.append(reverse_complement(genome[contig][end - 50 : end]))
                read = code if not multiome else "".join(generator.choices("ACGT", k=8)) + code
                codes.append(read if not multiome else read[:8] + reverse_complement(read[8:]))
                fragments[(contig, start + 4, end - 5, barcode)] += 1
        fastq(work / "inputs" / f"{library}.R1.fastq.gz", r1)
        fastq(work / "inputs" / f"{library}.R2.fastq.gz", r2)
        fastq(work / "inputs" / f"{library}.barcode.fastq.gz", codes)
        expected[library] = fragments
    return expected


def check_real_tools(work: Path) -> None:
    """chromap and cell calling with real tools on simulated 10x reads."""
    (work / "inputs").mkdir()
    expected = simulate(work, random.Random(1555))
    limits = work / "limits.config"
    limits.write_text(
        executor_limits(docker=False) + f"conda.cacheDir = '{work.parent / 'conda_cache'}'\n"
    )
    command = [
        "run",
        str(ROOT / "tests/scatac_reads_harness.nf"),
        "-c",
        str(ROOT / "conf/scatac-reads.config"),
        "-profile",
        "local,test,conda",
        "-c",
        str(limits),
        "--inputs",
        str(work / "inputs"),
        "--fasta",
        str(work / "sc.fa.gz"),
        "--libraries",
        "v1:16:whitelist.txt,arc:24:whitelist.txt",
        "--run_dir",
        str(work / "published"),
        "--scatac_min_fragments",
        "30",
        "-work-dir",
        str(work / "work"),
    ]
    assert nextflow(work, command, "real").returncode == 0
    for library, fragments in expected.items():
        # Each output once: later tasks stage it as a symlink.
        work_files = [
            p
            for p in sorted(work.glob(f"work/*/*/{library}.source.fragments.gz"))
            if not p.is_symlink()
        ]
        assert len(work_files) == 1, work_files
        with gzip.open(work_files[0], "rt") as stream:
            observed = Counter()
            for line in stream:
                contig, start, end, barcode, count = line.rstrip("\n").split("\t")
                observed[(contig, int(start), int(end), barcode)] = int(count)
        # Exact 10x fragments: Tn5 offsets applied, duplicates merged, barcodes corrected.
        assert observed == fragments, (library, len(observed), len(fragments))
        cells = sorted(work.glob(f"work/*/*/{library}.source.cells.tsv"))
        rows = [line.split("\t") for line in cells[0].read_text().splitlines()[1:]]
        called = {barcode for barcode, group in rows}
        nuclear = Counter(b for (contig, _s, _e, b) in fragments if contig != "ChrM")
        assert called == {b for b, n in nuclear.items() if n >= 30}, (library, called)
        assert {group for _b, group in rows} == {"all"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--conda", action="store_true", help="Also run chromap and cell calling with real tools"
    )
    parser.add_argument(
        "--keep", action="store_true", help="Retain fixtures, logs, and work directories"
    )
    args = parser.parse_args()
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="scatac-reads-", dir=runs))
    passed = False
    try:
        sheet = create_stub_fixture(work)
        assert run_workflow(work, sheet).returncode == 0
        check_stub(work)
        check_failures(work, sheet)
        if args.conda:
            real = work / "real"
            real.mkdir()
            check_real_tools(real)
        label = "scATAC reads regression" + (" and real-tool chromap" if args.conda else "")
        print(f"PASS: {label}", flush=True)
        passed = True
    finally:
        if passed and not args.keep:
            shutil.rmtree(work)
        else:
            print(f"Kept {work}", flush=True)


if __name__ == "__main__":
    main()
