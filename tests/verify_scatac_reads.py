"""Regression checks for scATAC from raw 10x reads (workflows/scatac-reads.nf).

Without options: stub runs of the whole workflow and its sheet checks. With --docker: real chromap
and cell calling, in the pinned images, on simulated 10x reads (tests/scatac_reads_harness.nf);
the SRA download needs NCBI and is checked on real data instead.
"""

from __future__ import annotations

import argparse
import csv
import functools
import gzip
import hashlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Callable
from pathlib import Path

from verify_download import FAKE_MD5SUM, publish
from verify_pipeline import (
    available_resources,
    executor_limits,
    failed_task_errors,
    in_parallel,
    launch_folder,
    run_logged,
    serve,
    trace_rows,
)

ROOT = Path(__file__).resolve().parents[1]
CHROMAP_PIN = re.search(
    r"chromap=([^\s,<>]+)", (ROOT / "environments/GENESIS_CHROMAP.yaml").read_text()
)
assert CHROMAP_PIN, "environments/GENESIS_CHROMAP.yaml does not pin a chromap version"
CHROMAP_VERSION = CHROMAP_PIN.group(1)
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
    work: Path,
    command: list[str],
    suffix: str,
    *,
    expect_failure: bool = False,
    tools: Path | None = None,
    extra_environment: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run Nextflow in work; tools is a folder of stand-in executables put first on PATH."""
    first = f"{tools}:" if tools else ""
    environment = {
        **os.environ,
        **(extra_environment or {}),
        # Nextflow adds the main script's bin/ to PATH, and harnesses are in tests/.
        "PATH": f"{first}{ROOT / 'workflows/bin'}:{Path(sys.executable).parent}:"
        + os.environ["PATH"],
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
    }
    log = str(work / f"nextflow-{suffix}.log")
    full = ["nextflow", "-log", log, *command, "-ansi-log", "false"]
    print(f"Running scATAC reads ({suffix}) in {work}", flush=True)
    return run_logged(
        full,
        cwd=work,
        env=environment,
        log=work / f"output-{suffix}.log",
        timeout=3600 if "local,test,docker" in full else 300,
        expect_failure=expect_failure,
    )


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
        # A known 10x list missing from the references folder is downloaded (GET_WHITELIST).
        {
            "sample_id": "known",
            "species": "Plant_synthetic",
            "sra_runs": "SRR0000004",
            "chemistry": "10x-atac-v2",
            "cells_url": "-",
            "reference_fasta": "sc.fa.gz",
            "reference_genes": "genes.bed",
            "organelle_contigs": "ChrM",
            "barcode_whitelist": "737K-cratac-v1.txt",
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
        "GET_REFERENCES:BWA_MEM3_INDEX": 1,
        "DOWNLOAD_SCATAC_READS": 4,
        "GET_CHROMAP_INDEX:CHROMAP_INDEX": 1,
        "GET_WHITELIST": 1,
        "SCATAC_CHROMAP": 3,
        "SCATAC_CELLS": 3,
        "SCATAC_PSEUDOBULKS:SCATAC_FRAGMENTS": 3,
        "SCATAC_PSEUDOBULKS:SCATAC_PEAKS": 3,
        "SCATAC_PSEUDOBULKS:SCATAC_TRACKS": 3,
        "SCATAC_PSEUDOBULKS:SCATAC_QC": 3,
        "MULTIQC_REPORT": 1,
    }, counts
    # Runs are downloaded and converted one task each, then mapped per library.
    process = "DOWNLOAD_SCATAC_READS"
    names = sorted(row["name"] for row in rows if row["name"].startswith(f"{process} "))
    assert names == [
        f"{process} (known:SRR0000004)",
        f"{process} (lanes:SRR0000001)",
        f"{process} (lanes:SRR0000002)",
        f"{process} (multiome:SRR0000003)",
    ], names
    # The library with the known list maps with the downloaded copy.
    (task,) = [row for row in rows if row["name"] == "SCATAC_CHROMAP (known)"]
    staged = work / "work" / task["hash"]
    assert any(
        path.name == "737K-cratac-v1.txt" for path in staged.parent.glob(f"{staged.name}*/*")
    )
    # Stub runs never publish into the references folder.
    assert sorted(p.name for p in (work / "references").iterdir()) == [
        "atac.whitelist.txt",
        "genes.bed",
        "sc.fa.gz",
    ]


def chromap_record(references: Path, ref_id: str, version: str = CHROMAP_VERSION) -> str:
    """The <ref_id>.chromap.sha256 CHROMAP_INDEX writes: FASTA hash, release, index file sizes."""
    fasta = f"{ref_id}.fa.gz"
    lines = [f"{hashlib.sha256((references / fasta).read_bytes()).hexdigest()}  {fasta}"]
    lines.append(f"{version}  chromap")
    index = references / f"{ref_id}.chromap"
    lines += [
        f"{path.stat().st_size}  {index.name}/{path.name}"
        for path in sorted(index.iterdir())
        if path.is_file()
    ]
    return "\n".join(lines) + "\n"


def check_index_reuse(work: Path, sheet: Path) -> None:
    """A published chromap index is reused only while its record matches the FASTA and release."""

    def case(name: str, version: str, rebuilt: bool) -> None:
        folder = work / "index-reuse" / name
        references = folder / "references"
        references.mkdir(parents=True)
        for path in (work / "references").iterdir():
            shutil.copy2(path, references / path.name)
        (references / "sc.chromap").mkdir()
        (references / "sc.chromap/genome.index").write_text("index\n")
        (references / "sc.chromap.sha256").write_text(chromap_record(references, "sc", version))
        assert run_workflow(folder, sheet).returncode == 0, name
        names = [row["name"] for row in trace_rows(folder)]
        assert ("GET_CHROMAP_INDEX:CHROMAP_INDEX (sc)" in names) == rebuilt, (name, names)

    in_parallel(
        [
            functools.partial(case, "current", CHROMAP_VERSION, False),
            functools.partial(case, "other-release", "0.0.0", True),
        ]
    )


# Streams FASTQ records as fastq-dump [-N FIRST] -X LAST --split-spot -Z does from the stand-in
# .sra file (JSON: "lengths", one entry per read of the spot, "n" for every read n bases, or
# "first/rest"): each spot's reads in turn, skipping empty ones, and the spot count on stderr.
# Headers are as --defline-seq gives them ($si and $ri only), or fastq-dump's default, which has
# no read number.
FAKE_FASTQ_DUMP = """#!{python}
import json, sys
from pathlib import Path
args = sys.argv[1:]
option = lambda name, default: args[args.index(name) + 1] if name in args else default
specs = json.loads(Path(args[-1]).read_text())["lengths"]
first, last = int(option("-N", 1)), min(int(option("-X", 0)), {spots})
defline = option("--defline-seq", "@SRR1.$si name length=$rl")
for spot in range(first, last + 1):
    for number, spec in enumerate(specs, 1):
        head, _, rest = spec.partition("/")
        length = int(head) if spot == 1 else int(rest or head)
        if length:
            header = defline.replace("$si", str(spot)).replace("$ri", str(number))
            header = header.replace("$ac", "SRR1").replace("$sn", "name")
            header = header.replace("$rl", str(length))
            bases, qualities = "A" * length, "I" * length
            sys.stdout.write(f"{{header}}\\n{{bases}}\\n+\\n{{qualities}}\\n")
print(f"Read {{last - first + 1}} spots for SRR1", file=sys.stderr)
"""
# vdb-dump --info FILE: the stand-in run's spot count, as SEQ. vdb-dump --table_enum FILE: its
# tables, with PRIMARY_ALIGNMENT when the stand-in is "aligned".
FAKE_VDB_DUMP = """#!{python}
import json, sys
from pathlib import Path
if "--table_enum" in sys.argv:
    aligned = json.loads(Path(sys.argv[-1]).read_text()).get("aligned")
    print("tbl : PRIMARY_ALIGNMENT" if aligned else "tbl : SEQUENCE")
else:
    print("SEQ    : {spots}")
"""
# pigz -p N -4, stdin to stdout, as sra-split calls it.
FAKE_PIGZ = "#!/bin/sh\nexec gzip -c\n"
STAND_IN_SPOTS = 6


def check_read_selection(work: Path) -> None:
    """DOWNLOAD_SCATAC_READS picks the barcode and genomic reads by length, from the first spots,
    needs all three in every spot, and tries an unavailable run again."""
    tools = work / "read-selection" / "tools"
    tools.mkdir(parents=True)
    stand_ins = {
        "fastq-dump": FAKE_FASTQ_DUMP.format(python=sys.executable, spots=STAND_IN_SPOTS),
        "vdb-dump": FAKE_VDB_DUMP.format(python=sys.executable, spots=STAND_IN_SPOTS),
        "pigz": FAKE_PIGZ,
    }
    if shutil.which("md5sum") is None:
        stand_ins["md5sum"] = FAKE_MD5SUM.format(python=sys.executable)
    for name, text in stand_ins.items():
        (tools / name).write_text(text)
        (tools / name).chmod(0o755)
    root = work / "read-selection" / "server"
    root.mkdir()

    def case(base: str, name: str, lengths: str, expected: list[set[int]] | str) -> None:
        folder = work / "read-selection" / name
        folder.mkdir()
        run = f"SRR{int(hashlib.md5(name.encode()).hexdigest(), 16) % 10**7:07d}"
        contents = json.dumps({"lengths": lengths.split(), "aligned": name == "aligned"}).encode()
        link = f"{base}/sra/{run}"
        publish(root, base, link, contents)
        md5 = "0" * 32 if name == "unavailable" else hashlib.md5(contents).hexdigest()
        located = {
            "result": [
                {
                    "files": [
                        {
                            "type": "sra",
                            "size": len(contents),
                            "md5": md5,
                            "locations": [{"link": link}],
                        }
                    ]
                }
            ]
        }
        publish(root, base, f"{base}/sdl/{run}", json.dumps(located).encode())
        limits = folder / "limits.config"
        limits.write_text(executor_limits(docker=False))
        command = [
            "run",
            str(ROOT / "tests/download_scatac_reads_harness.nf"),
            "-c",
            str(ROOT / "conf/scatac-reads.config"),
            "-profile",
            "local,test",
            "-c",
            str(limits),
            "--run",
            run,
            "--barcode_length",
            "16",
            "--sra_locate_url",
            f"{base}/sdl/",
            "--run_dir",
            str(folder / "published"),
            "-work-dir",
            str(folder / "work"),
        ]
        result = nextflow(
            folder, command, "run", expect_failure=isinstance(expected, str), tools=tools
        )
        attempts = list(folder.glob("work/*/*/.command.sh"))
        if isinstance(expected, str):
            assert result.returncode != 0, name
            errors = failed_task_errors(folder / "work")
            assert expected in errors, (name, errors[-3000:])
            # Only an unavailable run is tried again: the first attempt and its two retries.
            assert len(attempts) == (3 if name == "unavailable" else 1), (name, attempts)
            return
        assert result.returncode == 0, name
        assert len(attempts) == 1, (name, attempts)
        observed = []
        for read in ("R1", "R2", "barcode"):
            (path,) = folder.glob(f"work/*/*/{run}.{read}.fastq.gz")
            with gzip.open(path, "rt") as stream:
                lines = list(stream)
            assert len(lines) == 4 * STAND_IN_SPOTS, (name, read, len(lines))
            observed.append({len(line.rstrip("\n")) for line in lines[1::4]})
        # Read lengths in R1, R2 and the barcode FASTQ.
        assert observed == expected, (name, observed)
        # The .sra file is deleted once converted, and nothing uncompressed is left.
        (task,) = (path.parent for path in attempts)
        assert not list(task.glob("*.sra")) and not list(task.glob("*.fastq")), name

    # Lengths per read of the spot, in SRA's order: index, R1, barcode, R2 as in 8/50/16/49.
    with serve(root, docker=False) as base:
        in_parallel(
            [
                functools.partial(case, base, "untrimmed", "8 50 16 49", [{50}, {49}, {16}]),
                # A trimmed first read no longer hides a genomic read: its longest read counts.
                functools.partial(
                    case, base, "trimmed-first-mate", "8 20/50 16 49", [{20, 50}, {49}, {16}]
                ),
                functools.partial(
                    case,
                    base,
                    "trimmed-barcode",
                    "8 50 15/16 49",
                    "has barcode reads of 15 to 16 bases",
                ),
                functools.partial(
                    case, base, "no-barcode", "8 50 49", "has no single 16-base barcode read"
                ),
                # A barcode read only in the first spot: the others lack it.
                functools.partial(
                    case, base, "missing-barcode", "8 50 16/0 49", "read 3 in 1 of 6 spots"
                ),
                # Loaded from BAM: no barcode read to find.
                functools.partial(
                    case, base, "aligned", "8 50 16 49", "stores aligned reads (loaded from BAM)"
                ),
                # An MD5 that is never NCBI's: tried again, then the run stops.
                functools.partial(
                    case, base, "unavailable", "8 50 16 49", "SRA run unavailable (the .sra"
                ),
            ]
        )


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
        ("sci-list", [{**rows[0], "chemistry": "sci-atac"}], "set barcode_whitelist to -"),
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
    """10x and sci-ATAC reads for three libraries; returns each library's expected fragments."""
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
    # sci-ATAC: a fixed 6-base sample index, then a combinatorial cell index; the 28-base read is
    # the barcode as is, with no list to correct it against.
    sci = ["ATCACG" + "".join(generator.choices("ACGT", k=22)) for _ in range(9)]
    for library, multiome in (("v1", False), ("arc", True), ("sci", False)):
        r1, r2, codes = [], [], []
        fragments: Counter = Counter()
        cells = {"v1": whitelist[:6], "arc": whitelist[10:16], "sci": sci[:6]}[library]
        noise = {"v1": whitelist[30:33], "arc": whitelist[40:43], "sci": sci[6:]}[library]
        for barcode, count in [(b, 40) for b in cells] + [(b, 2) for b in noise]:
            chosen = []
            for _ in range(count):
                contig = "ChrM" if generator.random() < 0.1 else "chr1"
                start = generator.randrange(1000, len(genome[contig]) - 1000)
                chosen.append((contig, start, start + generator.randrange(150, 400)))
            # Some fragments are read twice: PCR duplicates, kept once.
            for contig, start, end in chosen + chosen[:5]:
                code = barcode
                if library != "sci" and generator.random() < 0.2:
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
    # chromap's --low-mem mode (part of --preset atac) allocates over 10 GB to merge its temporary
    # files, however few reads there are: a fixed 10 GiB buffer for BED output (chromap 0.3.2,
    # src/mapping_writer.h). So SCATAC_CHROMAP gets 12 GB, above the test profile's 2 GB cap,
    # and the executor at least as much, so that it can run.
    chromap = 12 * 2**30
    _cpus, memory = available_resources(docker=True)
    assert memory > chromap, f"chromap needs 12 GB; Docker has {memory / 2**30:.1f} GB"
    limits.write_text(
        executor_limits(docker=True)
        + f"executor.memory = '{max(chromap, memory // 2) // 2**20} MB'\n"
        + "process {\n"
        + "    withName: 'SCATAC_CHROMAP' {\n"
        + "        resourceLimits = [cpus: 1, memory: 12.GB, time: 5.m]\n"
        + "    }\n"
        + "}\n"
        # Linux Docker only resolves host.docker.internal, where the cell lists are served, with
        # an explicit host-gateway mapping.
        + "docker.runOptions = '-u $(id -u):$(id -g) "
        "--add-host=host.docker.internal:host-gateway'\n"
    )

    def harness(
        suffix: str,
        libraries: str = "v1:16:0:whitelist.txt,arc:24:8:whitelist.txt,sci:28:0:-",
        *,
        expect_failure: bool = False,
    ) -> list[str]:
        # The FASTA's folder is the references folder, where the chromap index is published.
        result = nextflow(
            work,
            [
                "run",
                str(ROOT / "tests/scatac_reads_harness.nf"),
                "-c",
                str(ROOT / "conf/scatac-reads.config"),
                "-profile",
                "local,test,docker",
                "-c",
                str(limits),
                "--inputs",
                str(work / "inputs"),
                "--fasta",
                str(work / "sc.fa.gz"),
                "--references",
                str(work),
                "--libraries",
                libraries,
                "--run_dir",
                str(work / f"published-{suffix}"),
                "--scatac_min_fragments",
                "30",
                "-work-dir",
                str(work / ("work" if suffix == "real" else f"work-{suffix}")),
                "-with-trace",
                str(work / f"trace-{suffix}.tsv"),
            ],
            suffix,
            expect_failure=expect_failure,
        )
        assert (result.returncode != 0) == expect_failure, (suffix, result.returncode)
        return [row["name"] for row in trace_rows(work, suffix)]

    assert "GET_CHROMAP_INDEX:CHROMAP_INDEX (sc)" in harness("real")
    # The index is published with a record of the FASTA, the pinned release and its size.
    record = work / "sc.chromap.sha256"
    assert record.read_text() == chromap_record(work, "sc"), record.read_text()
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
        # Exact 10x fragments: Tn5 offsets applied, duplicates merged, barcodes corrected (sci-ATAC
        # barcodes kept whole).
        assert observed == fragments, (library, len(observed), len(fragments))
        cells = sorted(work.glob(f"work/*/*/{library}.source.cells.tsv"))
        rows = [line.split("\t") for line in cells[0].read_text().splitlines()[1:]]
        called = {barcode for barcode, group in rows}
        nuclear = Counter(b for (contig, _s, _e, b) in fragments if contig != "ChrM")
        assert called == {b for b, n in nuclear.items() if n >= 30}, (library, called)
        assert {group for _b, group in rows} == {"all"}
    # A later run reuses the published index; a record from another release rebuilds it.
    assert "GET_CHROMAP_INDEX:CHROMAP_INDEX (sc)" not in harness("reuse")
    record.write_text(record.read_text().replace(f"{CHROMAP_VERSION}  chromap", "0.0.0  chromap"))
    assert "GET_CHROMAP_INDEX:CHROMAP_INDEX (sc)" in harness("rebuild")
    assert record.read_text() == chromap_record(work, "sc")
    check_cell_lists(work, expected["v1"], harness)


def check_cell_lists(work: Path, fragments: Counter, harness: Callable[..., list[str]]) -> None:
    """Authors' cell lists for the v1 library: its cells pass; another orientation or library's
    list stops the library."""
    nuclear = Counter(b for (contig, _s, _e, b) in fragments if contig != "ChrM")
    cells = sorted(b for b, n in nuclear.items() if n >= 30)
    background = sorted(b for b, n in nuclear.items() if n < 30)
    assert len(cells) > len(background) > 1, (cells, background)
    lists = work / "lists"
    lists.mkdir()
    cases = {
        # Cell Ranger's -1 suffix is dropped; a barcode of another library is not found.
        "listed": ([f"{b}-1" for b in cells] + ["ATCACG" + "A" * 22], None),
        "flipped": (
            [reverse_complement(b) for b in cells],
            "No listed barcode is in the fragments",
        ),
        "other": (background + cells[:1], "Most listed barcodes found have too few fragments"),
    }
    for name in cases:
        for suffix in ("R1", "R2", "barcode"):
            shutil.copyfile(
                work / "inputs" / f"v1.{suffix}.fastq.gz",
                work / "inputs" / f"{name}.{suffix}.fastq.gz",
            )
    with serve(lists, docker=True) as base_url:
        for name, (listed, diagnostic) in cases.items():
            (lists / f"{name}.tsv").write_text(
                "barcode\tgroup\n" + "".join(f"{b}\tg\n" for b in listed)
            )
            harness(
                f"cells-{name}",
                f"{name}:16:0:whitelist.txt:{base_url}/{name}.tsv",
                expect_failure=diagnostic is not None,
            )
            folder = work / f"work-cells-{name}"
            if diagnostic is not None:
                assert diagnostic in failed_task_errors(folder), name
                continue
            [published] = [
                p for p in folder.glob(f"*/*/{name}.source.cells.tsv") if not p.is_symlink()
            ]
            rows = [line.split("\t") for line in published.read_text().splitlines()[1:]]
            assert sorted(b for b, _group in rows) == sorted(b.removesuffix("-1") for b in listed)
            log = (published.parent / ".command.err").read_text()
            found = f"{len(cells)} of {len(listed)} listed barcodes found, {len(cells)} of them"
            assert found in log, log


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--docker",
        action="store_true",
        help="Also run chromap and cell calling in the pinned images",
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
        check_index_reuse(work, sheet)
        check_read_selection(work)
        check_failures(work, sheet)
        if args.docker:
            real = work / "real"
            real.mkdir()
            check_real_tools(real)
        label = "scATAC reads regression" + (" and real-tool chromap" if args.docker else "")
        print(f"PASS: {label}", flush=True)
        passed = True
    finally:
        if passed and not args.keep:
            shutil.rmtree(work)
        else:
            print(f"Kept {work}", flush=True)


if __name__ == "__main__":
    main()
