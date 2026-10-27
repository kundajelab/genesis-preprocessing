"""scATAC fragment-workflow regressions and optional real-tool validation; no test framework."""

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
from verify_pipeline import serve, trace_rows

ROOT = Path(__file__).resolve().parents[1]
HEADER = [
    "sample_id",
    "species",
    "fragments_url",
    "cells_url",
    "reference_fasta",
    "reference_genes",
    "organelle_contigs",
]
# Each group's cells are accessible at their own planted sites; genes start at each site.
SITES = {
    "A": [50000 + 40000 * index for index in range(5)],
    "B": [450000 + 40000 * index for index in range(5)],
}
CELLS = {"A": 20, "B": 10}


def write_sheet(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def library(generator: random.Random, sample: str) -> tuple[list[tuple], dict[str, str]]:
    """10x-style rows: listed cells with site, background and organelle fragments, plus noise."""
    cells = {
        f"{sample.upper()}{group}{index:03d}-1": group
        for group in CELLS
        for index in range(CELLS[group])
    }
    # One row per cell and interval, as producers write them.
    rows: dict[tuple[str, int, int, str], int] = {}
    for barcode, group in cells.items():
        for _ in range(generator.randrange(150, 250)):
            roll = generator.random()
            if roll < 0.5:
                center = generator.choice(SITES[group]) + generator.randrange(-100, 100)
                contig, start = "chr1", center - generator.randrange(40, 120)
            elif roll < 0.6:
                # Low, even accessibility around the group's sites gives the TSS flanks some cuts.
                center = generator.choice(SITES[group]) + generator.randrange(-3000, 3000)
                contig, start = "chr1", center - generator.randrange(40, 120)
            elif roll < 0.7:
                contig, start = "ChrM", generator.randrange(19000)
            else:
                contig, start = "chr1", generator.randrange(1000, 990000)
            end = start + generator.randrange(60, 400)
            rows[(contig, start, end, barcode)] = generator.randrange(1, 4)
    # Barcodes that were not called as cells.
    for index in range(200):
        start = generator.randrange(1000, 990000)
        rows[("chr1", start, start + 150, f"NOISE{index:03d}-1")] = 1
    return sorted((*key, support) for key, support in rows.items()), cells


def create_fixture(work: Path, base_url: str) -> Path:
    generator = random.Random(1555)
    references = work / "references"
    references.mkdir()
    with gzip.open(references / "sc.fa.gz", "wt") as stream:
        stream.write(">chr1\n" + "".join(generator.choices("ACGT", k=1000000)) + "\n")
        stream.write(">ChrM\n" + "".join(generator.choices("ACGT", k=20000)) + "\n")
    (references / "genes.bed").write_text(
        "".join(
            f"chr1\t{site}\t{site + 2000}\tgene{index}\t0\t+\n"
            if index % 2 == 0
            else f"chr1\t{site - 2000}\t{site}\tgene{index}\t0\t-\n"
            for index, site in enumerate(SITES["A"] + SITES["B"])
        )
    )
    rows = []
    inputs = work / "inputs"
    for sample in ("rep1", "rep2"):
        fragments, cells = library(generator, sample)
        # Producers write a comment header; the order here is the producer's, not sorted by barcode.
        lines = ["# id=synthetic\n"] + ["\t".join(map(str, row)) + "\n" for row in fragments]
        with gzip.open(inputs / f"{sample}.fragments.tsv.gz", "wt") as stream:
            stream.writelines(lines)
        (inputs / f"{sample}.cells.tsv").write_text(
            "barcode\tgroup\n"
            + "".join(f"{barcode}\t{group}\n" for barcode, group in cells.items())
        )
        rows.append(
            {
                "sample_id": sample,
                "species": "Plant_synthetic",
                "fragments_url": f"{base_url}/{sample}.fragments.tsv.gz",
                "cells_url": f"{base_url}/{sample}.cells.tsv",
                "reference_fasta": "sc.fa.gz",
                "reference_genes": "genes.bed",
                "organelle_contigs": "ChrM",
            }
        )
    sheet = work / "samples.tsv"
    write_sheet(sheet, rows)
    return sheet


def run_scatac(
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
    limits = work / "limits.config"
    if not limits.exists():
        text = "executor.cpus = 4\nexecutor.memory = '8 GB'\n"
        if engine == "docker":
            # Linux Docker only resolves host.docker.internal with an explicit host-gateway mapping.
            text += (
                "docker.runOptions = '-u $(id -u):$(id -g) "
                "--add-host=host.docker.internal:host-gateway'\n"
            )
        if engine == "conda":
            # Keep the environments built for earlier runs.
            text += f"conda.cacheDir = '{work.parent / 'conda_cache'}'\n"
        limits.write_text(text)
    command = [
        "nextflow",
        "-log",
        str(work / f"nextflow-{suffix}.log"),
        "run",
        str(ROOT / "workflows/scatac-seq.nf"),
        "-c",
        str(ROOT / "conf/scatac-seq.config"),
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
        str(limits),
        "-with-trace",
        str(work / f"trace-{suffix}.tsv"),
        "-ansi-log",
        "false",
    ]
    if not engine:
        command.append("-stub-run")
    if resume:
        command.append("-resume")
    print(f"Running {engine or 'stub'} scATAC ({suffix}) in {work}", flush=True)
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


def read_table(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def check_outputs(work: Path, *, engine: str | None) -> None:
    """Assert the process wiring, publication and recorded provenance."""
    rows = trace_rows(work)
    assert all(row["status"] == "COMPLETED" for row in rows), rows
    counts = Counter(row["name"].split(" ")[0] for row in rows)
    # Stub libraries have a single group; real ones have A and B.
    groups = 4 if engine else 2
    assert counts == {
        "GET_REFERENCES:CHROM_SIZES": 1,
        "GET_REFERENCES:BWA_MEM2_INDEX": 1,
        "SCATAC_DOWNLOAD": 2,
        "SCATAC_FRAGMENTS": 2,
        "SCATAC_PEAKS": groups,
        "SCATAC_TRACKS": groups,
        "SCATAC_QC": groups,
        "MULTIQC_REPORT": 1,
    }, counts
    published = work / "published" / "samples" / "output"
    for path in published.rglob("*"):
        assert not path.name.endswith((".bam", ".bai", ".fragments.bed")), path
    for sample in ("rep1", "rep2"):
        folder = published / "Plant_synthetic" / sample
        for name in (
            f"{sample}.fragments.tsv.gz",
            f"{sample}.fragments.tsv.gz.tbi",
            f"qc/{sample}.cells.tsv",
            f"qc/{sample}.scatac_library_mqc.tsv",
        ):
            assert (folder / name).is_file(), folder / name
        for group in ("A", "B") if engine else ("all",):
            pseudobulk = folder / "groups" / group
            for name in (
                f"{sample}.{group}.macs3_peaks.narrowPeak.gz",
                f"{sample}.{group}.cuts.counts.bw",
                f"{sample}.{group}.scatac_qc_mqc.tsv",
            ):
                assert (pseudobulk / name).is_file(), pseudobulk / name
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
    for sample in ("rep1", "rep2"):
        with gzip.open(work / "inputs" / f"{sample}.fragments.tsv.gz", "rt") as stream:
            source = [line.rstrip("\n").split("\t") for line in stream if not line.startswith("#")]
        groups = dict(
            line.split("\t")
            for line in (work / "inputs" / f"{sample}.cells.tsv").read_text().splitlines()[1:]
        )
        kept = sorted(
            (row for row in source if row[3] in groups and row[0] != "ChrM"),
            key=lambda row: (row[0], int(row[1]), int(row[2]), row[3]),
        )
        folder = published / "Plant_synthetic" / sample
        with gzip.open(folder / f"{sample}.fragments.tsv.gz", "rt") as stream:
            written = [line.rstrip("\n").split("\t") for line in stream]
        # Listed cells on nuclear contigs only, unchanged coordinates and support, sorted.
        assert written == kept, (len(written), len(kept))
        cells = {row["barcode"]: row for row in read_table(folder / f"qc/{sample}.cells.tsv")}
        assert set(cells) == set(groups)
        for barcode, row in cells.items():
            mine = [fields for fields in kept if fields[3] == barcode]
            organelle = sum(1 for fields in source if fields[3] == barcode and fields[0] == "ChrM")
            assert row["group"] == groups[barcode], row
            assert int(row["fragments"]) == len(mine), row
            assert int(row["read_pairs"]) == sum(int(fields[4]) for fields in mine), row
            assert int(row["organelle_fragments"]) == organelle, row
        (summary,) = read_table(folder / f"qc/{sample}.scatac_library_mqc.tsv")
        assert int(summary["source_fragments"]) == len(source), summary
        assert int(summary["unlisted_barcode_fragments"]) == sum(
            1 for row in source if row[3] not in groups
        )
        assert int(summary["cell_fragments"]) == len(kept), summary
        for group, sites in SITES.items():
            pseudobulk = folder / "groups" / group
            (qc,) = read_table(pseudobulk / f"{sample}.{group}.scatac_qc_mqc.tsv")
            fragments = [fields for fields in kept if groups[fields[3]] == group]
            assert int(qc["cells"]) == CELLS[group], qc
            assert int(qc["fragments"]) == len(fragments), qc
            # Two Tn5 insertions per fragment, most at the group's own sites.
            assert int(qc["insertions"]) == 2 * len(fragments), qc
            assert float(qc["FRiP"]) > 0.4, qc
            assert float(qc["TSS_enrichment"]) > 3, qc
            with gzip.open(
                pseudobulk / f"{sample}.{group}.macs3_peaks.narrowPeak.gz", "rt"
            ) as stream:
                peaks = [line.split("\t") for line in stream]
            found = {
                site
                for site in sites
                for fields in peaks
                if int(fields[1]) <= site < int(fields[2])
            }
            assert len(found) >= 4, (sample, group, sorted(found))
            other = [site for name, places in SITES.items() if name != group for site in places]
            assert not any(
                int(fields[1]) <= site < int(fields[2]) for site in other for fields in peaks
            )


def check_cut_sites(work: Path) -> None:
    """Each fragment [s, e) gives Tn5 cuts at bases s and e - 1, the bases bulk ATAC counts."""
    tasks = [row for row in trace_rows(work) if row["name"].startswith("SCATAC_TRACKS ")]
    assert len(tasks) == 4, tasks
    for row in tasks:
        folder = work / "work" / row["hash"]
        (folder,) = folder.parent.glob(folder.name + "*")
        (fragments,) = folder.glob("*.fragments.bed")
        expected: Counter[tuple[str, int]] = Counter()
        for line in fragments.read_text().splitlines():
            contig, start, end = line.split("\t")[:3]
            expected.update([(contig, int(start)), (contig, int(end) - 1)])
        (bam,) = folder.glob("*.cuts.bam")
        with pysam.AlignmentFile(str(bam)) as stream:
            cuts = Counter((read.reference_name, read.reference_start) for read in stream)
        assert cuts == expected, row["name"]


def check_resume(work: Path, sheet: Path, *, engine: str | None) -> None:
    """Resumed runs reuse the downloads, and then every task."""
    assert run_scatac(work, sheet, engine=engine, resume=True).returncode == 0
    rows = trace_rows(work, "resume")
    tasks = [row for row in rows if row["name"].split(" ")[0] == "SCATAC_DOWNLOAD"]
    assert tasks and all(row["status"] == "CACHED" for row in tasks), tasks
    # Newly published chromosome sizes change staged input paths on the first resume.
    # A second resume must cache all downstream computation at those stable paths.
    assert run_scatac(work, sheet, engine=engine, resume=True, suffix="stable").returncode == 0
    rows = trace_rows(work, "stable")
    assert rows and all(row["status"] == "CACHED" for row in rows), rows


def check_rejected_fragments(work: Path, sheet: Path, base_url: str, *, engine: str) -> None:
    """Cell barcodes that match no fragment, e.g. without the -1 suffix, stop the library."""
    inputs = work / "inputs"
    cells = (inputs / "rep1.cells.tsv").read_text()
    (inputs / "suffix.cells.tsv").write_text(cells.replace("-1\t", "\t"))
    with sheet.open() as stream:
        (row, _) = list(csv.DictReader(stream, delimiter="\t"))
    candidate = work / "suffix.tsv"
    changed = {"cells_url": f"{base_url}/suffix.cells.tsv", "sample_id": "suffix"}
    write_sheet(candidate, [{**row, **changed}])
    result = run_scatac(work, candidate, engine=engine, suffix="suffix", expect_failure=True)
    assert result.returncode != 0
    assert "No fragments for any listed barcode" in result.stdout, result.stdout[-3000:]


def check_duplicate_barcodes(work: Path, sheet: Path, base_url: str, *, engine: str) -> None:
    """A barcode listed twice, even in two groups, stops the library."""
    inputs = work / "inputs"
    lines = (inputs / "rep1.cells.tsv").read_text().splitlines(keepends=True)
    barcode = lines[1].split("\t")[0]
    (inputs / "twice.cells.tsv").write_text("".join(lines) + f"{barcode}\tB\n")
    with sheet.open() as stream:
        (row, _) = list(csv.DictReader(stream, delimiter="\t"))
    candidate = work / "twice.tsv"
    changed = {"cells_url": f"{base_url}/twice.cells.tsv", "sample_id": "twice"}
    write_sheet(candidate, [{**row, **changed}])
    result = run_scatac(work, candidate, engine=engine, suffix="twice", expect_failure=True)
    assert result.returncode != 0
    message = f"barcode {barcode} is listed again (first on line 2)"
    assert message in result.stdout, result.stdout[-3000:]


def check_failures(work: Path, sheet: Path) -> None:
    """A rejected sheet stops before any task."""
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    lines = sheet.read_text().splitlines(keepends=True)
    cases: list[tuple[str, str, str]] = []
    for suffix, changed, diagnostic in (
        ("not-url", [{**rows[0], "cells_url": "cells.tsv"}], "cells_url must be an HTTP"),
        ("missing-genes", [{**rows[0], "reference_genes": "absent.bed"}], "annotation missing"),
        ("organelles", [rows[0], {**rows[1], "organelle_contigs": "-"}], "disagree on"),
    ):
        candidate = work / f"{suffix}.tsv"
        write_sheet(candidate, changed)
        cases.append((suffix, str(candidate), diagnostic))
    extra = work / "extra-column.tsv"
    extra.write_text("".join(line.rstrip("\n") + "\tnote\n" for line in lines))
    cases.append(("extra-column", str(extra), "expected header"))
    for suffix, candidate, diagnostic in cases:
        assert run_scatac(work, Path(candidate), suffix=suffix, expect_failure=True).returncode != 0
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
    work = Path(tempfile.mkdtemp(prefix=f"scatac-{engine or 'regression'}-", dir=runs))
    passed = False
    try:
        (work / "inputs").mkdir()
        with serve(work / "inputs", docker=engine == "docker") as base_url:
            sheet = create_fixture(work, base_url)
            assert run_scatac(work, sheet, engine=engine).returncode == 0
            check_outputs(work, engine=engine)
            if engine:
                check_cut_sites(work)
            check_resume(work, sheet, engine=engine)
            if engine:
                check_rejected_fragments(work, sheet, base_url, engine=engine)
                check_duplicate_barcodes(work, sheet, base_url, engine=engine)
            else:
                check_failures(work, sheet)
                # Stub runs never publish into the references folder.
                assert sorted(p.name for p in (work / "references").iterdir()) == [
                    "genes.bed",
                    "sc.fa.gz",
                ]
        label = f"Real-tool ({engine}) scATAC" if engine else "scATAC Nextflow regression"
        print(f"{label} validation passed.", flush=True)
        passed = True
    finally:
        if args.keep or not passed:
            print(f"Artifacts retained: {work}", flush=True)
        else:
            shutil.rmtree(work)


if __name__ == "__main__":
    main()
