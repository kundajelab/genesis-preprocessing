"""Real-tool two-library ATAC regression with an unchanged second resume."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch

from atac_fixture import create
from genesis_tools.atac.inputs import FIELDS, digest
from genesis_tools.atac.run import execute, prepare

ROOT = Path(__file__).resolve().parents[1]


def verify_adapters(folder: Path) -> None:
    """Check paired adapter removal, raw-read preservation and report publication."""
    folder.mkdir()
    raw = {}
    for mate, base in ((1, "A"), (2, "C")):
        sequence = base * 40 + "AGATCGGAAGAGC"
        path = folder / f"raw_R{mate}.fastq.gz"
        path.write_bytes(
            gzip.compress(f"@read/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n".encode(), mtime=0)
        )
        raw[path] = digest(path)
    image = json.loads((ROOT / "conf/atac-images.json").read_text())["cutadapt"]
    config = folder / "test.config"
    config.write_text(
        f"params.images = [cutadapt: '{image}']\n"
        f"params.outdir = '{folder}/output'\n"
        "process.executor = 'local'\nprocess.cpus = 1\nprocess.memory = '1 GB'\n"
        "docker.enabled = true\ndocker.runOptions = '-u $(id -u):$(id -g)'\n"
        f"trace.enabled = true\ntrace.file = '{folder}/trace.tsv'\n"
        "trace.fields = 'name,status,exit,realtime,peak_rss,workdir'\n"
    )
    workflow = folder / "test.nf"
    workflow.write_text(
        "nextflow.enable.dsl = 2\n"
        f"include {{ ATAC_ADAPTERS }} from '{ROOT}/modules/atac/processes'\n"
        "workflow {\n  ATAC_ADAPTERS(Channel.of(tuple("
        "[library_id:'adapter-test', adapter_r1:'AGATCGGAAGAGC', adapter_r2:'AGATCGGAAGAGC'],"
        f"file('{folder}/raw_R1.fastq.gz'),file('{folder}/raw_R2.fastq.gz'))))\n}}\n"
    )
    subprocess.run(
        ["nextflow", "-C", str(config), "run", str(workflow), "-w", str(folder / "work")],
        cwd=folder,
        check=True,
    )
    trace = list(csv.DictReader((folder / "trace.tsv").open(), delimiter="\t"))
    assert len(trace) == 1 and trace[0]["status"] == "COMPLETED"
    work = Path(trace[0]["workdir"])
    for mate, base in ((1, "A"), (2, "C")):
        with gzip.open(work / f"analysis_R{mate}.fastq.gz", "rt") as stream:
            assert stream.read().splitlines() == [f"@read/{mate}", base * 40, "+", "I" * 40]
    report = json.loads((folder / "output/adapter-test/qc/adapters.json").read_text())
    assert report["read_counts"]["input"] == report["read_counts"]["output"] == 1
    assert raw == {path: digest(path) for path in raw}


def create_study_fixture(folder: Path) -> tuple[Path, Path, Path]:
    """The same real, two-library/multi-lane fixture for standalone and study validation."""
    fixture = folder / "fixture"
    create(fixture, assay="bulk-ATAC")
    ref = {
        "reference_id": "toy",
        "species": "synthetic",
        "assembly": "synthetic-v1",
        "annotation_release": "synthetic-v1",
        "genome_size": 2002000,
        "genome_size_method": "synthetic total bases",
        "fasta": str(fixture / "reference.fa.gz"),
        "fasta_sha256": digest(fixture / "reference.fa.gz"),
        "tss": str(fixture / "tss.bed"),
        "tss_sha256": digest(fixture / "tss.bed"),
        "mitochondrial_contigs": ["mitochondria"],
        "plastid_contigs": ["plastid"],
    }
    registry = folder / "references.json"
    registry.write_text(json.dumps({"schema_version": 1, "references": [ref]}))
    # Divide a library into technical lanes; both use the same basenames in different folders.
    for mate in (1, 2):
        with gzip.open(fixture / f"pe_treatment.R{mate}.fastq.gz", "rt") as source:
            for lane in (1, 2):
                target = fixture / f"lane{lane}"
                target.mkdir(exist_ok=True)
                with gzip.open(target / f"R{mate}.fastq.gz", "wt") as output:
                    for _ in range(16000):
                        for _ in range(4):
                            output.write(source.readline())
    sheet = folder / "samples.tsv"
    rows = []
    for library, rep in [("first", "1"), ("second", "2")]:
        for lane in (1, 2):
            r1, r2 = [fixture / f"lane{lane}/R{mate}.fastq.gz" for mate in (1, 2)]
            rows.append(
                {
                    "library_id": library,
                    "sample_id": "plant",
                    "biological_replicate": rep,
                    "lane_id": str(lane),
                    "reference_id": "toy",
                    "read1": str(r1),
                    "read2": str(r2),
                    "read1_sha256": digest(r1),
                    "read2_sha256": digest(r2),
                    "mapq": 30,
                    "duplicates": "exclude",
                    "adapter_r1": "-",
                    "adapter_r2": "-",
                }
            )
    with sheet.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=sorted(FIELDS), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return sheet, registry, fixture


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args()
    base = ROOT / "tests/.runs"
    base.mkdir(exist_ok=True)
    folder = Path(tempfile.mkdtemp(prefix="supported-atac-", dir=base))
    print(f"ATAC regression evidence: {folder}", flush=True)
    sheet, registry, fixture = create_study_fixture(folder)
    before = {str(p): digest(p) for p in fixture.rglob("*") if p.is_file()}
    run = folder / "run"
    checkout = folder / "checkout"
    for name in ("conf", "workflows", "modules/atac", "genesis_tools/src"):
        shutil.copytree(ROOT / name, checkout / name, ignore=shutil.ignore_patterns("__pycache__"))

    def invoke(*, resume: bool = False) -> int:
        with patch("genesis_tools.atac.run.ROOT", checkout):
            return execute(sheet, registry, run, resume=resume)

    assert invoke() == 0
    output = run / "output"
    trace = list(csv.DictReader((output / "trace.tsv").open(), delimiter="\t"))
    assert len(trace) == 21 and all(r["status"] == "COMPLETED" for r in trace), trace
    shutil.copy2(output / "trace.tsv", folder / "fresh-trace.tsv")
    scientific = {
        str(p): digest(p)
        for p in output.rglob("*")
        if p.is_file() and p.suffix in {".gz", ".bw", ".narrowPeak"}
    }
    for library in ("first", "second"):
        reports = output / library / "qc"
        assert json.loads((reports / "adapters.json").read_text()) == {"policy": "none"}
        assert "macs3 3.0.5" in (output / library / "peaks/macs3.version.txt").read_text()
        assert len(list((reports / "fastqc").glob("*.zip"))) == 2
        modules = json.loads((reports / "multiqc_data/required-modules.json").read_text())
        assert set(modules.values()) == {2} and len(modules) == 4
        assert (output / library / "peaks/atac_peaks.narrowPeak").stat().st_size > 0
        acquisition = json.loads((output / library / "provenance/acquisition.json").read_text())
        assert acquisition["templates"] == 32000 and len(acquisition["sources"]) == 4
    assert not list(output.rglob("*.bam")) and not list(output.rglob("*.fastq.gz"))
    comparisons = json.loads((output / "replicate-comparisons.json").read_text())["comparisons"]
    assert len(comparisons) == 1 and comparisons[0]["peak_bp_jaccard"] == 1
    assert not comparisons[0]["same_biological_replicate"]
    for repetition in (1, 2):
        assert invoke(resume=True) == 0
        trace = list(csv.DictReader((output / "trace.tsv").open(), delimiter="\t"))
        assert len(trace) == 21 and all(r["status"] == "CACHED" for r in trace), trace
        shutil.copy2(output / "trace.tsv", folder / f"resume-{repetition}-trace.tsv")
        assert scientific == {p: digest(Path(p)) for p in scientific}
    # Only the isolated checkout is changed. Production helpers and inputs stay untouched.
    previous_identity = json.loads((run / "provenance.json").read_text())["source_bundle_sha256"]
    helper = checkout / "genesis_tools/src/genesis_tools/atac/stages.py"
    old_source = helper.read_bytes()
    helper.write_bytes(old_source + b"\n# Isolated helper revision for cache invalidation.\n")
    assert invoke(resume=True) == 0
    changed_identity = json.loads((run / "provenance.json").read_text())["source_bundle_sha256"]
    assert previous_identity != changed_identity
    old_bundle = run / "inputs" / f"code-{previous_identity}"
    assert (old_bundle / "genesis_tools/atac/stages.py").read_bytes() == old_source
    trace = list(csv.DictReader((output / "trace.tsv").open(), delimiter="\t"))
    helpers = {"ATAC_ACQUIRE", "ATAC_FRAGMENTS", "ATAC_ENRICHMENT", "ATAC_REPORT"}
    assert all(r["status"] == "COMPLETED" for r in trace if r["name"].split()[0] in helpers)
    assert sum(r["name"].split()[0] in helpers for r in trace) == 8
    shutil.copy2(output / "trace.tsv", folder / "changed-helper-trace.tsv")
    assert scientific == {p: digest(Path(p)) for p in scientific}
    assert invoke(resume=True) == 0
    trace = list(csv.DictReader((output / "trace.tsv").open(), delimiter="\t"))
    assert len(trace) == 21 and all(r["status"] == "CACHED" for r in trace)
    shutil.copy2(output / "trace.tsv", folder / "changed-helper-stable-trace.tsv")
    assert before == {p: digest(Path(p)) for p in before}
    original = (run / "inputs/first.json").read_bytes()
    sheet.write_text(sheet.read_text().replace("\t30\t", "\t10\t"))
    try:
        prepare(sheet, registry, run)
    except ValueError:
        pass
    else:
        raise AssertionError("Changed run inputs accepted")
    assert (run / "inputs/first.json").read_bytes() == original
    verify_adapters(folder / "adapter-test")
    print(
        "PASS: two libraries, technical lanes, 21 cached tasks twice; "
        "raw/scientific outputs unchanged; changed helpers invalidate tasks"
    )
    if not args.keep:
        shutil.rmtree(folder)


if __name__ == "__main__":
    main()
