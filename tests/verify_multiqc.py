"""Deterministic identity/coverage checks and opt-in real MultiQC aggregation cases."""

from __future__ import annotations

import argparse
import base64
import copy
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from collections.abc import Callable
from pathlib import Path

from genesis_tools.multiqc_reporting import canonicalize_html, inventory, verify_parsed

ROOT = Path(__file__).resolve().parents[1]


def rejected(function: Callable[[], object], diagnostic: str) -> None:
    try:
        function()
    except ValueError as error:
        assert diagnostic in str(error), str(error)
    else:
        raise AssertionError(f"Expected rejection: {diagnostic}")


def unit_checks(work: Path) -> None:
    """Test overwrite and missing-parser safeguards without MultiQC or Docker."""
    inputs = work / "inputs"
    inputs.mkdir(parents=True)
    sample = dict(
        id="a", species="Plant", ref_id="ref", reference_fasta="ref.fa.gz", layout="SE", control="-"
    )
    manifest = {"samples": [sample]}
    (inputs / "a.metadata.tsv").write_text(
        "sample_id\tspecies\treference_fasta\tlayout\na\tPlant\tref.fa.gz\tSE\n"
    )
    for suffix in ("main.stats.txt", "read1.stats.txt", "flagstat.txt", "idxstats.tsv"):
        (inputs / f"a.qc.report.{suffix}").write_text("unit fixture\n")
    with zipfile.ZipFile(inputs / "a.read1_fastqc.zip", "w") as archive:
        archive.writestr(
            "a.read1_fastqc/fastqc_data.txt", "##FastQC\t0.12.1\nFilename\ta.read1.fastq.gz\n"
        )
    expected, metadata, _sources = inventory(manifest, inputs)
    assert list(metadata) == ["a"] and len(expected["stats"]) == 2
    rejected(lambda: inventory({"samples": [sample, sample]}, inputs), "duplicate sample ID")
    duplicate = inputs / "duplicate"
    duplicate.mkdir()
    shutil.copyfile(inputs / "a.qc.report.flagstat.txt", duplicate / "a.qc.report.flagstat.txt")
    rejected(lambda: inventory(manifest, inputs), "Duplicate input basename")
    shutil.rmtree(duplicate)
    (inputs / "a.qc.report.flagstat.txt").unlink()
    rejected(lambda: inventory(manifest, inputs), "Missing or empty required flagstat")
    data = work / "data"
    data.mkdir()
    (data / "multiqc_data.json").write_text('{"report_saved_raw_data": {}}\n')
    rejected(lambda: verify_parsed(data, expected), "Parsed fastqc sample mismatch")
    payload = b'{"known_numeric_metric": 123}'
    rendered = []
    for index in range(2):
        html_path = work / f"report-{index}.html"
        encoded = base64.b64encode(gzip.compress(payload, mtime=index + 1)).decode()
        html_path.write_text(
            '<script type="text/plain" id="mqc_compressed_plotdata">'
            + encoded
            + "</script>\n"
            + f'reportUuid = "different-{index}";\nconfigCreationDate = "different-{index}";\n'
        )
        canonicalize_html(html_path, {"run": "unit"})
        rendered.append(html_path.read_bytes())
        match = re.search(r'id="mqc_compressed_plotdata">([^<]+)', html_path.read_text())
        assert match is not None
        assert gzip.decompress(base64.b64decode(match[1])) == payload
    assert rendered[0] == rendered[1]
    print("MultiQC input identity, missing-parser and lossless HTML determinism checks passed.")


def make_case(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True)
    for path in source.rglob("*"):
        if path.is_file() and (
            path.name.endswith("_fastqc.zip")
            or ".qc.report." in path.name
            or path.name.endswith(".metadata.tsv")
        ):
            target = destination / path.name
            assert not target.exists(), target
            shutil.copyfile(path, target)


def clone_sample(inputs: Path, original: str, replacement: str) -> None:
    """Change identities only in derived test fixtures, retaining numerical QC values."""
    for source in list(inputs.glob(f"{original}.*")):
        target = source.with_name(replacement + source.name[len(original) :])
        if source.suffix == ".zip":
            with zipfile.ZipFile(source) as before, zipfile.ZipFile(target, "w") as after:
                for item in before.infolist():
                    data = before.read(item.filename)
                    if item.filename.endswith(("fastqc_data.txt", "summary.txt")):
                        data = data.replace(original.encode(), replacement.encode())
                    after.writestr(item.filename.replace(original, replacement), data)
        elif source.suffix in (".txt", ".tsv"):
            target.write_text(source.read_text().replace(original, replacement))
        else:
            shutil.copyfile(source, target)


def docker_cases(source: Path, manifest_path: Path, work: Path) -> None:
    """Run production aggregation with optional omissions and hostile sample identities."""
    image = next(
        line.split("=", 1)[1]
        for line in (ROOT / "conf/genesis-images.lock").read_text().splitlines()
        if line.startswith("GENESIS_MULTIQC_IMAGE=")
    )
    baseline = json.loads(manifest_path.read_text())
    observations = {}
    for name in (
        "complete",
        "optional-missing",
        "collisions",
        "cleaning-collision",
        "duplicate",
        "missing-required",
        "missing-module",
        "complete-repeat",
    ):
        case = work / name
        inputs = case / "inputs"
        make_case(source, inputs)
        manifest = copy.deepcopy(baseline)
        config = (ROOT / "multiqc_config.yaml").read_text()
        if name == "optional-missing":
            spp = next(inputs.glob("*.spp.tsv"))
            spp.unlink()
        elif name in ("collisions", "cleaning-collision"):
            original = next(s for s in manifest["samples"] if s["layout"] == "SE")
            replacements = ("collision", "collision.bam", "collision.report")
            for replacement in replacements:
                clone_sample(inputs, original["id"], replacement)
            manifest["samples"] = [
                {**original, "id": replacement, "control": "-"} for replacement in replacements
            ]
            for path in list(inputs.iterdir()):
                if not path.name.startswith("collision."):
                    path.unlink()
            if name == "cleaning-collision":
                config = config.replace(
                    "fn_clean_sample_names: false", "fn_clean_sample_names: true"
                )
        elif name == "duplicate":
            manifest["samples"].append(manifest["samples"][0])
        elif name == "missing-required":
            next(inputs.glob("*.flagstat.txt")).unlink()
        elif name == "missing-module":
            config = config.replace("  - fastqc\n", "")
        (case / "manifest.json").write_text(json.dumps(manifest, sort_keys=True) + "\n")
        (case / "config.yaml").write_text(config)
        before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs.iterdir()}
        command = [
            "docker",
            "run",
            "--rm",
            "--network",
            "none",
            "-u",
            f"{os.getuid()}:{os.getgid()}",
            "-v",
            f"{ROOT}:/repo:ro",
            "-v",
            f"{case}:/work",
            "-w",
            "/work",
            image,
            "python",
            "/repo/genesis_tools/src/genesis_tools/multiqc_reporting.py",
            "--manifest",
            "manifest.json",
            "--inputs",
            "inputs",
            "--config",
            "config.yaml",
            "--output",
            ".",
        ]
        result = subprocess.run(command, capture_output=True, text=True, timeout=180, check=False)
        (case / "stdout.log").write_text(result.stdout)
        (case / "stderr.log").write_text(result.stderr)
        observations[name] = {"command": command, "exit": result.returncode, "input_sha256": before}
        (work / "observations.json").write_text(json.dumps(observations, indent=2) + "\n")
        assert before == {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs.iterdir()
        }
        if name in ("complete", "optional-missing", "collisions", "complete-repeat"):
            assert result.returncode == 0, result.stderr
            assert (case / "multiqc_report.html").stat().st_size > 0
            assert len(list(case.glob("*.html"))) == 1
            provenance = json.loads((case / "multiqc_data/genesis_provenance.json").read_text())
            assert provenance["run"] == baseline["run"]
            if name == "collisions":
                stats = json.loads((case / "multiqc_data/genesis_stats.json").read_text())
                assert len(stats) == 6
                for sample in ("collision", "collision.bam", "collision.report"):
                    assert sample + ".qc.report.main.stats.txt" in stats
                    assert sample + ".qc.report.read1.stats.txt" in stats
        else:
            assert result.returncode != 0, name
            diagnostic = {
                "duplicate": "duplicate sample ID",
                "missing-required": "required flagstat",
                "missing-module": "Parsed fastqc sample mismatch",
                "cleaning-collision": "Parsed fastqc sample mismatch",
            }[name]
            assert diagnostic in result.stderr, result.stderr
            if name in ("missing-module", "cleaning-collision"):
                assert "MultiQC complete" in result.stderr
    for kind in ("fastqc", "stats", "flagstat", "idxstats"):
        data = f"multiqc_data/genesis_{kind}.json"
        assert (work / "complete" / data).read_bytes() == (
            work / "complete-repeat" / data
        ).read_bytes()
        assert (work / "complete" / data).read_bytes() == (
            work / "optional-missing" / data
        ).read_bytes()
    assert (work / "complete/multiqc_report.html").read_bytes() == (
        work / "complete-repeat/multiqc_report.html"
    ).read_bytes()
    observations["fresh_report_byte_identity"] = {
        "html": (work / "complete/multiqc_report.html").read_bytes()
        == (work / "complete-repeat/multiqc_report.html").read_bytes(),
        "parsed_metrics": True,
    }
    (work / "observations.json").write_text(json.dumps(observations, indent=2) + "\n")
    print("MultiQC optional inputs, collisions, missing modules and deterministic metrics passed.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    if args.source is not None:
        assert args.manifest is not None and args.work is not None
        docker_cases(args.source.resolve(), args.manifest.resolve(), args.work.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix="genesis-multiqc-") as folder:
            unit_checks(Path(folder))


if __name__ == "__main__":
    main()
