"""Read frozen benchmark evidence without executing scientific workflows."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import tarfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OLD = "validation/evidence/method-comparison-20261013"
NEW = "validation/evidence/method-comparison-20261014"
MACS = "validation/evidence/macs3-optimization-20261016"
COHORTS = [
    f"{assay}-{depth}m"
    for assay in ("dap-se", "dap-pe", "atac-arabidopsis")
    for depth in (1, 5)
]


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value: str) -> float | None:
    if value == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Non-finite reported value")
    return result


def verify_sources(
    root: Path = ROOT, catalog: Path = HERE / "sources.json"
) -> dict[str, Any]:
    data = json.loads(catalog.read_text())
    seen: set[str] = set()
    for record in data["sources"]:
        name = record["path"]
        if name in seen:
            raise ValueError(f"Duplicate source: {name}")
        seen.add(name)
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError("Source escapes checkout")
        if digest(path) != record["sha256"]:
            raise ValueError(f"Stale source: {name}")
    return data


def archive_json(path: Path, member: str) -> dict[str, Any]:
    # Never extract archive paths onto the filesystem.
    with tarfile.open(path) as archive:
        entries = [item for item in archive.getmembers() if item.name == member]
        if (
            len(entries) != 1
            or not entries[0].isfile()
            or entries[0].size > 100_000_000
        ):
            raise ValueError(f"Invalid evidence member: {member}")
        handle = archive.extractfile(entries[0])
        if handle is None:
            raise ValueError(member)
        return json.load(handle)


def unique_rows(rows: list[dict[str, str]], keys: tuple[str, ...]) -> None:
    seen: set[tuple[str, ...]] = set()
    for row in rows:
        identity = tuple(row[key] for key in keys)
        if identity in seen:
            raise ValueError(f"Duplicate observation: {identity}")
        seen.add(identity)


class Evidence:
    def __init__(self) -> None:
        self.catalog = verify_sources()
        self.quality = read_table(ROOT / NEW / "figures/quality-summary.tsv")
        unique_rows(self.quality, ("cohort", "method", "library", "metric"))
        definitions: dict[tuple[str, str], tuple[str, str, str]] = {}
        for row in self.quality:
            signature = (row["unit"], row["definition"], row["version"])
            key = (row["metric"], row["library"])
            if key in definitions and definitions[key] != signature:
                raise ValueError(f"Incompatible metric: {row['metric']}")
            definitions[key] = signature
            number(row["median"])
        self.cost = read_table(ROOT / NEW / "figures/cost-summary.tsv")
        self.macs = read_table(ROOT / MACS / "timings.tsv")
        self.motifs = read_table(ROOT / NEW / "dap-se-5m-motifs/motif-metrics.tsv")
        self.annotation = json.loads(
            (ROOT / NEW / "annotation-sensitivity-1m.json").read_text()
        )
        self.reference = json.loads(
            (ROOT / NEW / "reference-concordance-1m.json").read_text()
        )
        self.local = json.loads((HERE / "local-audit.json").read_text())
        self.samples: list[dict[str, Any]] = []
        self.cases: list[dict[str, Any]] = []
        self.counts: list[dict[str, str]] = []
        self.resources: list[dict[str, str]] = []
        for cohort in COHORTS:
            base = ROOT / OLD / cohort
            manifest = json.loads((base / "manifest.json").read_text())
            archive = base / "evidence.tar.gz"
            if digest(archive) != manifest["archive_sha256"]:
                raise ValueError("Archive does not match publication manifest")
            experiment = archive_json(archive, "experiment.json")
            result = archive_json(archive, "results.json")
            if result["status"] != "COMPLETE" or not result["inputs_unchanged"]:
                raise ValueError("Incomplete or changed benchmark inputs")
            for library in experiment["dataset"]["libraries"]:
                self.samples.append(
                    {
                        "cohort": cohort,
                        **library,
                        "assembly": experiment["dataset"]["reference"]["assembly"],
                        "experiment_id": experiment["experiment_id"],
                        "sampling": experiment["dataset"]["description"],
                    }
                )
            processing = read_table(base / "processing.tsv")
            unique_rows(processing, ("arm", "repetition"))
            for row in processing:
                if row["status"] != "COMPLETE" or row["resume"] != "False":
                    raise ValueError("Incomplete/resumed timing case")
                self.cases.append({"cohort": cohort, **row})
            for row in read_table(base / "metrics.tsv"):
                if row["metric"] == "independent_count_mismatches":
                    self.counts.append({"cohort": cohort, **row})
            for row in read_table(base / "resources.tsv"):
                self.resources.append({"cohort": cohort, **row})
        unique_rows(self.cases, ("cohort", "arm", "repetition"))

    def q(self, cohort: str, method: str, library: str, metric: str) -> float | None:
        rows = [
            r
            for r in self.quality
            if (r["cohort"], r["method"], r["library"], r["metric"])
            == (cohort, method, library, metric)
        ]
        if len(rows) != 1:
            raise ValueError(
                f"Missing/nonunique metric: {cohort}/{method}/{library}/{metric}"
            )
        return number(rows[0]["median"])

    def c(self, cohort: str, method: str) -> dict[str, str]:
        return next(
            r
            for r in self.cost
            if r["cohort"] == cohort
            and r["method"] == method
            and r["measurement"] == "processing_seconds"
        )


def tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Empty data table")
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, keys, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_ALL
        )
        writer.writeheader()
        writer.writerows(rows)
