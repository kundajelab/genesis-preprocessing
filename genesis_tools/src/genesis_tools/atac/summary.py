"""Summarize libraries without pooling biological replicates."""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path
from typing import Any

from .inputs import digest
from .metrics import intersection, intervals, merged


def summarize(manifest: dict[str, Any], output: Path) -> None:
    rows = []
    for library in manifest["libraries"]:
        folder = output / library["library_id"]
        metrics = json.loads((folder / "fragments/metrics.json").read_text())
        enrichment = json.loads((folder / "qc/enrichment.json").read_text())
        counts = metrics["counts"]
        rows.append(
            {
                key: library[key]
                for key in ("library_id", "sample_id", "biological_replicate", "reference_id")
            }
            | {
                "templates": counts["total_templates"],
                "mapped_pairs": counts.get("mapped_pairs", 0),
                "usable_fragments": counts.get("usable_fragments", 0),
                "mitochondrial_pairs": counts.get("raw_mitochondrial_pairs", 0),
                "plastid_pairs": counts.get("raw_plastid_pairs", 0),
                "duplicate_fraction": metrics["duplicate_fraction"],
                "peak_count": enrichment["peak_count"],
                "FRiP": enrichment["FRiP"],
            }
        )
    with (output / "libraries.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    comparisons = []
    for first, second in itertools.combinations(manifest["libraries"], 2):
        if (
            first["sample_id"] != second["sample_id"]
            or first["reference_id"] != second["reference_id"]
        ):
            continue
        a, _ = intervals(output / first["library_id"] / "peaks/atac_peaks.narrowPeak")
        b, _ = intervals(output / second["library_id"] / "peaks/atac_peaks.narrowPeak")
        shared = intersection(a, b)
        union = sum(e - s for _, s, e in merged(a)) + sum(e - s for _, s, e in merged(b)) - shared
        comparisons.append(
            {
                "first": first["library_id"],
                "second": second["library_id"],
                "same_biological_replicate": first["biological_replicate"]
                == second["biological_replicate"],
                "peak_bp_jaccard": shared / union if union else None,
            }
        )
    (output / "replicate-comparisons.json").write_text(
        json.dumps(
            {
                "comparisons": comparisons,
                "definition": "peak base-pair intersection / union",
                "pooled": False,
                "plant_thresholds": "UNSPECIFIED",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    files = {
        str(p.relative_to(output)): digest(p)
        for p in sorted(output.rglob("*"))
        if p.is_file() and p.suffix not in {".html", ".tsv"} and p.name != "checksums.json"
    }
    (output / "checksums.json").write_text(json.dumps(files, indent=2, sort_keys=True) + "\n")
