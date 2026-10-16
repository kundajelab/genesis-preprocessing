"""Freeze a compact, path-redacted inventory of completed local reference evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence import digest


def audit(root: Path) -> dict:
    records = []
    seen = set()
    for path in sorted(root.glob("*frozen-results/results.json")):
        experiment = path.with_name("experiment.json")
        identity = digest(experiment)
        if identity in seen:
            continue
        seen.add(identity)
        data = json.loads(path.read_text())
        records.append(
            {
                "name": path.parent.name,
                "status": data["status"],
                "experiment_sha256": identity,
                "results_sha256": digest(path),
                "scope": "Completed record inspected; scientific outputs not independently rerun",
            }
        )
    path = (
        root / "expanded-panel/assessments/atac-arabidopsis-5m/result/comparison.json"
    )
    data = json.loads(path.read_text())
    status_path = root / "expanded-panel/assessments/status.json"
    status = json.loads(status_path.read_text())
    if (
        status["jobs"]["atac-arabidopsis-5m"]["result_sha256"] != digest(path)
        or data["status"] != "COMPLETE"
    ):
        raise ValueError("Assessment receipt mismatch")
    libraries = {}
    for library, row in data["libraries"].items():
        libraries[library] = {"metrics": row["metrics"]}
        for side in ("left_mapping", "right_mapping"):
            libraries[library][side] = {
                key: row[side][key]
                for key in ("counts", "input_peaks", "unique_reciprocal_fraction")
            }
    return {
        "snapshot": "2026-10-16",
        "frozen_records": records,
        "assessment_statuses": {k: v["status"] for k, v in status["jobs"].items()},
        "reference_5m": {
            "source_sha256": digest(path),
            "status_receipt_sha256": digest(status_path),
            "status": data["status"],
            "experiments": data["experiments"],
            "scope": data["scope"],
            "mapping_sha256": data["mapping_sha256"],
            "script_sha256": data["script_sha256"],
            "libraries": libraries,
        },
        "limitation": "Compact extraction from local records; underlying local files are not distributed here. Not a new experiment or independently rerun analysis.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with args.out.open("x") as handle:
        json.dump(audit(args.root), handle, indent=2, sort_keys=True)
        handle.write("\n")
