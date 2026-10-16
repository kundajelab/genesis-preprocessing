"""Create a new source catalog after deliberate evidence review; never called by build."""

from __future__ import annotations

import json
import subprocess

from evidence import COHORTS, HERE, MACS, NEW, OLD, ROOT, digest


def main() -> None:
    names = []
    for cohort in COHORTS:
        names.extend(
            f"{OLD}/{cohort}/{name}"
            for name in (
                "manifest.json",
                "evidence.tar.gz",
                "metrics.tsv",
                "processing.tsv",
                "resources.tsv",
                "library-qc.tsv",
            )
        )
    names.extend(
        f"{NEW}/{name}"
        for name in (
            "figures/quality-summary.tsv",
            "figures/cost-summary.tsv",
            "figures/cost-observations.tsv",
            "figures/index-observations.tsv",
            "dap-se-5m-motifs/motif-metrics.tsv",
            "annotation-sensitivity-1m.json",
            "reference-concordance-1m.json",
            "caller-diagnostics-colcen-5m.json",
        )
    )
    names.extend(
        f"{MACS}/{name}"
        for name in (
            "timings.tsv",
            "results.json",
            "input-context.json",
            "builds.json",
            "verification.json",
        )
    )
    names.extend(
        [
            "benchmarks/reports/local-audit.json",
            "validation/reports/expanded-plant-comparison.md",
            "validation/reports/pipeline-method-comparison.md",
            "validation/reports/macs3-optimization.md",
            "benchmarks/experiments/plant-methods-20261013/README.md",
        ]
    )
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    sources = []
    for name in names:
        path = ROOT / name
        history = subprocess.check_output(
            ["git", "log", "-1", "--format=%H", "--", name], cwd=ROOT, text=True
        ).strip()
        sources.append(
            {
                "path": name,
                "sha256": digest(path),
                "bytes": path.stat().st_size,
                "source_revision": history
                or "local extraction; original hashes in local-audit.json",
                "evidence_kind": "reported historical narrative"
                if name.endswith(".md")
                else "structured recorded evidence",
            }
        )
    with (HERE / "sources.json").open("x") as handle:
        json.dump(
            {
                "schema_version": 1,
                "snapshot_date": "2026-10-16",
                "reviewed_revision": revision,
                "sources": sources,
                "scope": "Existing evidence rechecked; no biological pipeline rerun",
            },
            handle,
            indent=2,
        )
        handle.write("\n")


if __name__ == "__main__":
    main()
