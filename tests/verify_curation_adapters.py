"""Legacy DAP/ATAC import boundaries and controller-local provenance observations."""

from __future__ import annotations

import contextlib
import hashlib
import json
import sqlite3
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from genesis_tools.contracts.adapters import DAP_UNRECORDED, dap_run_record, from_atac, from_dap
from genesis_tools.contracts.records import dump, fingerprint, record
from genesis_tools.contracts.validation import bundle
from genesis_tools.curation.status_page import summary
from genesis_tools.execution.controller import initialize
from genesis_tools.registry import store
from verify_curation import altered, fixture, rejects


def legacy(root: Path) -> None:
    run = root / "atac"
    run.mkdir()
    library = {
        "library_id": "lib",
        "sample_id": "sample",
        "biological_replicate": "rep1",
        "reference": {"reference_id": "ref", "species": "plant"},
        "lanes": [{"lane_id": "one"}, {"lane_id": "two"}],
        "mapq": 30,
        "duplicates": "exclude",
        "adapter_r1": "-",
        "adapter_r2": "-",
    }
    dump(run / "manifest.json", {"schema_version": 1, "libraries": [library]})
    digest = hashlib.sha256((run / "manifest.json").read_bytes()).hexdigest()
    dump(
        run / "provenance.json",
        {
            "manifest_sha256": digest,
            "containers": {"synthetic": "test-only"},
            "source_sha256": {"fixture": "0" * 64},
            "git_sha": "0" * 40,
            "exit_code": 0,
        },
    )
    dump(
        run / "output/lib/fragments/metrics.json",
        {
            "counts": {
                "total_templates": 20,
                "eligible_before_dedup": 10,
                "usable_fragments": 9,
                "secondary_or_supplementary_records": 2,
            },
            "NRF": 0.9,
            "duplicate_fraction": 0.1,
            "PBC1": None,
            "PBC2": None,
            "metric_version": 1,
        },
    )
    dump(
        run / "output/lib/qc/enrichment.json",
        {
            "FRiP": 0,
            "peak_count": 0,
            "usable_fragments": 9,
            "fragments_in_peaks": 0,
            "metric_version": 1,
        },
    )
    before = (run / "manifest.json").read_bytes()
    manifests = from_atac(run, "accession-study", "remote-a", "study")
    data = manifests[0]["data"]
    assert len(data["lanes"]) == 2 and data["run"]["execution_state"] == "SUCCEEDED"
    measured = {m["data"]["metric"]: m["data"] for m in data["measurements"]}
    assert measured["FRiP"]["value"] == 0 and measured["FRiP"]["denominator_value"] == 9
    assert measured["secondary_or_supplementary_records"]["units"] == "records"
    assert measured["usable_fragments"]["units"] == "fragments"
    assert any(a["role"] == "fastqc" and a["required"] for a in data["artifacts"])
    assert (run / "manifest.json").read_bytes() == before
    dump(run / "provenance.json", {"manifest_sha256": "f" * 64})
    rejects(lambda: from_atac(run, "study", "local", None), "does not match")

    dap = root / "dap"
    provenance = {
        "run": {"git_sha": "0" * 40},
        "samples": [
            {
                "id": name,
                "species": "plant",
                "reference_fasta": "plant.fa.gz",
                "ref_id": "plant_ref",
                "layout": layout,
                "control": control,
            }
            for name, layout, control in (
                ("control", "PE", "-"),
                ("TF1", "PE", "control"),
                ("TF2", "SE", "control"),
            )
        ],
    }
    dump(dap / "output/multiqc/multiqc_data/genesis_provenance.json", provenance)
    dump(
        dap / "output/multiqc/multiqc_data/genesis_flagstat.json",
        {
            "TF1.qc.report.flagstat.txt": {"flagstat_total": 100},
        },
    )
    manifests = from_dap(dap, "dap-study", "local", "study", None, None)
    assert len({m["id"] for m in manifests}) == 3
    control, tf1, tf2 = (m["data"] for m in manifests)
    assert not any(a["role"] == "peaks" for a in control["artifacts"])
    assert tf1["metadata"]["control"] == tf2["metadata"]["control"] == "control"
    assert sum(a["role"] == "fastqc" for a in tf1["artifacts"]) == 2
    assert sum(a["role"] == "fastqc" for a in tf2["artifacts"]) == 1
    assert tf1["measurements"][0]["data"]["value"] == 100
    # No git-status.txt or image-lock copy: incomplete, with the reasons recorded.
    assert all(m["data"]["run"]["provenance_complete"] is False for m in manifests)
    assert "DAP run has no git-status.txt" in tf1["metadata"]["provenance_limitations"]
    # A clean checkout run in pinned containers keeps its evidence, but stays incomplete:
    # code and image identity say nothing about the reads, reference or parameters used.
    output = dap / "output"
    commit = "a" * 40
    digest = "kundajelab/dap_seq_peaks@sha256:" + "b" * 64
    (output / "dap-images.lock").write_text(f"DAP_SEQ_PEAKS_IMAGE={digest}\n")
    (output / "git-status.txt").write_text(
        f"branch: main\ncommit: {commit}\ncontainer engine: docker\ndirty files: none\n"
    )
    clean = from_dap(dap, "dap-study", "local", "study", None, None)
    assert all(m["data"]["run"]["provenance_complete"] is False for m in clean)
    assert all(
        m["data"]["metadata"]["provenance_limitations"] == sorted(DAP_UNRECORDED) for m in clean
    )
    run = clean[0]["data"]["run"]
    assert run["git_sha"] == commit
    assert run["parameters"]["images"] == {"DAP_SEQ_PEAKS_IMAGE": digest}
    assert run["parameters"]["container_engine"] == "docker"
    dap_run_cases(output, commit, digest)


def dap_run_cases(output: Path, commit: str, digest: str) -> None:
    """Each departure from a clean, containerized, digest-pinned run is reported."""
    clean = f"branch: main\ncommit: {commit}\ncontainer engine: docker\ndirty files: none\n"
    cases = [
        (
            f"branch: main\ncommit: {commit}\ncontainer engine: docker\n"
            "dirty files:\n M main.nf\n?? notes.txt\n",
            f"DAP_SEQ_PEAKS_IMAGE={digest}\n",
            ["DAP run used a checkout with uncommitted or untracked files"],
        ),
        (
            clean.replace("container engine: docker", "container engine: none"),
            f"DAP_SEQ_PEAKS_IMAGE={digest}\n",
            ["DAP run did not use a container engine"],
        ),
        (
            clean,
            "DAP_SEQ_PEAKS_IMAGE=kundajelab/dap_seq_peaks:latest\n",
            ["DAP run image pins are not all immutable digests"],
        ),
        (clean, "", ["DAP run has no recorded image pins"]),
        (
            "git status unavailable: /x is not a Git checkout, or git is not installed\n"
            "container engine: docker\n",
            f"DAP_SEQ_PEAKS_IMAGE={digest}\n",
            ["DAP run lacks a recorded Git commit"],
        ),
    ]
    for status, lock, expected in cases:
        (output / "git-status.txt").write_text(status)
        (output / "dap-images.lock").write_text(lock)
        recorded, evidence, limitations = dap_run_record(output)
        assert limitations == expected, (status, lock, limitations)
        assert recorded == (None if "unavailable" in status else commit)
    (output / "git-status.txt").write_text(cases[0][0])
    assert dap_run_record(output)[1]["dirty_files"] == [" M main.nf", "?? notes.txt"]


def observations(root: Path) -> None:
    directory = root / "registry"
    store.initialize(directory)
    base = fixture(root / "source")

    # Equal local paths on two nodes remain separate physical locations.
    def remote(worker: str) -> dict[str, Any]:
        data = json.loads(json.dumps(base["data"]))
        data["worker"] = worker
        for item in data["artifacts"]:
            item["worker"] = worker
        return record("manifest", data, base["id"])

    first = bundle([remote("node-a")])
    second = bundle([remote("node-b")])
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: store.import_bundle(directory, first), range(2)))
    assert sorted(r["imported"] for r in results) == [0, 1]
    store.import_bundle(directory, second)
    with contextlib.closing(store.connect(directory)) as db:
        assert db.execute("SELECT COUNT(DISTINCT worker) FROM locations").fetchone()[0] == 2
        assert db.execute("SELECT COUNT(DISTINCT path) FROM locations").fetchone()[0] == 2
    pruned = altered(
        second["data"]["manifests"][0],
        lambda d: d["artifacts"][0].update(
            required=False, availability="intentionally_pruned", retention="intermediate"
        ),
    )
    store.import_bundle(directory, bundle([pruned]))
    with store.write(directory) as db:
        assert store.search(db, filters={"availability": "intentionally_pruned"})
        db.execute("UPDATE meta SET value='false' WHERE key='fts5'")
        assert store.search(db, text="library") and not store.search(db, text="%")

    plan = {
        "schema_version": 1,
        "workers": {
            "local": {
                "transport": "local",
                "slots": 1,
                "repo": str(root),
                "root": str(root / "worker"),
                "python": sys.executable,
            }
        },
        "jobs": [
            {
                "id": "job",
                "worker": "local",
                "payload": {
                    "git_sha": "0" * 40,
                    "argv": ["true"],
                    "inputs": [{"path": str(root / "input"), "sha256": "0" * 64}],
                },
            }
        ],
    }
    dump(root / "plan.json", plan)
    campaign = root / "campaign"
    initialize(campaign, root / "plan.json")  # Initialize only; never launch the job.
    assert store.import_campaign(directory, campaign)["imported"] == 1
    assert store.import_campaign(directory, campaign)["imported"] == 0

    def link(data: dict[str, Any]) -> None:
        data.update(name="<script>alert(1)</script>")
        data["run"].update(campaign_id=fingerprint(plan), job_id="job", attempt=1)

    linked = altered(base, link)
    store.import_bundle(directory, bundle([linked]))
    page = summary(directory, campaign)
    assert "&lt;script&gt;" in page and "<script>" not in page
    with sqlite3.connect(campaign / "state.sqlite") as db:
        assert db.execute("SELECT state,attempt FROM jobs").fetchone() == ("QUEUED", 0)


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="genesis-curation-adapters-") as temporary:
        root = Path(temporary)
        legacy(root)
        observations(root)
    print("PASS: legacy adapters, worker paths, concurrent imports, campaign evidence and HTML")


if __name__ == "__main__":
    main()
