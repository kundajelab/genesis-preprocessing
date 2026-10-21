"""Local comparison recipe with verified shared inputs, not fresh timing repeats.

Reuses identical reference indexes and alignments across single-repetition cases.
The receipt identifies this wrapper separately from the core benchmark source.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from genesis_tools.benchmark import methods
from genesis_tools.benchmark.runner import run
from genesis_tools.benchmark.spec import load, sha256, write_json


def identity(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def execute(spec: Path, out: Path, cache_path: Path) -> None:
    plan = load(spec, verify=True)
    if plan["status"] != "READY":
        raise ValueError(plan["blockers"])
    if plan["spec"]["performance"]["repetitions"] != 1:
        raise ValueError("Reuse recipe cannot measure fresh timing repetitions")
    if out.exists():
        raise ValueError("Use a new output path")
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    events: list[dict[str, Any]] = []
    original = methods.Methods
    wrapper_hash = sha256(Path(__file__))
    common = {
        "recipe": plan["source"]["recipe_sha256"],
        "wrapper": wrapper_hash,
        "environment": plan["spec"]["environment"],
        "assay": plan["spec"]["assay"],
        "cpus": plan["spec"]["limits"]["cpus"],
        "reference": plan["dataset"]["reference"],
    }

    def retrieve(key: str) -> Path | None:
        row = cache.get(key)
        if row and all(
            Path(p).is_file() and sha256(Path(p)) == h for p, h in row["outputs"].items()
        ):
            events.append({"status": "REUSED", "key": key, **row})
            return Path(row["path"])
        return None

    def publish(key: str, path: Path, outputs: dict[str, str]) -> None:
        cache[key] = {"path": str(path), "outputs": outputs}
        write_json(cache_path, cache)
        events.append({"status": "GENERATED", "key": key, **cache[key]})

    class ReusingMethods(original):
        def reference(self, aligner: str) -> Path:
            key = identity({"kind": "index", "common": common, "aligner": aligner})
            cached = retrieve(key)
            if cached:
                return cached
            fasta = super().reference(aligner)
            receipt = json.loads((fasta.parent / "task.json").read_text())
            publish(key, fasta, receipt["outputs"])
            return fasta

        def align(
            self, folder: Path, lib: dict[str, Any], method: dict[str, Any], fasta: Path
        ) -> Path:
            policy = {
                k: method.get(k)
                for k in ("aligner", "bowtie2_mode", "trimmer", "adapter_r1", "adapter_r2")
            }
            reads = [plan["assets"][lib[k]]["sha256"] for k in ("read1", "read2") if k in lib]
            key = identity(
                {
                    "kind": "alignment",
                    "common": common,
                    "library": lib["id"],
                    "layout": lib["layout"],
                    "reads": reads,
                    "policy": policy,
                }
            )
            cached = retrieve(key)
            if cached:
                target = folder / lib["id"] / "alignment" / "aligned.bam"
                target.parent.mkdir(parents=True, exist_ok=True)
                target.symlink_to(cached)
                return target
            bam = super().align(folder, lib, method, fasta)
            publish(key, bam, {str(bam): sha256(bam)})
            return bam

    methods.Methods = ReusingMethods
    try:
        run(plan, out)
        for event in events:
            if any(
                not Path(p).is_file() or sha256(Path(p)) != h for p, h in event["outputs"].items()
            ):
                raise ValueError("Shared artifact changed during experiment")
    finally:
        methods.Methods = original
        out.mkdir(exist_ok=True)
        write_json(
            out / "reuse.json",
            {
                "wrapper_sha256": wrapper_hash,
                "core_recipe_sha256": plan["source"]["recipe_sha256"],
                "scope": (
                    "Scientific comparison with verified reuse; "
                    "elapsed time is not cold workflow timing"
                ),
                "events": events,
            },
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    execute(args.experiment, args.out, args.cache)
