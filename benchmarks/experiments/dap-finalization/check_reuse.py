"""Opt-in real-tool SE/PE fixture: reused alignments must preserve all measured metrics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from genesis_tools.benchmark.methods import fixture
from genesis_tools.benchmark.spec import write_json
from pilot import inspect
from run_native import execute


def main(out: Path, environment_checkout: Path) -> None:
    if out.exists():
        raise ValueError("Choose a new fixture directory")
    spec = fixture(out / "input", assay="DAP-seq", small=False, checkout=environment_checkout)
    value = json.loads(spec.read_text())
    baseline = value["methods"][0]
    value["methods"] = [baseline, dict(baseline, id="same-alignment")]
    value["limits"]["max_arms"] = 2
    write_json(spec, value)
    execute(spec, out / "first", out / "cache.json")
    execute(spec, out / "second", out / "cache.json")
    first = json.loads((out / "first/results.json").read_text())
    second = json.loads((out / "second/results.json").read_text())
    for left, right in zip(first["results"], second["results"], strict=True):
        for library in left["libraries"]:
            assert left["libraries"][library]["metrics"] == right["libraries"][library]["metrics"]
    reused = json.loads((out / "second/reuse.json").read_text())
    assert reused["events"] and all(row["status"] == "REUSED" for row in reused["events"])
    write_json(out / "inspection.json", inspect(out / "second"))
    print("PASS real-tool SE/PE reuse: all metrics unchanged; indexes and alignments reused")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--environment-checkout", type=Path, required=True)
    args = parser.parse_args()
    main(args.out.resolve(), args.environment_checkout.resolve())
