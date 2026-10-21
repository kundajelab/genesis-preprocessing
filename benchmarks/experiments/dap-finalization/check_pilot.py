"""Bounded offline checks for pilot evidence and protected delivery time."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from pilot import admit, digest, forecast, inspect, verify_binding


def main() -> None:
    observed = [{"templates": 50000, "wall_seconds": 100, "setup_seconds": 20}]
    value = forecast(observed, 1000000)
    assert value["conservative_seconds"] == 3240 and value["second_pilot_required"]
    value = forecast(
        observed + [{"templates": 200000, "wall_seconds": 500, "setup_seconds": 20}], 1000000
    )
    assert value["conservative_seconds"] == 4840 and not value["second_pilot_required"]
    for bad in (float("nan"), float("inf"), -1, True):
        try:
            forecast([{"templates": 50000, "wall_seconds": bad, "setup_seconds": 0}], 1000000)
        except TypeError, ValueError:
            pass
        else:
            raise AssertionError("Invalid measurement accepted")
    assert admit(now=100, delivery=25000, estimate=1000, remaining_required=2000)
    assert not admit(now=100, delivery=24000, estimate=1000, remaining_required=2000)
    assert not admit(now=100, delivery=25000, estimate=1000, remaining_required=3000)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        assert inspect(root)["status"] == "NEEDS_CORRECTION"
        artifact = root / "zero-peaks.bed"
        artifact.write_bytes(b"")
        result = {
            "status": "COMPLETE",
            "inputs_unchanged": True,
            "wall_seconds": 1,
            "results": [
                {
                    "arm": "baseline",
                    "status": "COMPLETE",
                    "libraries": {
                        "sample": {
                            "artifacts": {str(artifact): digest(artifact)},
                            "metrics": {"independent_count_mismatches": {"value": 0}},
                        }
                    },
                }
            ],
        }
        p = root / "results.json"
        p.write_text(json.dumps(result))
        before = p.read_bytes()
        assert inspect(root)["status"] == "READY", "Zero peaks are not structural failure"
        assert p.read_bytes() == before
        artifact.write_text("changed")
        assert inspect(root)["status"] == "NEEDS_CORRECTION"
        binding = {}
        for name in ("target_spec", "target_dataset", "pilot_spec", "pilot_dataset"):
            binding[name] = str(p)
            binding[name + "_sha256"] = digest(p)
        verify_binding(binding)
        p.write_text("{}")
        try:
            verify_binding(binding)
        except ValueError:
            pass
        else:
            raise AssertionError("Stale pilot admitted")
    print(
        "PASS pilot forecasts, reserves, invalid measurements, stale evidence, "
        "zero peaks and no mutation"
    )


if __name__ == "__main__":
    main()
