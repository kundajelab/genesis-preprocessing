"""Read-only pilot inspection and conservative admission for local DAP recipes."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def positive(value: object, name: str, *, zero: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be numeric")
    if not math.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"{name} is out of range")
    return float(value)


def forecast(observations: list[dict[str, Any]], target_templates: int) -> dict[str, Any]:
    """Use the slowest measured per-template rate; do not assume linear certainty."""
    target = positive(target_templates, "target_templates")
    if not observations:
        raise ValueError("A measured pilot is required")
    rates, overheads, depths = [], [], []
    for observation in observations:
        count = positive(observation["templates"], "templates")
        total = positive(observation["wall_seconds"], "wall_seconds")
        setup = positive(observation["setup_seconds"], "setup_seconds", zero=True)
        if setup >= total:
            raise ValueError("Setup time must be below total wall time")
        rates.append((total - setup) / count)
        overheads.append(setup)
        depths.append(count)
    seconds = 2 * (max(overheads) + target * max(rates))
    return {
        "conservative_seconds": math.ceil(seconds),
        "second_pilot_required": max(depths) < min(target, 200000) and seconds > 1800,
        "assumption": (
            "Twice measured setup plus slowest per-template processing/scoring rate; "
            "not a runtime guarantee"
        ),
    }


def admit(
    *,
    now: float,
    delivery: float,
    estimate: float,
    remaining_required: float,
    report_seconds: float = 10800,
    contingency_seconds: float = 10800,
) -> bool:
    for name, value in [("now", now), ("delivery", delivery), ("estimate", estimate)]:
        positive(value, name)
    for name, value in [
        ("remaining_required", remaining_required),
        ("report_seconds", report_seconds),
        ("contingency_seconds", contingency_seconds),
    ]:
        positive(value, name, zero=True)
    return now + estimate + remaining_required + report_seconds + contingency_seconds <= delivery


def inspect(root: Path, *, approved_warnings: tuple[str, ...] = ()) -> dict[str, Any]:
    """Inspect existing complete results; never create, repair or approve scientific data."""
    result_path = root / "results.json"
    if not result_path.is_file():
        return {"status": "NEEDS_CORRECTION", "reasons": ["No completed result receipt"]}
    result = json.loads(result_path.read_text())
    reasons: list[str] = []
    if result.get("status") != "COMPLETE" or result.get("inputs_unchanged") is not True:
        reasons.append("Incomplete run or changed inputs")
    records = result.get("results", [])
    if not records:
        reasons.append("No measured candidate records")
    metrics = []
    for record in records:
        if record.get("status") != "COMPLETE":
            reasons.append(f"Candidate incomplete: {record.get('arm')}")
        for library, value in record.get("libraries", {}).items():
            for name, expected in value.get("artifacts", {}).items():
                path = Path(name)
                if not path.is_file() or digest(path) != expected:
                    reasons.append(f"Artifact absent or changed for {library}")
            measured = value.get("metrics", {})
            mismatches = measured.get("independent_count_mismatches", {}).get("value")
            if mismatches != 0:
                reasons.append(f"Independent counting is missing or differs for {library}")
            metrics.append({"arm": record["arm"], "library": library, "metrics": measured})
    warnings = sorted(
        {
            line.strip()
            for path in (root / "commands").glob("*/attempt-*/stderr")
            for line in path.read_text(errors="replace").splitlines()
            if "warning" in line.lower()
        }
    )
    unreviewed = [line for line in warnings if line not in approved_warnings]
    if unreviewed:
        reasons.append("Unreviewed tool warnings")
    return {
        "status": "READY" if not reasons else "NEEDS_CORRECTION",
        "scope": "Software/output inspection only; not scientific approval",
        "result_sha256": digest(result_path),
        "wall_seconds": result.get("wall_seconds"),
        "reasons": reasons,
        "warnings": warnings,
        "unreviewed_warnings": unreviewed,
        "measurements": metrics,
    }


def verify_binding(binding: dict[str, str]) -> None:
    for name in ("target_spec", "target_dataset", "pilot_spec", "pilot_dataset"):
        if digest(Path(binding[name])) != binding[name + "_sha256"]:
            raise ValueError(f"Stale pilot: {name} changed")
