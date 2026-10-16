from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from evidence import (
    Evidence,
    digest,
    number,
    read_table,
    tsv,
    unique_rows,
    verify_sources,
)
from reports import escape, fmt, reports


class ReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.e = Evidence()

    def test_completed_cases_and_actual_depths(self) -> None:
        self.assertEqual(len(self.e.cases), 64)
        sample = next(
            s
            for s in self.e.samples
            if s["cohort"] == "dap-pe-5m" and s["id"] == "SRR31238589"
        )
        self.assertEqual(sample["raw_templates"], 1657793)
        self.assertEqual(sample["control"], "SRR31238993")
        sample = next(
            s
            for s in self.e.samples
            if s["cohort"] == "dap-se-5m" and s["id"] == "SRR2926068"
        )
        self.assertEqual(sample["raw_templates"], 4241389)

    def test_export_preserves_empty_fields_without_trailing_whitespace(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "data.tsv"
            rows = [
                {"value": "0", "reason": ""},
                {"value": "", "reason": "not assessed"},
            ]
            tsv(path, rows)
            self.assertEqual(read_table(path), rows)
            for line in path.read_text().splitlines():
                self.assertEqual(line, line.rstrip())

    def test_missing_is_not_zero(self) -> None:
        self.assertIsNone(number(""))
        self.assertEqual(number("0"), 0)
        self.assertEqual(fmt(None), "nd")
        with self.assertRaises(ValueError):
            number("nan")
        self.assertIsNone(
            self.e.q(
                "atac-arabidopsis-1m", "genrich", "SRR4000468", "replicate_peak_jaccard"
            )
        )
        self.assertEqual(
            self.e.q("atac-arabidopsis-1m", "genrich", "SRR4000468", "peak_count"), 0
        )

    def test_counting_denominators_and_no_mismatches(self) -> None:
        self.assertTrue(self.e.counts)
        for row in self.e.counts:
            if row["status"] == "OK":
                self.assertEqual(float(row["value"]), 0)
                self.assertGreater(int(row["denominator"]), 0)

    def test_duplicate_observations_rejected(self) -> None:
        with self.assertRaises(ValueError):
            unique_rows([{"id": "a"}, {"id": "a"}], ("id",))

    def test_missing_and_stale_sources_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "a"
            source.write_text("original")
            catalog = root / "catalog.json"
            catalog.write_text(
                json.dumps({"sources": [{"path": "a", "sha256": digest(source)}]})
            )
            verify_sources(root, catalog)
            source.write_text("changed")
            with self.assertRaises(ValueError):
                verify_sources(root, catalog)
            source.unlink()
            with self.assertRaises(FileNotFoundError):
                verify_sources(root, catalog)

    def test_incompatible_metric_rejected(self) -> None:
        import evidence

        original = evidence.read_table

        def changed(path: Path) -> list[dict[str, str]]:
            rows = original(path)
            if path.name == "quality-summary.tsv":
                rows[0]["unit"] = "incompatible units"
            return rows

        with (
            patch("evidence.read_table", side_effect=changed),
            self.assertRaisesRegex(ValueError, "Incompatible metric"),
        ):
            Evidence()

    def test_absent_arial_rejected(self) -> None:
        from plots import configure

        with (
            tempfile.TemporaryDirectory() as tmp,
            self.assertRaisesRegex(ValueError, "Actual Arial"),
        ):
            configure(Path(tmp))

    def test_later_reference_evidence_preserves_exclusions(self) -> None:
        record = self.e.local["reference_5m"]
        self.assertEqual(record["status"], "COMPLETE")
        for row in record["libraries"].values():
            for side in ("left_mapping", "right_mapping"):
                mapping = row[side]
                self.assertEqual(
                    sum(mapping["counts"].values()), mapping["input_peaks"]
                )
                self.assertLess(
                    mapping["counts"]["UNIQUE_RECIPROCAL"], mapping["input_peaks"]
                )

    def test_output_scope_and_escaping(self) -> None:
        before = copy.deepcopy(self.e.samples)
        bodies = reports(self.e)
        self.assertEqual(len(bodies), 7)
        self.assertEqual(self.e.samples, before)
        self.assertIn("not biological accuracy", bodies["interval-counting"])
        self.assertIn("nd", bodies["interval-counting"])
        self.assertIn("existing recorded results", bodies["references"])
        self.assertEqual(escape("x_%&"), r"x\_\%\&")
        for body in bodies.values():
            self.assertNotIn("/data/homes/", body)
            self.assertNotIn("\u2014", body)
            self.assertNotIn("\u2013", body)


if __name__ == "__main__":
    unittest.main()
