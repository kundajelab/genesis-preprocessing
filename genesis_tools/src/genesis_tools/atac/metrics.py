"""Interval and TSS-enrichment helpers for bulk ATAC (shared with the benchmark comparisons)."""

from __future__ import annotations

import bisect
import collections
import gzip
import hashlib
import itertools
import math
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pysam

Interval = tuple[str, int, int]


def metric(
    value: object, definition: str, unit: str, denominator: object = None, reason: str | None = None
) -> dict[str, Any]:
    return {
        "value": value,
        "definition": definition,
        "version": 1,
        "unit": unit,
        "denominator": denominator,
        "status": "OK" if value is not None else "UNDEFINED",
        "reason": reason if value is None else None,
    }


def intervals(path: Path) -> tuple[list[Interval], list[float]]:
    opener = gzip.open if path.suffix == ".gz" else open
    rows, scores = [], []
    with opener(path, "rt") as stream:
        for line in stream:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError("Interval input requires at least three tab-separated columns")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError("Invalid BED0 half-open interval")
            rows.append((chrom, start, end))
            scores.append(float(fields[8]) if len(fields) >= 9 else 0.0)
    if any(not math.isfinite(x) for x in scores):
        raise ValueError("Non-finite peak scores")
    return rows, scores


def merged(rows: list[Interval]) -> list[Interval]:
    result: list[Interval] = []
    for chrom, start, end in sorted(rows):
        if result and result[-1][0] == chrom and result[-1][2] >= start:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    return result


def intersection(left: list[Interval], right: list[Interval]) -> int:
    a, b = merged(left), merged(right)
    i = j = total = 0
    while i < len(a) and j < len(b):
        c, s, e = a[i]
        d, x, y = b[j]
        if c == d:
            total += max(0, min(e, y) - max(s, x))
        if c < d or (c == d and e <= y):
            i += 1
        else:
            j += 1
    return total


class Regions:
    """Index intervals; preserve individual region counts and merge only for union metrics."""

    def __init__(self, rows: list[Interval]) -> None:
        self.rows = rows
        self.by_chrom: dict[str, list[tuple[int, int, int]]] = collections.defaultdict(list)
        for i, (chrom, start, end) in enumerate(rows):
            self.by_chrom[chrom].append((start, end, i))
        self.starts = {}
        self.prefix_end = {}
        for chrom, entries in self.by_chrom.items():
            entries.sort()
            self.starts[chrom] = [e[0] for e in entries]
            self.prefix_end[chrom] = list(itertools.accumulate((e[1] for e in entries), max))

    def hits(self, chrom: str, start: int, end: int) -> list[int]:
        entries = self.by_chrom.get(chrom, [])
        hi = bisect.bisect_left(self.starts.get(chrom, []), end)
        lo = bisect.bisect_right(self.prefix_end.get(chrom, []), start, hi=hi)
        return [i for s, e, i in entries[lo:hi] if s < end and e > start]


def groups(bam: Path, temporary: Path) -> Iterator[list[pysam.AlignedSegment]]:
    """Name collate on disk; memory bounded by one template rather than the library."""
    collated = temporary / "names.bam"
    pysam.sort("-n", "-m", "128M", "-o", str(collated), str(bam))
    try:
        with pysam.AlignmentFile(str(collated)) as stream:
            for _, records in itertools.groupby(stream, key=lambda r: r.query_name):
                yield list(records)
    finally:
        collated.unlink(missing_ok=True)


def primary(records: list[pysam.AlignedSegment]) -> list[pysam.AlignedSegment]:
    return [r for r in records if not (r.is_unmapped or r.is_secondary or r.is_supplementary)]


def fragment(records: list[pysam.AlignedSegment], layout: str) -> Interval | None:
    if layout == "SE":
        if len(records) != 1:
            return None
    elif (
        len(records) != 2
        or not all(r.is_paired for r in records)
        or sum(r.is_read1 for r in records) != 1
        or sum(r.is_read2 for r in records) != 1
        or records[0].reference_id != records[1].reference_id
    ):
        return None
    starts = [r.reference_start for r in records]
    ends = [r.reference_end for r in records]
    if any(end is None for end in ends):
        return None
    chrom = records[0].reference_name
    if chrom is None:
        return None
    return chrom, min(starts), max(e for e in ends if e is not None)


def tn5_cuts(start: int, end: int) -> tuple[int, int]:
    """BED0 cut-base coordinates corresponding to shifted interval [start+4,end-5)."""
    if end - start < 10:
        raise ValueError("Fragment too short for +4/-5 interval")
    return start + 4, end - 6


def tss_enrichment(
    bam: Path, tss_file: Path, chrom_sizes: dict[str, int], temporary: Path
) -> dict[str, Any]:
    """Common ATAC profile: BED0 cut bases, +/-2000, 100-bp flank mean, center bin."""
    window = 2000
    tsses = []
    excluded = 0
    for line in tss_file.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        chrom, position, strand = line.split("\t")
        pos = int(position)
        if strand not in {"+", "-"} or chrom not in chrom_sizes or pos < 0:
            raise ValueError("Invalid three-column BED0 TSS point/strand")
        if pos < window or pos + window >= chrom_sizes[chrom]:
            excluded += 1
        else:
            tsses.append((chrom, pos, strand))
    windows = Regions([(chrom, pos - window, pos + window + 1) for chrom, pos, _ in tsses])
    profile = [0] * (2 * window + 1)
    short = 0
    for records in groups(bam, temporary):
        span = fragment(primary(records), "PE")
        if span is None:
            continue
        chrom, start, end = span
        if end - start < 10:
            short += 1
            continue
        for cut in tn5_cuts(start, end):
            for i in windows.hits(chrom, cut, cut + 1):
                _, pos, strand = tsses[i]
                offset = cut - pos if strand == "+" else pos - cut
                profile[offset + window] += 1
    background = (sum(profile[:100]) + sum(profile[-100:])) / 200
    value = profile[window] / background if tsses and background else None
    return {
        "metric": metric(
            value,
            "Tn5 BED0 cuts start+4/end-6; +/-2000; center bin divided by mean "
            "100bp flanks; full windows only",
            "fold_enrichment",
            background,
            "no full TSS windows or zero flank background",
        ),
        "profile": profile,
        "usable_TSS": len(tsses),
        "edge_excluded_TSS": excluded,
        "too_short_fragments": short,
        "annotation_sha256": hashlib.sha256(tss_file.read_bytes()).hexdigest(),
    }
