"""Paired FASTQ reading for bulk ATAC lane checks."""

from __future__ import annotations

import gzip
import itertools
from collections.abc import Iterator
from pathlib import Path
from typing import Protocol


class LineReader(Protocol):
    def readline(self, size: int = -1, /) -> bytes: ...


def fastq(path: Path) -> Iterator[tuple[bytes, bytes]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as stream:
        yield from fastq_stream(stream)


def fastq_stream(stream: LineReader) -> Iterator[tuple[bytes, bytes]]:
    """Validate complete records from a binary FASTQ stream."""
    while header := stream.readline():
        sequence, plus, quality = stream.readline(), stream.readline(), stream.readline()
        if (
            not header.startswith(b"@")
            or not header[1:].split()
            or not plus.startswith(b"+")
            or not quality.endswith(b"\n")
            or not sequence.strip()
            or len(sequence.rstrip(b"\r\n")) != len(quality.rstrip(b"\r\n"))
        ):
            raise ValueError("Corrupt or truncated FASTQ record")
        name = header[1:].split()[0]
        if name.endswith((b"/1", b"/2")):
            name = name[:-2]
        if not name:
            raise ValueError("FASTQ template name is empty")
        yield name, header + sequence + plus + quality


def paired(read1: Path, read2: Path) -> Iterator[tuple[bytes, bytes, bytes]]:
    found = False
    for first, second in itertools.zip_longest(fastq(read1), fastq(read2)):
        if first is None or second is None or first[0] != second[0]:
            raise ValueError("FASTQ mates differ in record count or template names")
        found = True
        yield first[0], first[1], second[1]
    if not found:
        raise ValueError("Empty FASTQ library")
