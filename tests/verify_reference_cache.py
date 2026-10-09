"""Known-answer reference provenance checks, without external tools or a framework."""

from __future__ import annotations

import gzip
import json
import os
import tempfile
from pathlib import Path

from genesis_tools.reference_cache import products, record, sha256, sidecar, verify


def check(work: Path) -> None:
    fasta = work / "tiny.fa.gz"
    a, b = b">chr1\nACGTACGT\n", b">chr1\nTTTTACGT\n"
    fasta.write_bytes(gzip.compress(a, compresslevel=0, mtime=0))
    sizes = work / "tiny.chrom.sizes"
    index = work / "tiny.bwa-mem2"
    container, recipe = "example@sha256:" + "a" * 64, "b" * 64

    def verify_all(mode: str = "stub") -> None:
        for artifact, kind in ((sizes, "chrom-sizes"), (index, "bwa-mem2")):
            verify(artifact, kind, sha256(fasta, decompress=True), mode, container, recipe)

    def rejected() -> None:
        try:
            verify_all()
        except ValueError as error:
            assert "fresh reference directory" in str(error)
        else:
            raise AssertionError("Unverified reference was accepted")

    verify_all()  # Absent products can be built.
    sizes.write_text("chr1\t8\n")
    rejected()  # Legacy product is not proof of identity.
    index.mkdir()
    for path in products(index, "bwa-mem2"):
        path.write_bytes(b"known artifact\n")
    for artifact, kind in ((sizes, "chrom-sizes"), (index, "bwa-mem2")):
        record(fasta, artifact, kind, "stub", container, recipe)
    verify_all()
    initial = {p: p.read_bytes() for p in work.rglob("*") if p.is_file()}
    before = fasta.stat()
    os.utime(fasta, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
    verify_all()  # An mtime is not an assembly identity.
    fasta.write_bytes(gzip.compress(a, mtime=1))
    verify_all()  # Gzip metadata is not FASTA content.
    fasta.write_bytes(gzip.compress(b, compresslevel=0, mtime=0))
    os.utime(fasta, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert fasta.stat().st_size == before.st_size
    rejected()  # Content differs even with the original timestamp.
    fasta.write_bytes(initial[fasta])
    renamed = work / "renamed.fa.gz"
    renamed.write_bytes(initial[fasta])
    verify(
        sizes,
        "chrom-sizes",
        sha256(renamed, decompress=True),
        "stub",
        container,
        recipe,
    )
    damaged = products(index, "bwa-mem2")[0]
    damaged.write_bytes(b"damaged\n")
    rejected()
    damaged.write_bytes(initial[damaged])
    damaged.unlink()
    rejected()
    damaged.write_bytes(initial[damaged])
    extra = index / "genome.alt"
    extra.write_text("unrecorded auxiliary index\n")
    rejected()
    extra.unlink()
    manifest = sidecar(sizes, "chrom-sizes")
    manifest.write_text("{not json")
    rejected()
    manifest.write_bytes(initial[manifest])
    for field, replacement in (("schema_version", 2), ("fasta_sha256", "0" * 64)):
        obj = json.loads(initial[manifest])
        obj[field] = replacement
        manifest.write_text(json.dumps(obj))
        rejected()
        manifest.write_bytes(initial[manifest])
    for altered_container, altered_recipe in (
        ("different", recipe),
        (container, "different"),
    ):
        try:
            verify(
                sizes,
                "chrom-sizes",
                sha256(fasta, decompress=True),
                "stub",
                altered_container,
                altered_recipe,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("Different generator accepted")
    try:
        verify_all("real")
    except ValueError:
        pass
    else:
        raise AssertionError("Stub reference accepted for a real analysis")
    verify_all()
    assert all(path.read_bytes() == content for path, content in initial.items())
    print("PASS: reference identity, metadata-only changes, corruption, generator and mode checks")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory(prefix="reference-cache-") as directory:
        check(Path(directory))
