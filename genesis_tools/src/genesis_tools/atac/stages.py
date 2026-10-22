"""ATAC task helpers executed inside pinned containers."""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
import os
import shutil
import subprocess
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

from .inputs import digest
from .reads import paired


class TemporaryDownloadError(RuntimeError):
    """A transport failure for which bounded acquisition retries are safe."""


def acquire(library: dict[str, Any]) -> None:
    evidence = []
    counts = 0
    with Path("raw_R1.fastq.gz").open("wb") as first, Path("raw_R2.fastq.gz").open("wb") as second:
        for lane in library["lanes"]:
            paths = []
            for mate in ("read1", "read2"):
                path = Path(f"{lane['lane_id']}.{mate}.fastq.gz")
                source = lane[mate]
                if source.startswith("https://"):
                    try:
                        with (
                            urllib.request.urlopen(source, timeout=60) as response,
                            path.open("wb") as target,
                        ):
                            shutil.copyfileobj(response, target, 1024 * 1024)
                    except urllib.error.HTTPError as error:
                        if error.code == 429 or 500 <= error.code < 600:
                            raise TemporaryDownloadError(
                                "Temporary download service failure"
                            ) from None
                        raise ValueError(
                            f"Download HTTP status {error.code}; verify the source"
                        ) from None
                    except urllib.error.URLError, TimeoutError:
                        raise TemporaryDownloadError(
                            "Download unavailable; retry acquisition after checking connectivity"
                        ) from None
                else:
                    shutil.copyfile(source, path)
                checksum = digest(path)
                if checksum != lane[mate + "_sha256"]:
                    raise ValueError("Downloaded FASTQ checksum does not match the manifest")
                paths.append(path)
                evidence.append({"lane": lane["lane_id"], "mate": mate, "sha256": checksum})
            lane_count = sum(1 for _ in paired(paths[0], paths[1]))
            if lane_count == 0:
                raise ValueError("FASTQ lane contains no paired records")
            counts += lane_count
            for source, target in zip(paths, (first, second), strict=True):
                with source.open("rb") as stream:
                    shutil.copyfileobj(stream, target, 1024 * 1024)
                source.unlink()
    Path("acquisition.json").write_text(
        json.dumps({"templates": counts, "sources": evidence}, indent=2) + "\n"
    )


def reference(fasta: Path, metadata: dict[str, Any]) -> None:
    import pysam

    if digest(fasta) != metadata["fasta_sha256"]:
        raise ValueError("Reference FASTA checksum changed")
    opener = gzip.open if fasta.suffix == ".gz" else open
    with opener(fasta, "rb") as source, Path("genome.fa").open("wb") as target:
        shutil.copyfileobj(source, target, 1024 * 1024)
    pysam.faidx("genome.fa")
    with pysam.FastaFile("genome.fa") as reader:
        sizes = dict(zip(reader.references, reader.lengths, strict=True))
    declared = set(metadata["mitochondrial_contigs"]) | set(metadata["plastid_contigs"])
    if declared - sizes.keys():
        raise ValueError("Reference is missing declared organellar contigs")
    Path("chrom.sizes").write_text("".join(f"{c}\t{n}\n" for c, n in sizes.items()))
    subprocess.run(["bwa-mem2", "index", "genome.fa"], check=True)
    Path("reference.json").write_text(json.dumps(metadata, indent=2) + "\n")


def finish_fragments(policy: dict[str, Any]) -> None:
    import pysam

    from .fragments import select

    select(Path("marked.bam"), policy, Path("."))
    for source, target in (
        ("fragments.unsorted.bed", "fragments.bed"),
        ("cuts.unsorted.bed", "cuts.bed"),
    ):
        with Path(target).open("w") as stream:
            subprocess.run(
                ["sort", "-S", "128M", "-T", ".", "-k1,1", "-k2,2n", "-k3,3n", source],
                stdout=stream,
                check=True,
                env={**os.environ, "LC_ALL": "C"},
            )
        Path(source).unlink()
    with Path("cuts.bed").open() as source, Path("cuts.bedgraph").open("w") as target:
        for line, items in itertools.groupby(source):
            target.write(line.rstrip() + "\t" + str(sum(1 for _ in items)) + "\n")
    pysam.tabix_compress("fragments.bed", "fragments.bed.gz", force=True)
    pysam.tabix_index("fragments.bed.gz", preset="bed", force=True)
    for tool, name in (
        (pysam.flagstat, "flagstat"),
        (pysam.stats, "stats"),
        (pysam.idxstats, "idxstats"),
    ):
        text = tool("usable.bam")
        if not isinstance(text, str):
            raise ValueError("Unexpected samtools output")
        Path(f"usable.{name}.txt").write_text(text)


def enrichment() -> None:
    from .metrics import Regions, tss_enrichment

    regions = []
    for line in Path("atac_peaks.narrowPeak").read_text().splitlines():
        fields = line.split("\t")
        regions.append((fields[0], int(fields[1]), int(fields[2])))
    index = Regions(regions)
    total = overlaps = 0
    with gzip.open("fragments.bed.gz", "rt") as source:
        for line in source:
            chrom, start, end, _ = line.rstrip().split("\t")
            total += 1
            overlaps += bool(index.hits(chrom, int(start), int(end)))
    sizes = {
        c: int(n)
        for c, n in (line.split() for line in Path("chrom.sizes").read_text().splitlines())
    }
    tss = tss_enrichment(Path("usable.bam"), Path("tss.bed"), sizes, Path("."))
    value = {
        "usable_fragments": total,
        "peak_count": len(regions),
        "fragments_in_peaks": overlaps,
        "FRiP": overlaps / total if total else None,
        "tss": tss,
        "plant_thresholds": "UNSPECIFIED",
        "metric_version": 1,
    }
    Path("enrichment.json").write_text(json.dumps(value, indent=2) + "\n")


def check_report() -> None:
    data = json.loads(Path("multiqc_data/multiqc_data.json").read_text())
    raw = data["report_saved_raw_data"]
    expected = {
        "multiqc_fastqc": 2,
        "multiqc_samtools_flagstat": 2,
        "multiqc_samtools_stats": 2,
        "multiqc_samtools_idxstats": 2,
    }
    found = {}
    for key, count in expected.items():
        if key not in raw or len(raw[key]) != count:
            raise ValueError(f"Missing or ambiguous MultiQC metrics: {key}")
        found[key] = len(raw[key])
    Path("multiqc_data/required-modules.json").write_text(json.dumps(found, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage", choices=["acquire", "reference", "fragments", "enrichment", "report"]
    )
    parser.add_argument("--manifest", type=Path, default=Path("library.json"))
    args = parser.parse_args()
    if args.stage == "enrichment":
        enrichment()
    elif args.stage == "report":
        check_report()
    else:
        value = json.loads(args.manifest.read_text())
        if args.stage == "acquire":
            try:
                acquire(value)
            except TemporaryDownloadError as error:
                parser.exit(75, f"{error}\n")
        elif args.stage == "reference":
            reference(Path("input.fa.gz"), value)
        else:
            finish_fragments(value)


if __name__ == "__main__":
    main()
