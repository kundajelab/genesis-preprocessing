"""Deterministic synthetic SE/PE reads, reference and annotations for the ATAC checks.

Adapted from the benchmark fixture generator (now in genesis-curation): data files only.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import random
from pathlib import Path
from typing import Any

import pysam


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def create(out: Path, *, assay: str = "bulk-ATAC", small: bool = False) -> Path:
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    generator = random.Random(60141)
    length = 2000000 if not small else 5000
    sequence = "".join(generator.choices("ACGT", k=length))
    sequence = sequence[: length - 1000] + sequence[1000:1500] + sequence[length - 500 :]
    fasta = out / "reference.fa.gz"
    content = f">chr1\n{sequence}\n>plastid\n{'A' * 1000}\n>mitochondria\n{'C' * 1000}\n"
    fasta.write_bytes(gzip.compress(content.encode(), mtime=0))
    sizes = out / "chrom.sizes"
    sizes.write_text(f"chr1\t{length}\nplastid\t1000\nmitochondria\t1000\n")
    tss = out / "tss.bed"
    tss.write_text("chr1\t104\t+\nchr1\t194\t-\n")
    genes = out / "genes.bed"
    genes.write_text("chr1\t100\t500\ttoy1\t0\t+\nchr1\t800\t1200\ttoy2\t0\t-\n")
    gtf = out / "genes.gtf"
    gtf.write_text(
        'chr1\ttoy\texon\t101\t500\t.\t+\t.\tgene_id "toy1"; transcript_id "toy1";\n'
        'chr1\ttoy\texon\t801\t1200\t.\t-\t.\tgene_id "toy2"; transcript_id "toy2";\n'
    )
    assets = {}

    def asset(name: str, path: Path) -> str:
        assets[name] = {"path": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
        return name

    reference = {
        "id": "synthetic-v1",
        "species": "synthetic",
        "assembly": "synthetic-v1",
        "cultivar": "NOT_APPLICABLE",
        "fasta": asset("fasta", fasta),
        "fasta_content_sha256": hashlib.sha256(content.encode()).hexdigest(),
        "chrom_sizes": asset("sizes", sizes),
        "genome_size": length + 2000,
        "mitochondrial": ["mitochondria"],
        "plastid": ["plastid"],
        "tss": asset("tss", tss),
        "annotation_release": "synthetic-v1",
    }
    asset("genes", genes)
    reference["gtf"] = asset("gtf", gtf)
    libraries: list[dict[str, Any]] = []
    names = (
        ["se_control", "se_treatment", "se_second", "pe_control", "pe_treatment"]
        if assay == "DAP-seq"
        else ["pe_treatment"]
    )
    for name in names:
        layout = "SE" if name.startswith("se") else "PE"
        control = name.endswith("control")
        total = 12 if small else (6000 if control else 32000)
        header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [
                    {"SN": "chr1", "LN": length},
                    {"SN": "plastid", "LN": 1000},
                    {"SN": "mitochondria", "LN": 1000},
                ],
                "RG": [{"ID": name, "SM": name}],
            }
        )
        records = []
        fastqs: list[list[str]] = [[], []]
        last_start = 100
        for i in range(total):
            if small:
                start = 100 + (i % 4) * 300
            elif control or i >= 30000:
                start = generator.randrange(500, length - 500)
            else:
                # Unequal peak depths deliberately avoid the old nearly flat RPM fixture.
                peak = int((i / 30000) ** 1.3 * 200)
                start = 20000 + peak * 8000 + generator.randrange(150)
            duplicate = i % 11 == 1
            if duplicate:
                start = last_start
            last_start = start
            chrom = 1 if i == 2 else (2 if i == 3 else 0)
            if chrom:
                start = 100
            span = 180 if small else 180 + generator.randrange(-20, 21)
            q = [0, 5, 10, 20, 30, 60][i % 6] if small or i < 300 else 60
            genome = sequence if chrom == 0 else ("A" if chrom == 1 else "C") * 1000
            for mate in range(2 if layout == "PE" else 1):
                read = pysam.AlignedSegment(header)
                read.query_name = f"{name}_{i:08d}"
                reverse = bool(mate) if layout == "PE" else bool(i % 2)
                read.flag = (99 if mate == 0 else 147) if layout == "PE" else (16 if reverse else 0)
                read.reference_id = chrom
                read.reference_start = start + (span - 75 if mate else 0)
                read.mapping_quality = q if mate == 1 or layout == "SE" else 60
                read.cigarstring = "75M"
                read.query_sequence = genome[read.reference_start : read.reference_start + 75]
                read.query_qualities = pysam.qualitystring_to_array("I" * 75)
                read.set_tag("RG", name)
                if layout == "PE":
                    read.next_reference_id = chrom
                    read.next_reference_start = start + (span - 75 if mate == 0 else 0)
                    read.template_length = span if mate == 0 else -span
                records.append((read, duplicate))
                raw = read.query_sequence
                assert raw is not None
                if reverse:
                    raw = raw.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                fastqs[mate].append(f"@{read.query_name}\n{raw}\n+\n{'I' * 75}\n")
        raw_path, marked_path = out / f"{name}.bam", out / f"{name}.marked.bam"
        records.sort(key=lambda entry: (entry[0].reference_id, entry[0].reference_start))
        with (
            pysam.AlignmentFile(str(raw_path), "wb", header=header) as raw_stream,
            pysam.AlignmentFile(str(marked_path), "wb", header=header) as marked_stream,
        ):
            for read, duplicate in records:
                raw_stream.write(read)
                if duplicate:
                    read.flag |= 1024
                marked_stream.write(read)
        for path in (raw_path, marked_path):
            pysam.index(str(path))
        lib = {
            "id": name,
            "layout": layout,
            "control": None if control or assay == "bulk-ATAC" else name[:2] + "_control",
            "raw_templates": total,
            "biological_replicate": "synthetic",
            "bam": asset(name + "-bam", raw_path),
            "marked_bam": asset(name + "-marked", marked_path),
            "marking_method": "synthetic injected duplicate labels; no PCR inference",
            "tissue": "NOT_APPLICABLE",
            "assay_subtype": assay,
        }
        for mate in range(2 if layout == "PE" else 1):
            path = out / f"{name}.R{mate + 1}.fastq.gz"
            path.write_bytes(gzip.compress("".join(fastqs[mate]).encode(), mtime=0))
            lib[f"read{mate + 1}"] = asset(name + f"-R{mate + 1}", path)
        libraries.append(lib)
    dataset = {
        "schema_version": 1,
        "id": "synthetic-policy-v1",
        "assay": assay,
        "synthetic": True,
        "status": "READY",
        "reference": reference,
        "assets": assets,
        "libraries": libraries,
        "thresholds": "UNSPECIFIED",
        "description": (
            "Assigned MAPQs/duplicate labels test policy semantics, not aligner accuracy. "
            "Seed 60141. Unequal peak depths, two organelles and a repeated block."
        ),
    }
    write_json(out / "dataset.json", dataset)
    return out
