"""Regression checks for TSS mapping (workflows/tss-seq.nf).

Without options: stub runs and sheet checks. With --conda: real tools on a synthetic genome, with
CAGE reads (barcode, the cap's non-templated G, 3' adapter), TSS-seq reads (start at the TSS),
RAMPAGE reads (template-switch Gs) and Smar2C2 read pairs (primer, UMI and Gs in either mate,
PCR duplicates, oligo-dT reads), checking the exact tag positions, clusters and QC.
"""

from __future__ import annotations

import argparse
import csv
import functools
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from verify_pipeline import executor_limits, in_parallel, launch_folder, serve, trace_rows

ROOT = Path(__file__).resolve().parents[1]
HEADER = [
    "sample_id",
    "species",
    "read1_url",
    "read2_url",
    "method",
    "trim_5p",
    "adapter_3p",
    "reference_fasta",
    "reference_annotation",
    "organelle_contigs",
]
BARCODE = "ACGCAGCAG"
ADAPTER = "TCGTATGCCGTCTTCTGCTTG"
SMART = "AAGCAGTGGTATCAACGCAGAGTAC"
# (sample, method): CAGE carries a barcode, the cap G and an adapter; TSS-seq starts at the TSS;
# RAMPAGE adds template-switch Gs; Smar2C2 hides its TSS reads among Tn5 fragments in either mate.
LIBRARIES = (("cage", "CAGE"), ("tssseq", "TSS-seq"), ("rampage", "RAMPAGE"), ("smar", "Smar2C2"))
# (gene, strand, transcript start, transcript end), 1-based on chr1.
GENES = [
    ("geneA", "+", 20001, 22000),
    ("geneB", "-", 60001, 62000),
    ("geneC", "+", 120001, 122000),
    ("geneD", "-", 180001, 182000),
]
# Tags per offset from the annotated start, downstream positive.
PROFILE = {0: 40, -2: 6, -1: 8, 1: 10, 3: 5}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def write_sheet(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def create_fixture(work: Path, base_url: str, *, conda: bool) -> tuple[Path, dict[str, Counter]]:
    """References, reads and a sheet; returns each library's expected tags."""
    generator = random.Random(6031)
    chromosome = list("".join(generator.choices("ACGT", k=300000)))
    # The 3 bases before each TSS (on the transcript's strand) are not G, so the cap's added G and
    # the template-switch Gs always mismatch; a genomic G there is a known ambiguity.
    for _gene, strand, start, end in GENES:
        tss = start if strand == "+" else end
        for offset in PROFILE:
            site = tss + offset if strand == "+" else tss - offset
            for back in (1, 2, 3):
                before = site - back if strand == "+" else site + back
                chromosome[before - 1] = "A" if strand == "+" else "T"
    plastid = list("".join(generator.choices("ACGT", k=20000)))
    for index in range(12):
        for back in (1, 2, 3):
            plastid[1000 + 300 * index - 1 - back] = "A"
    contigs = {"chr1": "".join(chromosome), "ChrC": "".join(plastid)}
    references = work / "references"
    references.mkdir()
    with gzip.open(references / "plant.fa.gz", "wt") as stream:
        stream.writelines(f">{name}\n{sequence}\n" for name, sequence in contigs.items())
    lines = ["##gff-version 3\n"]
    for gene, strand, start, end in GENES:
        lines += [
            f"chr1\tsynthetic\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene}\n",
            f"chr1\tsynthetic\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={gene}.1;Parent={gene}\n",
            f"chr1\tsynthetic\texon\t{start}\t{end}\t.\t{strand}\t.\tParent={gene}.1\n",
        ]
    (references / "plant.gff3").write_text("".join(lines))

    def random_bases(length: int) -> str:
        return "".join(generator.choices("ACGT", k=length))

    def internal() -> str:
        # A Tn5 fragment from inside a gene: no primer, so Smar2C2 extraction skips it.
        start = generator.randrange(30000, 50000)
        return contigs["chr1"][start : start + 70]

    sites: list[tuple[str, str, int, int]] = []
    for _gene, strand, start, end in GENES:
        tss = start if strand == "+" else end
        for offset, count in PROFILE.items():
            sites.append(("chr1", strand, tss + offset if strand == "+" else tss - offset, count))
    # Plastid TSSs count as organelle tags.
    sites += [("ChrC", "+", 1000 + 300 * index, 1) for index in range(12)]
    expected: dict[str, Counter] = {}
    rows = []
    for sample, method in LIBRARIES:
        first: list[str] = []
        second: list[str] = []
        tags: Counter = Counter()
        for contig, strand, site, count in sites:
            if strand == "+":
                transcript = contigs[contig][site - 1 : site - 1 + 40]
            else:
                transcript = reverse_complement(contigs[contig][site - 40 : site])
            for copy in range(count):
                tags[(contig, strand, site)] += 1
                if method == "CAGE":
                    first.append(BARCODE + "G" + transcript[:30] + ADAPTER)
                elif method == "TSS-seq":
                    first.append(transcript)
                elif method == "RAMPAGE":
                    first.append("G" * (3 + copy % 3) + transcript)
                else:
                    # Each molecule is read twice (a PCR duplicate, same UMI), in either mate.
                    # Tn5 also cuts inside the oligo, so the primer can start far into the read.
                    code = random_bases(12)
                    spacer = random_bases(5 if copy % 3 else 60)
                    read = spacer + SMART + code + "AT" + "GGG" + transcript
                    for _ in range(2):
                        pair = [read, internal()]
                        if copy % 2:
                            pair.reverse()
                        first.append(pair[0])
                        second.append(pair[1])
        if method == "Smar2C2":
            # Oligo-dT (3' end) reads and plain Tn5 fragments give no tags.
            for _ in range(20):
                first.append(random_bases(40) + SMART + "T" * 30 + random_bases(20))
                second.append(internal())
        folder = work / "inputs" / sample
        folder.mkdir(parents=True)
        for name, reads in (("reads", first), ("mates", second)):
            if not reads:
                continue
            with gzip.open(folder / f"{name}.fastq.gz", "wt") as stream:
                for index, read in enumerate(reads):
                    stream.write(f"@{sample}{index}\n{read}\n+\n{'I' * len(read)}\n")
        expected[sample] = tags
        rows.append(
            {
                "sample_id": sample,
                "species": "Plant_synthetic",
                "read1_url": f"{base_url}/{sample}/reads.fastq.gz",
                "read2_url": f"{base_url}/{sample}/mates.fastq.gz" if second else "-",
                "method": method,
                "trim_5p": "9" if method == "CAGE" else "-",
                "adapter_3p": ADAPTER if method == "CAGE" else "-",
                "reference_fasta": "plant.fa.gz",
                "reference_annotation": "plant.gff3",
                "organelle_contigs": "ChrC",
            }
        )
    sheet = work / "samples.tsv"
    write_sheet(sheet, rows)
    limits = executor_limits(docker=False)
    if conda:
        limits += f"conda.cacheDir = '{work.parent / 'conda_cache'}'\n"
    (work / "limits.config").write_text(limits)
    return sheet, expected


def run_tss(
    work: Path, sheet: Path, *, conda: bool = False, expect_failure: bool = False
) -> subprocess.CompletedProcess[str]:
    command = [
        "nextflow",
        "-log",
        str(work / "nextflow-run.log"),
        "run",
        str(ROOT / "workflows/tss-seq.nf"),
        "-c",
        str(ROOT / "conf/tss-seq.config"),
        "-profile",
        "local,test,conda" if conda else "local,test",
        "--input",
        str(sheet),
        "--references",
        str(work / "references"),
        "--workspace",
        str(work / "published"),
        "-work-dir",
        str(work / "work"),
        "-c",
        str(work / "limits.config"),
        "-with-trace",
        str(work / "trace-run.tsv"),
        "-ansi-log",
        "false",
    ]
    if not conda:
        command.append("-stub-run")
    print(f"Running {'conda' if conda else 'stub'} TSS mapping in {work}", flush=True)
    environment = {
        **os.environ,
        "PATH": f"{Path(sys.executable).parent}:{os.environ['PATH']}",
        "NXF_OFFLINE": "true",
        "NXF_DISABLE_CHECK_LATEST": "true",
    }
    result = subprocess.run(
        command,
        cwd=work,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=3600,
        check=False,
    )
    (work / "output-run.log").write_text(result.stdout)
    if result.returncode and not expect_failure:
        print(result.stdout[-12000:], flush=True)
    return result


def check_trace(work: Path) -> None:
    rows = trace_rows(work)
    assert all(row["status"] == "COMPLETED" for row in rows), rows
    counts = Counter(row["name"].split(" ")[0] for row in rows)
    assert counts == {
        "GET_REFERENCES:CHROM_SIZES": 1,
        "GET_REFERENCES:BWA_MEM2_INDEX": 1,
        "RNA_ANNOTATION": 1,
        "GET_STAR_INDEX:RNA_INDEX": 1,
        "DOWNLOAD": 1,
        "FASTQC": 5,
        "TSS_EXTRACT": 1,
        "TSS_TRIM": 4,
        "TSS_ALIGN": 4,
        "TSS_CTSS": 4,
        "TSS_CLUSTERS": 4,
        "TSS_QC": 4,
        "MULTIQC_REPORT": 1,
    }, counts


def check_results(work: Path, expected: dict[str, Counter]) -> None:
    """Exact tags after trimming, G handling and UMI deduplication, clusters, and the QC table."""
    output = work / "published" / "samples" / "output" / "Plant_synthetic"
    for sample, tags in expected.items():
        observed: Counter = Counter()
        with gzip.open(output / sample / f"{sample}.ctss.tsv.gz", "rt") as stream:
            for line in stream:
                contig, strand, position, count = line.rstrip("\n").split("\t")
                observed[(contig, strand, int(position))] = int(count)
        assert observed == tags, (sample, sorted(set(observed) ^ set(tags))[:10])
        for track in ("plus", "minus", "plus.cpm", "minus.cpm"):
            assert (output / sample / f"{sample}.ctss.{track}.bw").stat().st_size > 0
        clusters = [
            line.split("\t")
            for line in (output / sample / f"{sample}.tss_clusters.bed").read_text().splitlines()
        ]
        nuclear = [c for c in clusters if c[0] == "chr1"]
        assert len(nuclear) == len(GENES), (sample, nuclear)
        for (_gene, strand, start, end), cluster in zip(
            sorted(GENES, key=lambda g: g[2]), sorted(nuclear, key=lambda c: int(c[1])), strict=True
        ):
            tss = start if strand == "+" else end
            assert cluster[5] == strand, cluster
            # The dominant TSS (0-based) is the annotated start, and the cluster holds every tag.
            assert int(cluster[7]) + 1 == tss, (sample, cluster)
            assert int(cluster[4]) == sum(PROFILE.values()), (sample, cluster)
        with (output / sample / "qc" / f"{sample}.tss_qc_mqc.tsv").open() as stream:
            rows = [line for line in stream if not line.startswith("#")]
        qc = next(csv.DictReader(rows, delimiter="\t"))
        total = sum(tags.values())
        organelle = sum(n for (contig, _s, _p), n in tags.items() if contig == "ChrC")
        assert int(qc["tags"]) == total, qc
        assert float(qc["near_TSS_100"]) == round((total - organelle) / total, 4), qc
        assert float(qc["antisense_500"]) == 0, qc
        assert float(qc["organelle"]) == round(organelle / total, 4), qc
        assert qc["rRNA"] == "NA", qc


def check_failures(work: Path, sheet: Path) -> None:
    """A rejected sheet stops before any task."""
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    lines = sheet.read_text().splitlines(keepends=True)
    cases: list[tuple[str, str, str]] = []
    for suffix, changed, diagnostic in (
        ("method", [{**rows[0], "method": "TIF-seq"}], "method must be"),
        ("trim", [{**rows[0], "trim_5p": "0"}], "trim_5p must be"),
        ("adapter", [{**rows[0], "adapter_3p": "ACG"}], "adapter_3p must be"),
        ("annotation", [{**rows[0], "reference_annotation": "absent.gff3"}], "annotation missing"),
        ("organelles", [rows[0], {**rows[1], "organelle_contigs": "-"}], "disagree on"),
        ("smar-mates", [{**rows[3], "read2_url": "-"}], "takes both read1_url and read2_url"),
        ("cage-mates", [{**rows[0], "read2_url": rows[3]["read2_url"]}], "takes one FASTQ"),
    ):
        candidate = work / f"{suffix}.tsv"
        write_sheet(candidate, changed)
        cases.append((suffix, str(candidate), diagnostic))
    extra = work / "extra-column.tsv"
    extra.write_text("".join(line.rstrip("\n") + "\tnote\n" for line in lines))
    cases.append(("extra-column", str(extra), "expected header"))

    def reject(suffix: str, candidate: str, diagnostic: str) -> None:
        case = launch_folder(work, suffix)
        assert run_tss(case, Path(candidate), expect_failure=True).returncode != 0
        assert trace_rows(case) == [], suffix
        log = (case / "output-run.log").read_text()
        assert diagnostic in log, (suffix, log[-2000:])

    in_parallel([functools.partial(reject, *case) for case in cases])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conda", action="store_true", help="Use real tools from environments/")
    parser.add_argument(
        "--keep", action="store_true", help="Retain fixtures, logs, and work directories"
    )
    args = parser.parse_args()
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(
        tempfile.mkdtemp(prefix=f"tss-{'conda' if args.conda else 'regression'}-", dir=runs)
    )
    passed = False
    try:
        (work / "inputs").mkdir()
        with serve(work / "inputs", docker=False) as base_url:
            sheet, expected = create_fixture(work, base_url, conda=args.conda)
            assert run_tss(work, sheet, conda=args.conda).returncode == 0
            check_trace(work)
            if args.conda:
                check_results(work, expected)
            else:
                check_failures(work, sheet)
                # Stub runs never publish into the references folder.
                assert sorted(p.name for p in (work / "references").iterdir()) == [
                    "plant.fa.gz",
                    "plant.gff3",
                ]
        label = "Real-tool (conda) TSS mapping" if args.conda else "TSS Nextflow regression"
        print(f"PASS: {label}", flush=True)
        passed = True
    finally:
        if passed and not args.keep:
            shutil.rmtree(work)
        else:
            print(f"Kept {work}", flush=True)


if __name__ == "__main__":
    main()
