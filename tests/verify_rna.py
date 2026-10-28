"""RNA-seq Nextflow regressions and optional real-tool validation; no test framework."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Callable
from pathlib import Path

from verify_pipeline import INDEX_SUFFIXES, reference_record, serve, trace_rows, write_reads

ROOT = Path(__file__).resolve().parents[1]
HEADER = [
    "sample_id",
    "species",
    "read1_url",
    "read2_url",
    "strandedness",
    "reference_fasta",
    "reference_annotation",
    "organelle_contigs",
]
READ_LENGTH = 75
# Plant-like genes on chr1: short introns, one 3 kb intron (the longest annotated), both strands.
# Each is (gene, strand, exons as 0-based half-open intervals).
GENES = [
    ("AT1G00010", "+", [(10000, 10400), (10520, 10900), (11000, 11500)]),
    ("AT1G00020", "-", [(30000, 30500), (30600, 31100)]),
    ("AT1G00030", "+", [(50000, 50400), (53400, 54000)]),
    ("AT1G00040", "-", [(70000, 70700)]),
    # Identical sequence at both loci, like homeologs in a polyploid: every read multi-maps.
    ("HOMEOLOG_A", "+", [(90000, 90600)]),
    ("HOMEOLOG_B", "+", [(110000, 110600)]),
    ("RRNA1", "+", [(130000, 130800)]),
    ("ATCG00010", "+", [(2000, 2600)]),
]
# Liftoff-style extra copy: the same transcript ID at a second locus, with its own sequence.
COPIES = {"AT1G00050": [(150000, 150600), (170000, 170600)]}
# Reads spliced across this unannotated 8 kb intron must not align spliced: 8 kb is longer than
# any annotated intron.
NOVEL = [(200000, 200300), (208300, 208600)]
EXPRESSION = {
    "AT1G00010": 400,
    "AT1G00020": 200,
    "AT1G00030": 300,
    "AT1G00040": 100,
    "HOMEOLOG_A": 150,
    "RRNA1": 250,
    "ATCG00010": 120,
    "AT1G00050": 160,
    "NOVEL": 80,
}


def write_sheet(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_annotation(path: Path, *, gff3: bool) -> None:
    """GFF3 with Parent links, or GTF-style attributes, as plant annotations come."""
    lines = []
    genes = [(gene, strand, exons, "chr1") for gene, strand, exons in GENES]
    genes = [(g, s, e, "ChrC" if g.startswith("ATCG") else c) for g, s, e, c in genes]
    for gene, strand, exons, contig in genes:
        kind = "rRNA" if gene.startswith("RRNA") else "mRNA"
        start, end = exons[0][0] + 1, exons[-1][1]
        if gff3:
            lines.append(f"{contig}\tsyn\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene}\n")
            lines.append(
                f"{contig}\tsyn\t{kind}\t{start}\t{end}\t.\t{strand}\t.\tID={gene}.1;Parent={gene}\n"
            )
            lines += [
                f"{contig}\tsyn\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tParent={gene}.1\n"
                for s, e in exons
            ]
        else:
            ids = f'gene_id "{gene}"; transcript_id "{gene}.1";'
            lines.append(f"{contig}\tsyn\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{ids}\n")
            lines += [
                f"{contig}\tsyn\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{ids}\n" for s, e in exons
            ]
    for gene, loci in COPIES.items():
        if gff3:
            start, end = loci[0][0] + 1, loci[-1][1]
            lines.append(f"chr1\tsyn\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}\n")
        for copy, (s, e) in enumerate(loci, 1):
            if gff3:
                lines.append(
                    f"chr1\tsyn\tmRNA\t{s + 1}\t{e}\t.\t+\t.\tID={gene}.{copy};Parent={gene}\n"
                )
                lines.append(f"chr1\tsyn\texon\t{s + 1}\t{e}\t.\t+\t.\tParent={gene}.{copy}\n")
            else:
                # As Liftoff writes them (Col-CEN): each extra copy keeps the gene's transcript ID.
                ids = f'transcript_id "{gene}.1"; gene_id "{gene}"; Parent "{gene}";'
                gene_ids = f'transcript_id "{gene}"; gene_id "{gene}";'
                lines.append(f"chr1\tsyn\tgene\t{s + 1}\t{e}\t.\t+\t.\t{gene_ids}\n")
                lines.append(f'chr1\tsyn\texon\t{s + 1}\t{e}\t.\t+\t.\t{ids} ID "exon_1";\n')
                lines.append(f'chr1\tsyn\tmRNA\t{s + 1}\t{e}\t.\t+\t.\t{ids} ID "mRNA_1";\n')
    path.write_text("".join(lines))


def transcripts(contigs: dict[str, str]) -> dict[str, list[str]]:
    """Each expressed unit's spliced sense sequences."""
    result: dict[str, list[str]] = {}
    for gene, strand, exons in GENES:
        contig = "ChrC" if gene.startswith("ATCG") else "chr1"
        sequence = "".join(contigs[contig][s:e] for s, e in exons)
        result[gene] = [reverse_complement(sequence) if strand == "-" else sequence]
    for gene, loci in COPIES.items():
        result[gene] = [contigs["chr1"][s:e] for s, e in loci]
    result["NOVEL"] = ["".join(contigs["chr1"][s:e] for s, e in NOVEL)]
    return result


def create_fixture(work: Path, base_url: str, *, engine: str | None) -> Path:
    generator = random.Random(2718)
    chromosome = list("".join(generator.choices("ACGT", k=300000)))
    start, end = GENES[4][2][0]
    copy_start = GENES[5][2][0][0]
    chromosome[copy_start : copy_start + end - start] = chromosome[start:end]
    contigs = {
        "chr1": "".join(chromosome),
        "ChrC": "".join(generator.choices("ACGT", k=20000)),
    }
    references = work / "references"
    references.mkdir()
    for name in ("plant.fa.gz", "lifted.fa.gz"):
        with gzip.open(references / name, "wt") as stream:
            stream.writelines(f">{contig}\n{sequence}\n" for contig, sequence in contigs.items())
    write_annotation(references / "plant.gff3", gff3=True)
    write_annotation(references / "lifted.gff3", gff3=False)
    units = transcripts(contigs)
    rows = []
    libraries = (
        ("rep1", "reverse", "PE", "plant"),
        ("rep2", "forward", "SE", "plant"),
        ("lifted", "unstranded", "PE", "lifted"),
    )
    for sample, strandedness, layout, reference in libraries:
        first, second = [], []
        for unit, count in EXPRESSION.items():
            for index in range(count):
                sequence = units[unit][index % len(units[unit])]
                insert = generator.randrange(180, min(300, len(sequence)))
                position = generator.randrange(len(sequence) - insert + 1)
                fragment = sequence[position : position + insert]
                sense = fragment[:READ_LENGTH]
                antisense = reverse_complement(fragment)[:READ_LENGTH]
                # dUTP (reverse) libraries read the antisense strand first; unstranded ones either.
                flip = strandedness == "reverse" or (strandedness == "unstranded" and index % 2)
                first.append(antisense if flip else sense)
                second.append(sense if flip else antisense)
        folder = work / "inputs" / sample
        (folder / "mate").mkdir(parents=True)
        write_reads(folder / "reads.fastq.gz", first)
        if layout == "PE":
            write_reads(folder / "mate" / "reads.fastq.gz", second)
        rows.append(
            {
                "sample_id": sample,
                "species": "Plant_synthetic",
                "read1_url": f"{base_url}/{sample}/reads.fastq.gz",
                "read2_url": f"{base_url}/{sample}/mate/reads.fastq.gz" if layout == "PE" else "-",
                "strandedness": strandedness,
                "reference_fasta": f"{reference}.fa.gz",
                "reference_annotation": f"{reference}.gff3",
                "organelle_contigs": "ChrC" if reference == "plant" else "-",
            }
        )
    sheet = work / "samples.tsv"
    write_sheet(sheet, rows)
    limits = "executor.cpus = 4\nexecutor.memory = '8 GB'\n"
    if engine == "docker":
        # Linux Docker only resolves host.docker.internal with an explicit host-gateway mapping.
        limits += (
            "docker.runOptions = '-u $(id -u):$(id -g) "
            "--add-host=host.docker.internal:host-gateway'\n"
        )
    if engine == "conda":
        # Keep the environments built for earlier runs.
        limits += f"conda.cacheDir = '{work.parent / 'conda_cache'}'\n"
    (work / "limits.config").write_text(limits)
    return sheet


def run_rna(
    work: Path,
    sheet: Path,
    *,
    engine: str | None = None,
    resume: bool = False,
    suffix: str | None = None,
    expect_failure: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run the actual workflow; write full logs for actionable failure diagnostics."""
    suffix = suffix or ("resume" if resume else "run")
    command = [
        "nextflow",
        "-log",
        str(work / f"nextflow-{suffix}.log"),
        "run",
        str(ROOT / "workflows/rna-seq.nf"),
        "-c",
        str(ROOT / "conf/rna-seq.config"),
        "-profile",
        f"local,test,{engine}" if engine else "local,test",
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
        str(work / f"trace-{suffix}.tsv"),
        "--download_batch_size",
        "2",
        "-ansi-log",
        "false",
    ]
    if not engine:
        command.append("-stub-run")
    if resume:
        command.append("-resume")
    print(f"Running {engine or 'stub'} RNA-seq ({suffix}) in {work}", flush=True)
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
    (work / f"output-{suffix}.log").write_text(result.stdout)
    if result.returncode and not expect_failure:
        print(result.stdout[-12000:], flush=True)
    return result


def read_table(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def check_outputs(work: Path, *, engine: str | None) -> None:
    """Assert the process wiring, publication and recorded provenance."""
    rows = trace_rows(work)
    assert all(row["status"] == "COMPLETED" for row in rows), rows
    counts = Counter(row["name"].split(" ")[0] for row in rows)
    assert counts == {
        "GET_REFERENCES:CHROM_SIZES": 2,
        "GET_REFERENCES:BWA_MEM2_INDEX": 2,
        "RNA_ANNOTATION": 2,
        "GET_STAR_INDEX:RNA_INDEX": 2,
        "DOWNLOAD": 2,
        "FASTQC": 5,
        "RNA_ALIGN": 3,
        "RNA_COUNTS": 3,
        "RNA_TRACKS": 3,
        "RNA_QC": 3,
        "MULTIQC_REPORT": 1,
    }, counts
    published = work / "published" / "samples" / "output"
    for path in published.rglob("*"):
        assert not path.name.endswith((".fastq", ".fastq.gz", ".bam", ".bai")), path
    for sample in ("rep1", "rep2", "lifted"):
        library = published / "Plant_synthetic" / sample
        for name in (
            f"{sample}.featureCounts.tsv",
            f"{sample}.cpm.bw",
            f"qc/{sample}.featureCounts.tsv.summary",
            f"qc/{sample}.star.Log.final.out",
            f"qc/{sample}.star.SJ.out.tab",
            f"qc/{sample}.rna_qc_mqc.tsv",
            f"qc/fastqc/{sample}.read1_fastqc.zip",
        ):
            assert (library / name).is_file(), library / name
    assert (published / "multiqc/multiqc_report.html").stat().st_size > 0
    assert (published / "dap-images.lock").read_bytes() == (
        ROOT / "conf/dap-images.lock"
    ).read_bytes()
    status = (published / "git-status.txt").read_text().splitlines()
    assert status[0].startswith("branch: ") and status[1].startswith("commit: "), status
    log = (work / "output-run.log").read_text()
    assert "PIPELINE SUCCEEDED" in log, log[-2000:]
    if engine:
        check_results(published)


def check_results(published: Path) -> None:
    """Known answers from the synthetic libraries, through the real tools."""
    total = sum(EXPRESSION.values())
    for sample, forward in (("rep1", 0), ("rep2", 1), ("lifted", 0.5)):
        library = published / "Plant_synthetic" / sample
        lines = (library / f"{sample}.featureCounts.tsv").read_text().splitlines()[1:]
        table = csv.DictReader(lines, delimiter="\t")
        genes = {row["Geneid"]: int(row[f"{sample}.bam"]) for row in table}
        # Unique genes keep every fragment, including both Liftoff copies of AT1G00050.
        for gene in (
            "AT1G00010",
            "AT1G00020",
            "AT1G00030",
            "AT1G00040",
            "RRNA1",
            "ATCG00010",
        ):
            assert genes[gene] == EXPRESSION[gene], (sample, gene, genes[gene])
        assert genes["AT1G00050"] == EXPRESSION["AT1G00050"], (sample, genes)
        # Homeolog reads align to both copies, so neither counts them.
        assert genes["HOMEOLOG_A"] == genes["HOMEOLOG_B"] == 0, (sample, genes)
        (qc,) = read_table(library / f"qc/{sample}.rna_qc_mqc.tsv")
        assert abs(float(qc["forward_fraction"]) - forward) < 0.1, (sample, qc)
        assert int(qc["input_reads"]) == total, qc
        assert abs(float(qc["multimapped"]) - EXPRESSION["HOMEOLOG_A"] / total) < 0.01, qc
        assigned = sum(EXPRESSION.values()) - EXPRESSION["HOMEOLOG_A"] - EXPRESSION["NOVEL"]
        assert abs(float(qc["rRNA"]) - EXPRESSION["RRNA1"] / assigned) < 0.01, (
            sample,
            qc,
        )
        if sample == "lifted":
            assert float(qc["organelle"]) == 0, qc
        else:
            assert float(qc["organelle"]) > 0.05, qc
        # No junction longer than the longest annotated intron (3 kb), so the 8 kb one is absent.
        junctions = (library / f"qc/{sample}.star.SJ.out.tab").read_text().splitlines()
        lengths = [int(line.split("\t")[2]) - int(line.split("\t")[1]) + 1 for line in junctions]
        assert lengths and max(lengths) == 3000, (sample, sorted(lengths))
        assert qc["max_intron"] == "3000", qc


def check_resume(work: Path, sheet: Path, *, engine: str | None) -> None:
    """Resumed runs reuse the downloads, and then every task."""
    assert run_rna(work, sheet, engine=engine, resume=True).returncode == 0
    rows = trace_rows(work, "resume")
    tasks = [row for row in rows if row["name"].split(" ")[0] == "DOWNLOAD"]
    assert tasks and all(row["status"] == "CACHED" for row in tasks), tasks
    indexes = [row for row in rows if row["name"].split(" ")[0] == "GET_STAR_INDEX:RNA_INDEX"]
    if engine:
        # Published STAR indexes are reused from disk without any index task.
        assert not indexes, indexes
    else:
        # Stub runs publish no indexes, so they are rebuilt, from Nextflow's task cache.
        assert len(indexes) == 2 and all(row["status"] == "CACHED" for row in indexes), indexes
    # Newly published chromosome sizes and STAR indexes change staged input paths on the first
    # resume.
    # A second resume must cache all downstream computation at those stable paths.
    assert run_rna(work, sheet, engine=engine, resume=True, suffix="stable").returncode == 0
    rows = trace_rows(work, "stable")
    assert rows and all(row["status"] == "CACHED" for row in rows), rows


STAR_PIN = re.search(r"star=([^\s,<>]+)", (ROOT / "environments/GENESIS_STAR.yaml").read_text())
assert STAR_PIN, "environments/GENESIS_STAR.yaml does not pin a STAR version"
STAR_VERSION = STAR_PIN.group(1)
# Header lines of a STAR index record: FASTA, GTF, STAR version and three index parameters.
STAR_HEADER = 6


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def task_dir(work: Path, name: str) -> Path:
    """A traced task's work folder, from its hash prefix."""
    (row,) = [row for row in trace_rows(work) if row["name"] == name]
    prefix, rest = row["hash"].split("/")
    (folder,) = (work / "work" / prefix).glob(f"{rest}*")
    return folder


def star_record(references: Path, ref_id: str, annotation: str, gtf: Path) -> list[str]:
    """The <index>.sha256 lines RNA_INDEX writes, with index files sorted (glob order varies)."""
    name = f"{ref_id}.{annotation.removesuffix('.gz')}.star"
    lengths = [
        int(line.split("\t")[1])
        for line in (references / f"{ref_id}.chrom.sizes").read_text().splitlines()
    ]
    length = sum(lengths)
    parameters = {
        "sjdbOverhang": 100,
        "genomeSAindexNbases": min(14, math.floor(math.log2(length) / 2 - 1)),
        "genomeChrBinNbits": min(18, math.floor(math.log2(max(length / len(lengths), 100)))),
    }
    header = [
        f"{sha256(references / f'{ref_id}.fa.gz')}  {ref_id}.fa.gz",
        f"{sha256(gtf)}  {ref_id}.genes.gtf",
        f"{STAR_VERSION}  STAR",
    ] + [f"{value}  {key}" for key, value in parameters.items()]
    files = [path for path in (references / name).iterdir() if path.is_file()]
    return header + sorted(f"{path.stat().st_size}  {name}/{path.name}" for path in files)


def read_star_record(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    return lines[:STAR_HEADER] + sorted(lines[STAR_HEADER:])


def synthetic_products(references: Path, work: Path) -> None:
    """Reference products as a real run publishes them, from the first stub run's tasks.

    The STAR index and record are the ones the stub wrote, with the pinned STAR version for the
    stub's, so the reuse cases also exercise RNA_INDEX's record writer.
    """
    for ref_id in ("plant", "lifted"):
        sizes = task_dir(work, f"GET_REFERENCES:CHROM_SIZES ({ref_id})") / f"{ref_id}.chrom.sizes"
        shutil.copy2(sizes, references / sizes.name)
        (references / f"{ref_id}.bwa-mem2").mkdir()
        for suffix in INDEX_SUFFIXES:
            (references / f"{ref_id}.bwa-mem2/genome{suffix}").write_text(f"{ref_id}{suffix}\n")
        (references / f"{ref_id}.fa.gz.sha256").write_text(
            reference_record(references, ref_id, f"{ref_id}.fa.gz")
        )
        built = task_dir(work, f"GET_STAR_INDEX:RNA_INDEX ({ref_id})")
        name = f"{ref_id}.{ref_id}.gff3.star"
        shutil.copytree(built / name, references / name)
        record = (built / f"{name}.sha256").read_text()
        assert "STUB  STAR\n" in record, record
        (references / f"{name}.sha256").write_text(
            record.replace("STUB  STAR\n", f"{STAR_VERSION}  STAR\n")
        )
        gtf = task_dir(work, f"RNA_ANNOTATION ({ref_id})") / f"{ref_id}.genes.gtf"
        expected = star_record(references, ref_id, f"{ref_id}.gff3", gtf)
        assert read_star_record(references / f"{name}.sha256") == expected, ref_id


def check_star_rebuild(work: Path, sheet: Path, *, engine: str | None) -> None:
    """Rebuild a STAR index only when its record no longer matches its inputs and files."""
    parent = work / "star-rebuild"
    parent.mkdir()
    built = work / "references"
    record_name = "plant.plant.gff3.star.sha256"
    if engine:
        # The real run published indexes and records the reader must accept.
        gtf = task_dir(work, "RNA_ANNOTATION (plant)") / "plant.genes.gtf"
        expected = star_record(built, "plant", "plant.gff3", gtf)
        assert read_star_record(built / record_name) == expected

    def changed_fasta(references: Path) -> None:
        fasta = references / "plant.fa.gz"
        content = gzip.decompress(fasta.read_bytes())
        fasta.write_bytes(gzip.compress(content.replace(b"A", b"T", 1), mtime=0))

    def changed_annotation(references: Path) -> None:
        annotation = references / "plant.gff3"
        annotation.write_text(annotation.read_text().replace("\t10900\t", "\t10901\t", 1))

    def resized_index(references: Path) -> None:
        index = references / "plant.plant.gff3.star/SA"
        index.write_bytes(index.read_bytes() + b"x")  # e.g. replaced or partly copied

    def star_version(references: Path) -> None:
        record = references / record_name
        record.write_text(record.read_text().replace(f"{STAR_VERSION}  STAR", "2.7.10a  STAR"))

    def second_annotation(references: Path) -> None:
        shutil.copy2(references / "plant.gff3", references / "plant.v2.gff3")

    cases: list[tuple[str, Callable[[Path], object] | None, bool]] = [
        ("unchanged", None, False),
        ("changed-fasta", changed_fasta, True),
        ("changed-annotation", changed_annotation, True),
        ("missing-record", lambda r: (r / record_name).unlink(), True),
        ("resized-index", resized_index, True),
        (
            "extra-index-file",
            lambda r: (r / "plant.plant.gff3.star/extra").write_text("x"),
            True,
        ),
        ("star-version", star_version, True),
        # Another annotation for the same FASTA gets its own index beside the first.
        ("second-annotation", second_annotation, True),
    ]
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    for scenario, change, rebuilt in cases:
        case = parent / scenario
        case.mkdir()
        references = case / "references"
        shutil.copytree(built, references)
        if not engine:
            synthetic_products(references, work)
        shutil.copy2(work / "limits.config", case / "limits.config")
        if change:
            change(references)
        annotation = "plant.v2.gff3" if scenario == "second-annotation" else "plant.gff3"
        case_sheet = case / "samples.tsv"
        write_sheet(
            case_sheet,
            [
                {**row, "reference_annotation": annotation}
                if row["reference_fasta"] == "plant.fa.gz"
                else row
                for row in rows
            ],
        )
        before = {
            path: path.read_bytes() for path in sorted(references.rglob("*")) if path.is_file()
        }
        assert run_rna(case, case_sheet, engine=engine).returncode == 0, scenario
        names = [row["name"] for row in trace_rows(case)]
        # Only the altered plant index is rebuilt; the lifted index is always reused.
        assert ("GET_STAR_INDEX:RNA_INDEX (plant)" in names) == rebuilt, (
            scenario,
            names,
        )
        assert "GET_STAR_INDEX:RNA_INDEX (lifted)" not in names, (scenario, names)
        after = {
            path: path.read_bytes() for path in sorted(references.rglob("*")) if path.is_file()
        }
        if engine and rebuilt:
            # A rebuild publishes a fresh index and a record matching it.
            gtf = task_dir(case, "RNA_ANNOTATION (plant)") / "plant.genes.gtf"
            name = f"plant.{annotation}.star.sha256"
            expected = star_record(references, "plant", annotation, gtf)
            assert read_star_record(references / name) == expected, scenario
            if scenario == "second-annotation":
                # The first annotation's index is left as it was.
                assert all(after[path] == data for path, data in before.items()), scenario
        else:
            # Stub runs never publish into the references folder, and reuse changes nothing.
            assert after == before, scenario


def check_failures(work: Path, sheet: Path) -> None:
    """A rejected sheet stops before any task."""
    with sheet.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    lines = sheet.read_text().splitlines(keepends=True)
    cases: list[tuple[str, str, str]] = []
    for suffix, changed, diagnostic in (
        ("strand", [{**rows[0], "strandedness": "dUTP"}], "strandedness must be"),
        (
            "missing",
            [{**rows[0], "reference_annotation": "absent.gff3"}],
            "annotation missing",
        ),
        # Both rows use plant.fa.gz but name different organelles.
        ("organelles", [rows[0], {**rows[1], "organelle_contigs": "-"}], "disagree on"),
    ):
        candidate = work / f"{suffix}.tsv"
        write_sheet(candidate, changed)
        cases.append((suffix, str(candidate), diagnostic))
    extra = work / "extra-column.tsv"
    extra.write_text("".join(line.rstrip("\n") + "\tnote\n" for line in lines))
    cases.append(("extra-column", str(extra), "expected header"))
    for suffix, candidate, diagnostic in cases:
        assert run_rna(work, Path(candidate), suffix=suffix, expect_failure=True).returncode != 0
        assert trace_rows(work, suffix) == []
        log = (work / f"output-{suffix}.log").read_text()
        assert diagnostic in log, (suffix, log[-2000:])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    engines = parser.add_mutually_exclusive_group()
    engines.add_argument(
        "--docker",
        action="store_true",
        help="Use real tools in the pinned Docker images",
    )
    engines.add_argument(
        "--conda",
        action="store_true",
        help="Use real tools from the environments/*.yaml files",
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Retain fixtures, logs, and work directories",
    )
    args = parser.parse_args()
    engine = "docker" if args.docker else "conda" if args.conda else None
    runs = ROOT / "tests" / ".runs"
    runs.mkdir(exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f"rna-{engine or 'regression'}-", dir=runs))
    passed = False
    try:
        (work / "inputs").mkdir()
        with serve(work / "inputs", docker=engine == "docker") as base_url:
            sheet = create_fixture(work, base_url, engine=engine)
            assert run_rna(work, sheet, engine=engine).returncode == 0
            check_outputs(work, engine=engine)
            check_resume(work, sheet, engine=engine)
            check_star_rebuild(work, sheet, engine=engine)
            if not engine:
                check_failures(work, sheet)
        label = f"Real-tool ({engine}) RNA-seq" if engine else "RNA-seq Nextflow regression"
        print(f"{label} validation passed.", flush=True)
        passed = True
    finally:
        if args.keep or not passed:
            print(f"Artifacts retained: {work}", flush=True)
        else:
            shutil.rmtree(work)


if __name__ == "__main__":
    main()
