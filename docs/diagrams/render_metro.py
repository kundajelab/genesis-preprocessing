"""Deterministic, repository-native SVG maps. No network or rendering dependencies."""

import argparse
import hashlib
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COLORS = {
    "reads": "#087f8c", "reference": "#4169b1", "qc": "#9260b7",
    "peaks": "#ce6030", "tracks": "#438362", "spp": "#bd527d",
    "future": "#81909f",
}


Point = tuple[int, int]
Node = tuple[str, int, int, str, str, str]
Edge = tuple[str, list[Point], str, bool]


def route(points: list[Point]) -> str:
    result = f"M {points[0][0]} {points[0][1]}"
    for i in range(1, len(points) - 1):
        a, b, c = points[i - 1], points[i], points[i + 1]
        l1 = abs(b[0] - a[0]) + abs(b[1] - a[1])
        l2 = abs(c[0] - b[0]) + abs(c[1] - b[1])
        radius = min(14, l1 / 2, l2 / 2)
        p = (b[0] + (a[0] - b[0]) * radius / l1,
             b[1] + (a[1] - b[1]) * radius / l1)
        q = (b[0] + (c[0] - b[0]) * radius / l2,
             b[1] + (c[1] - b[1]) * radius / l2)
        result += f" L {p[0]:g} {p[1]:g} Q {b[0]} {b[1]} {q[0]:g} {q[1]:g}"
    return result + f" L {points[-1][0]} {points[-1][1]}"


def document(title: str, subtitle: str, nodes: list[Node], edges: list[Edge],
             notes: list[str], workflows: list[str], animated: bool = False) -> str:
    css = """
    svg{--bg:#fbfcfe;--text:#172c3d;--muted:#536779;--panel:#e2e8ef;background:var(--bg);font-family:Arial,Helvetica,sans-serif}
    .title{font-size:32px;font-weight:700;letter-spacing:1px;fill:var(--text)}
    .subtitle,.note{font-size:17px;fill:var(--muted)}
    .label{font-size:21px;font-weight:700;fill:var(--text)}
    .detail{font-size:16px;fill:var(--muted)}
    .station{fill:var(--bg);stroke-width:4}
    .line{fill:none;stroke-width:5;stroke-linecap:round;stroke-linejoin:round}
    .planned{stroke-dasharray:9 9}
    .flow{fill:none;stroke:white;stroke-width:2;stroke-dasharray:2 38;stroke-linecap:round;animation:travel 12s linear infinite;pointer-events:none}
    @keyframes travel{to{stroke-dashoffset:-400}}
    @media(prefers-reduced-motion:reduce){.flow{display:none}}
    @media(prefers-color-scheme:dark){svg{--bg:#14212d;--text:#eff6fc;--muted:#bdcbd6;--panel:#344957}.line,.station{filter:brightness(1.25)}}
    @media print{.flow{display:none}svg{--bg:#fff;--text:#172c3d;--muted:#536779;--panel:#e2e8ef}}
    """
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1480 870" role="img" aria-labelledby="title description">',
        f'<title id="title">{escape(title)}</title>',
        f'<desc id="description">{escape(subtitle + ". " + " ".join(notes))}</desc>',
        f"<style>{css}</style>",
        '<rect width="1480" height="870" fill="var(--bg)"/>',
        f'<text class="title" x="48" y="52">{escape(title)}</text>',
        f'<text class="subtitle" x="48" y="85">{escape(subtitle)}</text>',
    ]
    for i, (color, points, label, planned) in enumerate(edges):
        geometry = route(points)
        cls = "line planned" if planned else "line"
        parts.append(f'<path id="route-{i}" class="{cls}" stroke="{COLORS[color]}" d="{geometry}"><title>{escape(label)}</title></path>')
        if animated and not planned:
            parts.append(f'<path class="flow" style="animation-delay:-{i % 12}s" d="{geometry}"/>')
    for ident, x, y, label, detail, color in nodes:
        parts.extend([
            f'<g id="{ident}"><title>{escape(label + ": " + detail)}</title>',
            f'<circle class="station" cx="{x}" cy="{y}" r="9" stroke="{COLORS[color]}"/>',
            f'<text class="label" text-anchor="middle" x="{x}" y="{y - 40}">{escape(label)}</text>',
            f'<text class="detail" text-anchor="middle" x="{x}" y="{y - 17}">{escape(detail)}</text></g>',
        ])
    parts.append('<path d="M48 762 H1432" stroke="var(--panel)" stroke-width="2"/>')
    for i, note in enumerate(notes):
        parts.append(f'<text class="note" x="48" y="{798 + i * 25}">{escape(note)}</text>')
    fingerprints = "; ".join(
        f"{name} SHA256 {hashlib.sha256((ROOT / name).read_bytes()).hexdigest()}"
        for name in workflows
    )
    parts.append(f"<metadata>Genesis diagram schema 2; {fingerprints}; schematic animation.</metadata></svg>")
    return "\n".join(parts) + "\n"


def edge(color: str, points: list[Point], label: str, planned: bool = False) -> Edge:
    return color, points, label, planned


DAP_NODES = [
    ("samples", 140, 330, "Sample sheet", "SE / PE + controls", "reads"),
    ("validate", 370, 330, "Validate sheet", "Schema + controls", "reads"),
    ("DOWNLOAD", 600, 330, "Download", "Unmodified raw mates", "reads"),
    ("BWA_MEM2_ALIGN", 830, 330, "bwa-mem2", "Full-read alignment", "reads"),
    ("CALL_PEAKS", 1060, 330, "MACS3 peaks", "Assigned control", "peaks"),
    ("QUANTIFY", 1290, 330, "Quantify", "RPM + mean RPKM", "peaks"),
    ("FASTQC", 600, 180, "FastQC", "One report per mate", "qc"),
    ("QC", 1060, 180, "samtools QC", "Three metric families", "qc"),
    ("MULTIQC", 1290, 180, "MultiQC", "HTML + parsed data", "qc"),
    ("CHROM_SIZES", 140, 510, "Chromosome sizes", "Verified products", "reference"),
    ("METADATA", 600, 510, "Metadata", "Genome + read length", "reads"),
    ("spp", 830, 510, "SPP", "Separate 50-bp R1", "spp"),
    ("GET_REFERENCES", 140, 690, "Reference FASTA", "Hash-checked reuse", "reference"),
    ("BWA_MEM2_INDEX", 370, 690, "bwa-mem2 index", "Verified products", "reference"),
    ("TRACKS", 1060, 690, "deepTools tracks", "BigWig signal", "tracks"),
]
DAP_EDGES = [
    edge("reads", [(140,330),(600,330),(830,330)], "Validated libraries download and align"),
    edge("peaks", [(830,330),(1060,330),(1290,330)], "Treatment/control peaks and quantification"),
    edge("qc", [(600,330),(715,330),(715,180),(600,180)], "Raw mates independently feed FastQC"),
    edge("qc", [(600,180),(715,180),(715,230),(1405,230),(1405,180),(1290,180)], "FastQC reports feed MultiQC"),
    edge("qc", [(1060,180),(1290,180)], "Native samtools reports feed MultiQC"),
    edge("qc", [(830,330),(945,330),(945,180),(1060,180)], "Alignment QC feeds MultiQC"),
    edge("reads", [(600,330),(485,330),(485,510),(600,510)], "Read metadata"),
    edge("reads", [(600,510),(715,510),(715,330),(830,330)], "Metadata accompanies alignment"),
    edge("reference", [(140,690),(255,690),(255,330),(370,330)], "Reference identity is validated"),
    edge("reference", [(370,330),(255,330),(255,510),(140,510)], "Build or verify chromosome sizes"),
    edge("reference", [(370,330),(255,330),(255,690),(370,690)], "Build or verify index"),
    edge("reference", [(140,510),(600,510)], "Chromosome sizes feed metadata and index resources"),
    edge("reference", [(370,690),(715,690),(715,330),(830,330)], "Verified index feeds alignment"),
    edge("spp", [(830,330),(945,330),(945,510),(830,510)], "Separate shortened-read alignment feeds SPP"),
    edge("tracks", [(830,330),(945,330),(945,690),(1060,690)], "Full-read coverage"),
    edge("reads", [(830,330),(945,330),(945,390),(1405,390),(1405,330),(1290,330)], "Primary treatment BAM feeds quantification"),
    edge("reference", [(140,510),(255,510),(255,735),(1405,735),(1405,330),(1290,330)], "Chromosome sizes bound quantification"),
    edge("tracks", [(1060,690),(1405,690),(1405,330),(1290,330)], "Coverage and chromosome bounds support quantification"),
    edge("qc", [(600,510),(485,510),(485,105),(1405,105),(1405,180),(1290,180)], "Metadata joins run-level reporting"),
]
ATAC_NODES = [
    ("samples",140,330,"Sample sheet","Paired-end libraries","reads"),
    ("DOWNLOAD",370,330,"Download","Unmodified raw mates","reads"),
    ("ATAC_TRIM",600,330,"cutadapt","Nextera adapters","reads"),
    ("ATAC_ALIGN",830,330,"bwa-mem2","Align + mark duplicates","reads"),
    ("ATAC_FILTER",1060,330,"Filter","MAPQ, organelles, duplicates","reads"),
    ("ATAC_PEAKS",1290,330,"MACS3 peaks","BAMPE fragments","peaks"),
    ("FASTQC",370,180,"FastQC","One report per mate","qc"),
    ("MULTIQC_REPORT",1290,180,"MultiQC","HTML + parsed data","qc"),
    ("ATAC_QC",1290,510,"ATAC QC","FRiP, TSS, complexity","qc"),
    ("CHROM_SIZES",140,510,"Chromosome sizes","Verified products","reference"),
    ("GET_REFERENCES",140,690,"Reference FASTA","Hash-checked reuse","reference"),
    ("BWA_MEM2_INDEX",370,690,"bwa-mem2 index","Verified products","reference"),
    ("genes",600,690,"Gene starts","Optional BED6","reference"),
    ("ATAC_TRACKS",1060,690,"Tn5 cut sites","+4 / −5 shifted bigWig","tracks"),
]
ATAC_EDGES = [
    edge("reads",[(140,330),(1060,330)],"Libraries download, trim, align and filter"),
    edge("peaks",[(1060,330),(1290,330)],"Usable fragments feed BAMPE peaks"),
    edge("qc",[(370,330),(370,180)],"Raw mates independently feed FastQC"),
    edge("qc",[(370,180),(1290,180)],"FastQC reports feed MultiQC"),
    edge("qc",[(600,330),(715,330),(715,250),(1405,250),(1405,180),(1290,180)],"Trimming, alignment and filter reports feed MultiQC"),
    edge("tracks",[(1060,330),(1175,330),(1175,690),(1060,690)],"Usable reads give shifted cut sites"),
    edge("qc",[(1290,330),(1290,510)],"Peaks give FRiP"),
    edge("qc",[(1060,690),(1290,690),(1290,510)],"Cut sites give the TSS profile"),
    edge("qc",[(1290,510),(1405,510),(1405,180),(1290,180)],"ATAC QC table and TSS profile feed MultiQC"),
    edge("reference",[(140,690),(255,690),(255,510),(140,510)],"Build or verify chromosome sizes"),
    edge("reference",[(140,690),(370,690)],"Build or verify index"),
    edge("reference",[(370,690),(370,600),(830,600),(830,330)],"Verified index feeds alignment"),
    edge("reference",[(140,510),(945,510),(945,330),(1060,330)],"Chromosome sizes give nuclear contigs and genome size"),
    edge("reference",[(600,690),(715,690),(715,735),(1350,735),(1350,510),(1290,510)],"Gene starts define strand-aware TSS windows"),
]
RNA_NODES = [
    ("samples",140,330,"Sample sheet","SE / PE + strandedness","reads"),
    ("DOWNLOAD",370,330,"Download","Unmodified raw mates","reads"),
    ("RNA_ALIGN",830,330,"STAR","Annotated junctions","reads"),
    ("RNA_COUNTS",1060,330,"featureCounts","Unique reads, stranded","reads"),
    ("RNA_QC",1290,330,"RNA QC","Strand, rRNA, organelles","qc"),
    ("FASTQC",370,180,"FastQC","One report per mate","qc"),
    ("MULTIQC_REPORT",1290,180,"MultiQC","HTML + parsed data","qc"),
    ("CHROM_SIZES",140,510,"Chromosome sizes","Verified products","reference"),
    ("RNA_ANNOTATION",370,510,"Annotation","GTF + longest intron","reference"),
    ("GET_REFERENCES",140,690,"Reference FASTA","Hash-checked reuse","reference"),
    ("GET_STAR_INDEX",600,690,"STAR index","Hash-checked reuse","reference"),
    ("RNA_TRACKS",1060,690,"CPM track","Per-base bigWig","tracks"),
]
RNA_EDGES = [
    edge("reads",[(140,330),(830,330)],"Libraries download and align"),
    edge("reads",[(830,330),(1060,330)],"Uniquely aligned reads are counted per gene"),
    edge("qc",[(1060,330),(1290,330)],"Counts and assignment feed RNA QC"),
    edge("qc",[(370,330),(370,180)],"Raw mates independently feed FastQC"),
    edge("qc",[(370,180),(1290,180)],"FastQC reports feed MultiQC"),
    edge("qc",[(830,330),(945,330),(945,250),(1405,250),(1405,180),(1290,180)],"STAR and featureCounts logs feed MultiQC"),
    edge("qc",[(1290,330),(1290,180)],"RNA QC table feeds MultiQC"),
    edge("tracks",[(830,330),(875,330),(875,690),(1060,690)],"Aligned reads give a CPM track"),
    edge("reference",[(140,690),(255,690),(255,510),(140,510)],"Build or verify chromosome sizes"),
    edge("reference",[(140,510),(370,510)],"Annotation contigs must be in the FASTA"),
    edge("reference",[(140,690),(600,690)],"FASTA feeds the STAR index"),
    edge("reference",[(370,510),(485,510),(485,690),(600,690)],"Annotated junctions feed the STAR index"),
    edge("reference",[(600,690),(715,690),(715,330),(830,330)],"Verified index and intron limit feed alignment"),
    edge("reference",[(370,510),(945,510),(945,330),(1060,330)],"GTF defines the genes to count"),
    edge("reference",[(370,510),(370,600),(1350,600),(1350,330),(1290,330)],"rRNA genes and organelle contigs feed QC"),
]
SCATAC_NODES = [
    ("samples",140,330,"Sample sheet","Fragments + cell list","reads"),
    ("SCATAC_DOWNLOAD",370,330,"Download","Published files","reads"),
    ("SCATAC_FRAGMENTS",600,330,"Fragments","Listed cells, nuclear","reads"),
    ("groups",830,330,"Cell groups","One pseudobulk each","reads"),
    ("SCATAC_PEAKS",1060,330,"MACS3 peaks","Per group","peaks"),
    ("SCATAC_QC",1290,330,"scATAC QC","FRiP, TSS enrichment","qc"),
    ("MULTIQC_REPORT",1290,180,"MultiQC","HTML + parsed data","qc"),
    ("CHROM_SIZES",140,510,"Chromosome sizes","Verified products","reference"),
    ("GET_REFERENCES",140,690,"Reference FASTA","Hash-checked reuse","reference"),
    ("genes",600,690,"Gene starts","Optional BED6","reference"),
    ("SCATAC_TRACKS",1060,690,"Tn5 cut sites","Per-base bigWig","tracks"),
]
SCATAC_EDGES = [
    edge("reads",[(140,330),(830,330)],"Fragments and cells are checked, then split into groups"),
    edge("peaks",[(830,330),(1060,330)],"Each group's fragments give peaks"),
    edge("qc",[(1060,330),(1290,330)],"Peaks give FRiP"),
    edge("tracks",[(830,330),(945,330),(945,690),(1060,690)],"Each group's fragments give cut sites"),
    edge("qc",[(1060,690),(1290,690),(1290,330)],"Cut sites give the TSS profile"),
    edge("qc",[(1290,330),(1290,180)],"Group QC and TSS profiles feed MultiQC"),
    edge("qc",[(600,330),(715,330),(715,250),(1405,250),(1405,180),(1290,180)],"Library summary feeds MultiQC"),
    edge("reference",[(140,690),(255,690),(255,510),(140,510)],"Build or verify chromosome sizes"),
    edge("reference",[(140,510),(485,510),(485,330),(600,330)],"Every fragment must lie on a contig"),
    edge("reference",[(600,690),(715,690),(715,735),(1350,735),(1350,330),(1290,330)],"Gene starts define strand-aware TSS windows"),
]
SCATAC_READS_NODES = [
    ("samples",140,330,"Sample sheet","SRA runs + chemistry","reads"),
    ("SCATAC_SRA",370,330,"SRA download","Barcode read kept","reads"),
    ("SCATAC_CHROMAP",600,330,"chromap","10x fragments","reads"),
    ("SCATAC_CELLS",830,330,"Cells","Called or listed","reads"),
    ("SCATAC_FRAGMENTS",1060,330,"Fragments","Cells, nuclear","reads"),
    ("SCATAC_PEAKS",1290,330,"MACS3 peaks","Per group","peaks"),
    ("SCATAC_QC",1290,510,"scATAC QC","FRiP, TSS enrichment","qc"),
    ("MULTIQC_REPORT",1290,180,"MultiQC","HTML + parsed data","qc"),
    ("whitelist",370,510,"Barcode list","User-supplied","reference"),
    ("CHROM_SIZES",140,510,"Chromosome sizes","Verified products","reference"),
    ("GET_REFERENCES",140,690,"Reference FASTA","Hash-checked reuse","reference"),
    ("GET_CHROMAP_INDEX",370,690,"chromap index","Hash-checked reuse","reference"),
    ("genes",600,690,"Gene starts","Optional BED6","reference"),
    ("SCATAC_TRACKS",1060,690,"Tn5 cut sites","Per-base bigWig","tracks"),
]
SCATAC_READS_EDGES = [
    edge("reads",[(140,330),(1060,330)],"Runs are downloaded, mapped per library, then cells are called"),
    edge("peaks",[(1060,330),(1290,330)],"Each group's fragments give peaks"),
    edge("qc",[(1290,330),(1290,510)],"Peaks give FRiP"),
    edge("tracks",[(1060,330),(1060,690)],"Each group's fragments give cut sites"),
    edge("qc",[(1060,690),(1290,690),(1290,510)],"Cut sites give the TSS profile"),
    edge("qc",[(1290,510),(1405,510),(1405,180),(1290,180)],"Group QC and TSS profiles feed MultiQC"),
    edge("reference",[(370,510),(485,510),(485,330),(600,330)],"Barcodes are corrected against the list"),
    edge("reference",[(140,690),(255,690),(255,510),(140,510)],"Build or verify chromosome sizes"),
    edge("reference",[(140,690),(370,690)],"The FASTA gives the chromap index"),
    edge("reference",[(370,690),(540,690),(540,330),(600,330)],"Verified index feeds chromap"),
    edge("reference",[(140,510),(255,510),(255,600),(945,600),(945,330),(1060,330)],"Every fragment must lie on a contig"),
    edge("reference",[(600,690),(715,690),(715,780),(1350,780),(1350,510),(1290,510)],"Gene starts define strand-aware TSS windows"),
]
TSS_NODES = [
    ("samples",140,330,"Sample sheet","Method + trimming","reads"),
    ("DOWNLOAD",370,330,"Download","Unmodified raw reads","reads"),
    ("TSS_TRIM",600,330,"cutadapt","5' barcode, 3' adapter","reads"),
    ("TSS_ALIGN",830,330,"STAR","5' end kept","reads"),
    ("TSS_CTSS",1060,330,"TSS tags","5' ends, cap G corrected","reads"),
    ("TSS_CLUSTERS",1290,330,"paraclu","FANTOM clusters","peaks"),
    ("TSS_QC",1290,510,"TSS QC","Near annotated starts","qc"),
    ("FASTQC",370,180,"FastQC","Raw reads","qc"),
    ("MULTIQC_REPORT",1290,180,"MultiQC","HTML + parsed data","qc"),
    ("CHROM_SIZES",140,510,"Chromosome sizes","Verified products","reference"),
    ("RNA_ANNOTATION",370,510,"Annotation","GTF + longest intron","reference"),
    ("GET_REFERENCES",140,690,"Reference FASTA","Hash-checked reuse","reference"),
    ("GET_STAR_INDEX",600,690,"STAR index","Hash-checked reuse","reference"),
]
TSS_EDGES = [
    edge("reads",[(140,330),(1060,330)],"Reads are trimmed, aligned and reduced to their 5' ends"),
    edge("peaks",[(1060,330),(1290,330)],"Tags per strand give clusters"),
    edge("qc",[(1290,330),(1290,510)],"Clusters feed QC"),
    edge("qc",[(1060,330),(1060,600),(1175,600),(1175,510),(1290,510)],"Tags feed QC"),
    edge("qc",[(370,330),(370,180)],"Raw reads feed FastQC"),
    edge("qc",[(370,180),(1290,180)],"FastQC reports feed MultiQC"),
    edge("qc",[(1290,510),(1405,510),(1405,180),(1290,180)],"TSS QC, cutadapt and STAR logs feed MultiQC"),
    edge("reference",[(140,690),(255,690),(255,510),(140,510)],"Build or verify chromosome sizes"),
    edge("reference",[(140,510),(370,510)],"Annotation contigs must be in the FASTA"),
    edge("reference",[(370,510),(485,510),(485,690),(600,690)],"Annotated junctions feed the STAR index"),
    edge("reference",[(600,690),(715,690),(715,420),(830,420),(830,330)],"Verified index and intron limit feed alignment"),
    edge("reference",[(370,510),(370,780),(1350,780),(1350,510),(1290,510)],"Annotated transcript starts feed QC"),
]
CONTROL_NODES = [
    ("inputs",140,330,"Datasets","Exact references + policy","reads"),
    ("controller",370,330,"Controller","Queue + durable status","reads"),
    ("local",830,180,"Local Linux","Default Docker execution","reads"),
    ("scheduler",830,330,"Cluster scheduler","Slurm / PBS / LSF / SGE","reference"),
    ("ssh",830,510,"SSH workers","Existing verified access","tracks"),
    ("outputs",1290,330,"Results","QC + provenance","peaks"),
    ("monitor",1060,690,"Monitor","Read-only status page","qc"),
    ("review",370,690,"Human review","Investigate + retry","spp"),
    ("retry",600,180,"Bounded retries","Confirmed transient errors","qc"),
]
CONTROL_EDGES = [
    edge("reads",[(140,330),(370,330),(485,330),(485,180),(830,180)],"Default local execution"),
    edge("reference",[(370,330),(830,330)],"Native scheduler submits tasks to managed compute nodes"),
    edge("tracks",[(370,330),(485,330),(485,510),(830,510)],"Independent Linux workers run isolated datasets over existing SSH"),
    edge("reads",[(830,180),(1175,180),(1175,330),(1290,330)],"Local results"),
    edge("reference",[(830,330),(1290,330)],"Cluster results"),
    edge("tracks",[(830,510),(1175,510),(1175,330),(1290,330)],"Worker results"),
    edge("qc",[(1290,330),(1405,330),(1405,690),(1060,690),(370,690)],"Status and unresolved issues flow to human review"),
    edge("spp",[(370,690),(255,690),(255,330),(370,330)],"Explicit resolution returns confirmed failed datasets to the queue"),
]
ROAD_NODES = [
    ("manifest",140,330,"Manifest","Biological library IDs","reads"),
    ("acquisition",370,330,"Acquisition + QC","Original reads retained","reads"),
    ("dap",600,180,"DAP-seq","Implemented","reads"),
    ("dap_outputs",1060,180,"Peaks + quantification","Per treatment","reads"),
    ("bulk_atac",600,330,"Bulk ATAC","Paired-end libraries","tracks"),
    ("bulk_outputs",1060,330,"Fragments + ATAC QC","Separate biological replicates","tracks"),
    ("single_cell",600,510,"scATAC","Fragments or 10x reads","tracks"),
    ("cells",1060,510,"Cell groups","From the cell list","tracks"),
    ("pseudobulk",1290,510,"Pseudobulks","Group × library","tracks"),
    ("multiome",140,690,"Multiome identity","Tested linking contract","future"),
    ("rna",600,690,"Bulk RNA-seq","STAR + featureCounts","reads"),
    ("tss",830,780,"TSS mapping","CAGE, TSS-seq, csRNA","reads"),
    ("reference_registry",140,180,"Reference registry","Tested schema","reference"),
    ("models",1290,690,"Model inputs","Model-specific export","future"),
]
ROAD_EDGES = [
    edge("reads",[(140,330),(370,330),(485,330),(485,180),(600,180),(1060,180)],"Implemented DAP route"),
    edge("tracks",[(370,330),(600,330),(1060,330)],"Paired-end libraries bulk ATAC route"),
    edge("tracks",[(370,330),(485,330),(485,510),(600,510),(1060,510),(1290,510)],"scATAC from published fragments to group pseudobulks"),
    edge("reads",[(370,330),(485,330),(485,690),(600,690)],"Bulk RNA-seq route"),
    edge("reads",[(485,690),(485,780),(830,780)],"TSS mapping route"),
    edge("reference",[(140,180),(255,180),(255,330),(370,330)],"Versioned registry contract",True),
    edge("future",[(140,690),(600,690),(945,690),(945,510),(1060,510)],"Modality-specific QC with shared annotation provenance",True),
    edge("future",[(1060,330),(1175,330),(1175,690),(1290,690)],"Proposed bulk model export",True),
    edge("future",[(1290,510),(1405,510),(1405,690),(1290,690)],"Replicate-specific pseudobulk export",True),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs/images")
    out = parser.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=True)
    for assay, title, subtitle, nodes, edges, notes, workflows in [
        ("", "GENESIS / DAP-seq", "Implemented workflow · raw-read QC · explicit controls and reference identity",
         DAP_NODES, DAP_EDGES,
         ["QC is diagnostic. FastQC warnings do not reject a dataset; SPP is reported separately.",
          "Controls receive QC and tracks; treatments receive peaks and quantification. FASTQs and BAMs stay intermediate.",
          "Routes summarize data dependencies, not runtime. Crossings without a station do not join channels."],
         ["workflows/dap-seq.nf"]),
        ("_atac", "GENESIS / BULK ATAC", "Paired-end libraries · ENCODE-style filtering · separate biological replicates",
         ATAC_NODES, ATAC_EDGES,
         ["Adapters are trimmed first. Usable pairs need both mates at MAPQ ≥ 30 on nuclear contigs; duplicates are removed.",
          "BAMPE peaks use unshifted fragments. Cut sites use +4 / −5 shifts; the TSS profile uses strand-aware windows.",
          "QC is diagnostic, with no plant thresholds. Libraries and biological replicates are never pooled."],
         ["workflows/atac-seq.nf"]),
        ("_rna", "GENESIS / BULK RNA-seq", "SE / PE libraries · STAR with plant intron limits · unique-read gene counts",
         RNA_NODES, RNA_EDGES,
         ["The intron limit is the longest annotated intron; STAR's default suits mammals. Up to 20 loci per read.",
          "featureCounts counts uniquely aligned reads by the sheet's strandedness; rRNA and organelle genes stay in the counts.",
          "QC is diagnostic: the forward fraction checks the sheet's strandedness. Replicates are never pooled."],
         ["workflows/rna-seq.nf"]),
        ("_tss", "GENESIS / TSS MAPPING", "Single-read CAGE, TSS-seq and csRNA-seq · per-base TSS tags · paraclu clusters",
         TSS_NODES, TSS_EDGES,
         ["The 5' end of read 1 is the TSS: only fixed barcode bases are trimmed there; the 3' end is adapter- and quality-trimmed.",
          "Cap-trapping CAGE adds a G the genome lacks; such reads move 1 base downstream, as in FANTOM. Unique alignments only.",
          "paraclu clusters each strand; FANTOM's paraclu-cut keeps clusters ≤ 200 bases with a 2-fold density rise. No plant thresholds."],
         ["workflows/tss-seq.nf"]),
        ("_scatac", "GENESIS / scATAC PSEUDOBULKS", "Published fragments · one pseudobulk per cell group · separate libraries",
         SCATAC_NODES, SCATAC_EDGES,
         ["Starts from published fragments and a cell list. Organelle fragments and unlisted barcodes are counted, then dropped.",
          "Each group gets MACS3 peaks, a per-base cut-site bigWig (first and last fragment base), FRiP and TSS enrichment.",
          "QC is diagnostic, with no plant thresholds. Groups stay within one library; libraries are never pooled."],
         ["workflows/scatac-seq.nf"]),
        ("_scatac_reads", "GENESIS / scATAC FROM 10x READS", "NCBI SRA runs with the barcode read · chromap fragments · shared pseudobulk steps",
         SCATAC_READS_NODES, SCATAC_READS_EDGES,
         ["Runs come from NCBI SRA, which keeps the 10x barcode read; the barcode and genomic reads are told apart by length.",
          "chromap writes 10x fragments (Tn5 offsets applied), corrects barcodes against a user-supplied list and removes duplicates per cell.",
          "Without an authors' cell list, barcodes with enough nuclear fragments form one group. No plant TSS threshold is applied."],
         ["workflows/scatac-reads.nf"]),
    ]:
        for animated, suffix in [(False, ""), (True, "_animated")]:
            (out / f"genesis{assay}_metro_map{suffix}.svg").write_text(
                document(title, subtitle, nodes, edges, notes, workflows, animated))
    (out / "genesis_architecture_metro_map.svg").write_text(document(
        "GENESIS / ARCHITECTURE", "Available assay routes and future contracts · biological replicates remain separate",
        ROAD_NODES, ROAD_EDGES,
        ["Solid: available assay processing. Dashed: contracts or planned integration.",
         "RNA and ATAC retain separate QC and inclusion decisions. No implicit biological replicate pooling.",
         "Plant biological thresholds remain UNSPECIFIED. scATAC starts from published fragments or 10x reads on NCBI SRA; model export is not a production workflow."],
        ["workflows/dap-seq.nf", "workflows/atac-seq.nf", "workflows/rna-seq.nf", "workflows/scatac-seq.nf",
         "workflows/scatac-reads.nf", "workflows/tss-seq.nf"]))
    for animated, suffix in [(False, ""), (True, "_animated")]:
        (out / f"genesis_execution{suffix}.svg").write_text(document(
            "GENESIS / RUN AND REVIEW", "Start on one Linux server. Add managed clusters or independent workers when needed.",
            CONTROL_NODES, CONTROL_EDGES,
            ["Authentication stays under site control. No keys, permissions or cluster services are changed by Genesis.",
             "Healthy datasets continue. Failed datasets enter review; unreachable workers are never blindly resubmitted.",
             "Cluster profiles require site validation. Monitoring reports problems; people resolve data and scientific questions."],
            ["conf/execution.config"], animated))
    print("Rendered 13 Genesis SVG maps")


if __name__ == "__main__":
    main()
