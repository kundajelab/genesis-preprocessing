# Test sample sheets

One sheet per workflow, of real public libraries. They check that a change or a new site runs
end to end, and give enough measured tasks, across each source's range of sizes and library
types, to fit task limits and retries. Their rows are copied unchanged from the curated sheets
in the [Stanford-Genesis Hugging Face
bucket](https://huggingface.co/buckets/Stanford-Genesis/genesis-curation)
(`curated-sample-sheets/<Species>/<modality>-curated.tsv`), which list every public library
found for each workflow and the reasons others were left out. Sorghum comes first; where it has
few libraries, A. thaliana fills the rest. Read URLs point to ENA over HTTPS, or to GSA (China)
for CRA studies; the scATAC sheet lists SRA runs.

| Sheet | Libraries | FASTQ | What it covers |
| --- | ---: | ---: | --- |
| `dap-seq-test.tsv` | 24 | 70 GB | 8 Sorghum and 16 A. thaliana libraries from 10 studies: 12 treatments, each with its control; SE and PE, 0.05-9.7 GB per library |
| `atac-seq-test.tsv` | 24 | 123 GB | 6 Sorghum (leaf, 3 studies) and 18 A. thaliana (one per study, several tissues); spans size, TSS enrichment and organelle fraction; 5 libraries download from GSA; 0.05-38 GB |
| `rna-seq-test.tsv` | 20 | 208 GB | Sorghum, 16 projects; SE and PE; `reverse`, `forward` and `unstranded`; 0.21-72 GB |
| `tss-seq-test.tsv` | 24 | 36 GB | A. thaliana (Sorghum has no curated TSS libraries), all SE: 9 CAGE and 5 each of nAnT-iCAGE, TSS-seq and csRNA-seq; 0.39-7.8 GB |
| `scatac-reads-test.tsv` | 26 | 804 GB | Every curated library: 8 Sorghum (6 multiome ATAC, 2 sci-ATAC) and 18 A. thaliana (ATAC v1 and multiome ATAC); 1-9 runs per library, 1-4 files per run; 7-139 GB, up to 1.8 billion reads |

The scatac-seq workflow starts from fragment files, which no curated sheet provides, so it has
no test sheet.

Each sheet uses its
workflow's exact header (see the [DAP-seq](../docs/usage/dapseq-reference.md),
[bulk ATAC](../docs/usage/bulk-atac.md), [RNA-seq](../docs/usage/rna-seq.md),
[TSS mapping](../docs/usage/tss-seq.md) and [scATAC](../docs/usage/scatac-seq.md) guides).
Reference files are named in the sheet and looked up in the folder given by
`-R`/`--references`; they are not downloaded automatically. On Perlmutter, point it at the
project's shared references folder.

```bash
pixi run submit-pipeline --workflow rna-seq -R /path/to/references samplesheets/rna-seq-test.tsv
```

## Arabidopsis thaliana

The A. thaliana rows use Phytozome genome 447: `Athaliana_447_TAIR10.fa.gz` with the
Araport11 gene models `Athaliana_447_Araport11.gene_exons.gff3.gz`. TAIR10 already contains the
chloroplast (`ChrC`) and mitochondrion (`ChrM`), so nothing is appended, and
`organelle_contigs` is `ChrC,ChrM`. Phytozome downloads need a JGI login (ORCID works): get both
files from the [Phytozome portal](https://phytozome-next.jgi.doe.gov/info/Athaliana_Araport11),
check them against JGI's MD5s, and make the genes BED:

```bash
md5sum -c - <<'MD5'
c8b69730af117a236f0f2c2c51b5abe1  Athaliana_447_TAIR10.fa.gz
350389b5c9892cda4c72c2f351a6d2c8  Athaliana_447_Araport11.gene_exons.gff3.gz
MD5
gzip -dc Athaliana_447_Araport11.gene_exons.gff3.gz |
    awk -v OFS='\t' '$3 == "gene" {print $1, $4 - 1, $5, ".", 0, $7}' > Athaliana_447_Araport11.genes.bed
```

| File | SHA-256 |
| --- | --- |
| `Athaliana_447_TAIR10.fa.gz` | `984ed0c2d2c03eec22a12165e582060bac80c2d41901c4b9e3a5ba8b4964af9c` |
| `Athaliana_447_Araport11.gene_exons.gff3.gz` | `235769442d71b61ba869ff510529911045e0942ea88cc2f77bf110a14f317e70` |
| `Athaliana_447_Araport11.genes.bed` | `dd4332597f1d3f3db3431da5018944c24cb4c801319cd5d8280082c3ac82e263` |

The FASTA has 7 contigs (`Chr1`–`Chr5`, `ChrC` 154,478 bp, `ChrM` 366,924 bp); the BED has
27,655 genes, 88 on `ChrC` and 122 on `ChrM`. Phytozome's Araport11 gene_exons file has
protein-coding genes only and no rRNA features, so the RNA-seq `rRNA` metric reads `NA`.

## Sorghum reference

The Sorghum rows use Phytozome genome 730 (BTx623 v5.0 assembly, v5.1 genes), the latest
Phytozome release. All Sorghum libraries are aligned to this BTx623 reference, including those
from other cultivars; there is no cultivar-specific assembly, so reads from regions that differ
from BTx623 may map less well. The organelles are handled as agreed with Alex Borowsky and
Gabriel Ferreras Garrucho:

- **Mitochondrion:** 47 of v5's 55 unplaced scaffolds (2.5 Mb, 291 genes) are pieces of the
  mitochondrial genome, about 5 copies' worth: at least 90% of each matches RefSeq NC_008360 at
  at least 99% identity. They are removed, and the complete RefSeq mitochondrion is appended as
  `Mt`, so mitochondrial reads have one place to map. Kept with them, every such read would be a
  multi-mapper (MAPQ 0) and be filtered out.
- **Chloroplast:** v5 has none (one 580 bp nuclear insertion only); RefSeq NC_008602 is
  appended as `Pt`.
- **Organelle genes:** the RefSeq gene models are added, so RNA-seq counts chloroplast and
  mitochondrial genes. Trans-spliced genes (mitochondrial nad1, nad2, nad5; chloroplast
  rps12) become one transcript per piece under one gene, so no annotated intron joins pieces
  that are far apart or on opposite strands. The longest annotated intron stays 43 kb.
- The other 8 scaffolds are kept.

[`Sorghum_bicolor/reference/`](Sorghum_bicolor/reference/) builds it. Phytozome needs a JGI
login (ORCID works): download `Sbicolor_730_v5.0.fa.gz` and
`Sbicolor_730_v5.1.gene_exons.gff3.gz` from the
[Phytozome portal](https://phytozome-next.jgi.doe.gov/info/Sbicolor_v5_1) into an empty folder,
then run `samplesheets/Sorghum_bicolor/reference/build.sh` there. It checks JGI's MD5s,
downloads the RefSeq organelle FASTA and GFF3, and writes the three reference files.
`mito-scaffolds.txt` lists the removed scaffolds.

| File | SHA-256 |
| --- | --- |
| `Sbicolor_730_v5.0_PtMt.fa.gz`, decompressed | `73fa0491ad14c1747d1774f8d5fb6eb8cf90188e2ef576f9ffa71bb9efde130b` |
| `Sbicolor_730_v5.1_PtMt.gene_exons.gff3.gz`, decompressed | `1c03fb965fd5687d789f9ea36ce10b8e0dc3b8bc8675527fa806c5658e4ae93b` |
| `Sbicolor_730_v5.1_PtMt.genes.bed` | `bd27df170c5583d97bab10c8fb8b5cd638325a0899961f0b9dd8d33bb3f485b0` |

Check the compressed files with `gzip -dc FILE | sha256sum`: their own hash depends on the gzip
version. The FASTA has 20 contigs (`Chr01`–`Chr10`, 8 scaffolds, `Mt` 468,628 bp and `Pt`
140,754 bp; 717,997,770 bp in all), every GFF3 sequence is in it, and the annotation has
32,062 genes, 193 of them on `Mt` and `Pt`. Gene IDs carry a suffix, as in
`Sobic.001G000100.v5.1`. The BED lists each piece of a trans-spliced gene separately.

## scATAC barcode lists

The 10x ATAC barcode list (`737K-cratac-v1.txt`) is downloaded into the references folder when
it is missing. 10x does not allow redistributing it, so it is not in this repository. The
multiome ATAC libraries (6 Sorghum and 8 A. thaliana) need the multiome list
(`737K-arc-v1.txt`), which only ships inside Cell Ranger ARC, so their rows use `-`: barcodes
are not corrected, and reads with a barcode error are dropped. With Cell Ranger ARC's list in
the references folder, set `barcode_whitelist` to its name. The 2 sci-ATAC libraries (Sorghum)
have no barcode list: their combinatorial barcodes are used whole, uncorrected, so their rows
must use `-`.
