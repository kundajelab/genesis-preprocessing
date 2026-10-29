# Sample sheets

One folder per species, one sheet per study and modality:
`samplesheets/<Species>/<modality>-<study>.tsv`. Each sheet uses its workflow's exact header
(see the [DAP-seq](../docs/usage/dapseq-reference.md), [TSS mapping](../docs/usage/tss-seq.md), [bulk ATAC](../docs/usage/bulk-atac.md),
[scATAC](../docs/usage/scatac-seq.md) and [RNA-seq](../docs/usage/rna-seq.md) guides). Read URLs
point to ENA over HTTPS, or to GSA (China) for CRA studies. Reference files are named in the sheet and looked up in the
folder given by `--references`; they are not downloaded automatically.

| Sheet | Workflow | Libraries | Source |
| --- | --- | ---: | --- |
| `Arabidopsis_thaliana/dap-seq-GSE60141.tsv` | dap-seq | 936 (SE) | O'Malley 2016 |
| `Arabidopsis_thaliana/dap-seq-PRJNA1177481.tsv` | dap-seq | 800 (PE) | multiDAP |
| `Arabidopsis_thaliana/atac-seq-curated.tsv` | atac-seq | 439 (PE), 62 studies | Genesis "Bulk ATAC Fable Curated" tab: wild-type Col-0, untreated (below) |
| `Arabidopsis_thaliana/tss-seq-GSE136356.tsv` | tss-seq | 9 (SE CAGE) | GSE136356 (Thieffry 2020), Col-0 seedlings: wild type, hen2-4 and rrp4-2, 3 replicates each; reads start with a 3-base barcode and `CAGCAG` (`trim_5p` 9) |
| `Arabidopsis_thaliana/rna-seq-GSE206323.tsv` | rna-seq | 16 (PE) | GSE206323, Col-0 inflorescence stem, the 8 wild-type samples, 2 runs each; TruSeq stranded mRNA (`reverse`) |
| `Arabidopsis_lyrata/dap-seq-PRJNA1177479.tsv` | dap-seq | 405 (PE) | multiDAP |
| `Sorghum_bicolor/dap-seq-PRJNA1177471.tsv` | dap-seq | 142 (PE) | multiDAP, RTx430 plants; Genesis DAP-seq set 24 |
| `Sorghum_bicolor/dap-seq-PRJNA1177471-test.tsv` | dap-seq | 6 (PE) | Subset of the sheet above, for a quick test |
| `Sorghum_bicolor/atac-seq-curated.tsv` | atac-seq | 24 (PE) | Genesis "Bulk ATAC Fable Curated" tab: CRA002789, SRP110225, SRP188687, SRP308654, SRP674033 |
| `Sorghum_bicolor/rna-seq-GSE184725.tsv` | rna-seq | 6 (SE) | GSE184725, BTx623 shoot, control and sulfur-deficient, 3 each; NEBNext Ultra II Directional (`reverse`) |

Sample IDs in the curated Sorghum ATAC sheet are `<study>_<tissue>_<replicate>`. The Genesis
metadata sheet lists 25 Sorghum ATAC runs; SRR8742452 (SRP188687) is left out because it is the
study's genomic input control, not an ATAC library.

## Arabidopsis thaliana

The A. thaliana sheets use Phytozome genome 447: `Athaliana_447_TAIR10.fa.gz` with the
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

The curated ATAC sheet takes the tab's A. thaliana runs that are paired-end, with a wild-type
or Col-0 genotype and no treatment (or a control or mock), one row per run, named
`<study>_<tissue>_<run>`. ERP138859 (Srikant 2022) is left out: its "wildtype" runs are about
20 natural accessions, not Col-0, mixed with genomic DNA controls.

Runs whose reads are not public as a pair are left out (63): CRA015871 (26) and CRA015872 (20)
and one CRA024818 run have no files on the GSA download server; SRP606657 (12) and SRP548013
(1) have no FASTQ in ENA; and SRP098594 (3) has a single unsplit FASTQ per run.

The RNA-seq sheet lists each run separately because the workflow takes one pair of FASTQs per
row; sum the two runs of a sample for sample-level counts.

## Sorghum reference

All Sorghum sheets use Phytozome genome 730 (BTx623 v5.0 assembly, v5.1 genes), the latest
Phytozome release, with the organelles handled as agreed with Alex Borowsky and Gabriel
Ferreras Garrucho:

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
- The other 8 scaffolds are kept. The DAP plants are RTx430, aligned to the BTx623 assembly.

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

## Not yet available

A. thaliana scATAC has no sheet. GSE155304, which the workflow was validated on, has its
fragments only inside GEO tar archives, on Ensembl contig names (`1`–`5`, not Phytozome's
`Chr1`), and its cell clusters only in the authors' R object. The studies in the "scATAC Fable
Curated" tab list raw single-cell reads, not fragment files.

Sorghum scATAC has no sheet. The scATAC workflow starts from published fragment files, and
neither Sorghum study in the Genesis "scATAC Fable Curated" tab has them:

- Mendieta 2024 (PRJNA1063172, BTx623 leaf): raw FASTQs from a plate-based design only.
- Swift 2023 (GSE248919, seedling multiome): R objects only.
