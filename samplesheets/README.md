# Sample sheets

One folder per species, one sheet per study and modality:
`samplesheets/<Species>/<modality>-<study>.tsv`. Each sheet uses its workflow's exact header
(see the [DAP-seq](../docs/usage/dapseq-reference.md), [bulk ATAC](../docs/usage/bulk-atac.md),
[scATAC](../docs/usage/scatac-seq.md) and [RNA-seq](../docs/usage/rna-seq.md) guides). Read URLs
point to ENA over HTTPS, or to GSA (China) for CRA studies. Reference files are named in the sheet and looked up in the
folder given by `--references`; they are not downloaded automatically.

| Sheet | Workflow | Libraries | Source |
| --- | --- | ---: | --- |
| `Arabidopsis_thaliana/dap-seq-GSE60141.tsv` | dap-seq | 936 (SE) | O'Malley 2016 |
| `Arabidopsis_thaliana/dap-seq-PRJNA1177481.tsv` | dap-seq | 800 (PE) | multiDAP |
| `Arabidopsis_thaliana/atac-seq-curated.tsv` | atac-seq | 439 (PE), 62 studies | Genesis "Bulk ATAC Fable Curated" tab: wild-type Col-0, untreated (below) |
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

All Sorghum sheets use one reference: Ensembl Plants release 63 `Sorghum_bicolor_NCBIv3`
(BTx623), with the chloroplast (RefSeq NC_008602, renamed `Pt`) and mitochondrion (NC_008360,
renamed `Mt`) appended, because the Ensembl FASTA has neither. Without them, organelle reads
would align to nuclear copies of organelle DNA and the organelle QC would read zero. The DAP
plants are RTx430, so their reads are aligned to the BTx623 assembly.

Build the references folder with:

```bash
ens=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63
efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id='
curl -fsSLO "$ens/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz"
curl -fsSLO "$ens/gff3/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.63.gff3.gz"
curl -fsSL -o NC_008602.fa "${efetch}NC_008602"
curl -fsSL -o NC_008360.fa "${efetch}NC_008360"
{
    gzip -dc Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz
    sed '1s/^>.*/>Pt NC_008602 Sorghum bicolor chloroplast/' NC_008602.fa
    sed '1s/^>.*/>Mt NC_008360 Sorghum bicolor mitochondrion/' NC_008360.fa
} | gzip -n > Sorghum_bicolor_NCBIv3_PtMt.fa.gz
gzip -dc Sorghum_bicolor.Sorghum_bicolor_NCBIv3.63.gff3.gz |
    awk -v OFS='\t' '$3 == "gene" {print $1, $4 - 1, $5, ".", 0, $7}' > Sorghum_bicolor_NCBIv3.63.genes.bed
```

| File | SHA-256 |
| --- | --- |
| `Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz` (download) | `ad40aa753793907fee06b3fe18a57d97f2f72798f4bcef15417a87fe36df73fb` |
| `NC_008602.fa` (download) | `8c99da486f7f773c2978aafc6eb2ed4b71037d8ae631b4209810f377dcefffcd` |
| `NC_008360.fa` (download) | `35432359b7fe450b8c4df7fc8ddd62a440c4e5e24248939639a6597ae152a87d` |
| `Sorghum_bicolor_NCBIv3_PtMt.fa.gz`, decompressed | `1376e335a646121b519f659f8fd52c9cce41385b985bb80728356217708bbdd9` |
| `Sorghum_bicolor.Sorghum_bicolor_NCBIv3.63.gff3.gz` | `f178e26b4fa3f227c365f585d27ccbd1bd95e39e8023115e38f17e4a6e7d7339` |
| `Sorghum_bicolor_NCBIv3.63.genes.bed` | `3de15e99825dafa4567022a98e4760df05ece066139c7a9ea051cc4bb6db18ee` |

Check the FASTA with `gzip -dc Sorghum_bicolor_NCBIv3_PtMt.fa.gz | sha256sum`: the compressed
file's hash depends on the gzip version. The result has 869 contigs (chromosomes `1`–`10`, 857
scaffolds, `Pt` 140,754 bp and `Mt` 468,628 bp; 709,344,700 bp in all), every GFF3 sequence is
in it, and the BED has 34,118 genes.

## Not yet available

A. thaliana scATAC has no sheet. GSE155304, which the workflow was validated on, has its
fragments only inside GEO tar archives, on Ensembl contig names (`1`–`5`, not Phytozome's
`Chr1`), and its cell clusters only in the authors' R object. The studies in the "scATAC Fable
Curated" tab list raw single-cell reads, not fragment files.

Sorghum scATAC has no sheet. The scATAC workflow starts from published fragment files, and
neither Sorghum study in the Genesis "scATAC Fable Curated" tab has them:

- Mendieta 2024 (PRJNA1063172, BTx623 leaf): raw FASTQs from a plate-based design only.
- Swift 2023 (GSE248919, seedling multiome): R objects only.
