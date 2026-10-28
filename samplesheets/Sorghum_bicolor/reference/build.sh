#!/usr/bin/env bash
# Build the Genesis Sorghum reference in the current folder: Phytozome 730 (v5.0 assembly,
# v5.1 genes) without its 47 mitochondrial scaffolds (mito-scaffolds.txt), plus the RefSeq
# mitochondrion (Mt) and chloroplast (Pt) with their gene models. Phytozome needs a JGI login:
# download Sbicolor_730_v5.0.fa.gz and Sbicolor_730_v5.1.gene_exons.gff3.gz here first.
set -euo pipefail
here=$(cd -- "$(dirname -- "$0")" && pwd)
efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id='
gff='https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id='
md5sum -c - << 'MD5'
2808545a9ed3c471a9e3fe98b3f5cf9a  Sbicolor_730_v5.0.fa.gz
410c6b14c7094f0870714d2394ff02c8  Sbicolor_730_v5.1.gene_exons.gff3.gz
MD5
for id in NC_008360.1 NC_008602.1; do
    curl -fsSL -o "$id.fa" "$efetch$id"
    curl -fsSL -o "$id.gff3" "$gff$id"
done
{
    gzip -dc Sbicolor_730_v5.0.fa.gz |
        awk 'NR == FNR {drop[$1] = 1; next} /^>/ {keep = !(substr($1, 2) in drop)} keep' \
            "$here/mito-scaffolds.txt" -
    sed '1s/^>.*/>Mt NC_008360.1 Sorghum bicolor mitochondrion, complete genome/' NC_008360.1.fa
    sed '1s/^>.*/>Pt NC_008602.1 Sorghum bicolor chloroplast, complete genome/' NC_008602.1.fa
} | gzip -n > Sbicolor_730_v5.0_PtMt.fa.gz
{
    gzip -dc Sbicolor_730_v5.1.gene_exons.gff3.gz |
        awk -F '\t' 'NR == FNR {drop[$1] = 1; next} /^#/ || !($1 in drop)' \
            "$here/mito-scaffolds.txt" -
    python3 -I "$here/organelle_gff.py" NC_008360.1.gff3 Mt | grep -v '^#'
    python3 -I "$here/organelle_gff.py" NC_008602.1.gff3 Pt | grep -v '^#'
} | gzip -n > Sbicolor_730_v5.1_PtMt.gene_exons.gff3.gz
gzip -dc Sbicolor_730_v5.1_PtMt.gene_exons.gff3.gz |
    awk -v OFS='\t' '$3 == "gene" {print $1, $4 - 1, $5, ".", 0, $7}' > Sbicolor_730_v5.1_PtMt.genes.bed
