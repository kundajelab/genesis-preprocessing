"""Turn a RefSeq organelle GFF3 into gene -> transcript -> exon/CDS models on a renamed contig.

RefSeq gives protein-coding organelle genes as gene -> CDS without a transcript, and a
trans-spliced gene (e.g. nad1, rps12) as several `gene` rows sharing one ID, one per
piece (`part=N`), on either strand. Each piece becomes its own transcript (`<id>_partN`)
under one gene, so reads still count to the gene but no annotated intron joins pieces
that are far apart or on opposite strands. Pseudogenes and non-gene features are dropped.

    python3 -I organelle_gff.py REFSEQ.gff3 NEW_SEQID > out.gff3
"""

import sys
from collections import defaultdict
from urllib.parse import quote, unquote


def attributes(text):
    return dict(item.split("=", 1) for item in text.strip().split(";") if item)


def render(values):
    return ";".join(f"{k}={quote(unquote(v), safe=' ,:/()')}" for k, v in values.items())


def main(path, seqid):
    genes = defaultdict(list)
    children = defaultdict(list)
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        row = line.rstrip("\n").split("\t")
        row[0] = seqid
        row[3], row[4] = int(row[3]), int(row[4])
        row[8] = attributes(row[8])
        if row[2] == "gene":
            genes[row[8]["ID"]].append(row)
        elif row[2] in ("tRNA", "rRNA", "exon", "CDS") and "Parent" in row[8]:
            children[row[8]["Parent"]].append(row)
    out = [f"##gff-version 3\n"]
    for gene_id, parts in genes.items():
        parts.sort(key=lambda r: int(r[8].get("part", 1)))
        rnas = [r for r in children[gene_id] if r[2] in ("tRNA", "rRNA")]
        if rnas:
            models = [(r[2], r[8]["ID"], children[r[8]["ID"]]) for r in rnas]
        else:
            models = [("mRNA", "rna-" + gene_id.removeprefix("gene-"), children[gene_id])]
        keep = {k: v for k, v in parts[0][8].items() if k != "part"}
        for part_number, part in enumerate(parts, 1):
            gene = [*part[:8], {**keep, "ID": gene_id}]
            if len(parts) > 1:
                gene[8]["part"] = str(part_number)
            out.append(gene)
            for kind, rna_id, rows in models:
                inside = [
                    r for r in rows
                    if r[6] == part[6] and part[3] <= r[3] and r[4] <= part[4]
                ]
                if not inside:
                    continue
                # Organelle genes have no UTRs: a piece given only as CDS gets its exons from it.
                if not any(r[2] == "exon" for r in inside):
                    inside += [[*r[:2], "exon", *r[3:7], ".", {}] for r in inside if r[2] == "CDS"]
                tid = rna_id if len(parts) == 1 else f"{rna_id}_part{part_number}"
                start = min(r[3] for r in inside)
                end = max(r[4] for r in inside)
                note = {k: keep[k] for k in ("gene", "gene_biotype") if k in keep}
                out.append([seqid, "RefSeq", kind, start, end, ".", part[6], ".",
                            {"ID": tid, "Parent": gene_id, **note}])
                for kind_rows in ("exon", "CDS"):
                    rows_of = sorted((r for r in inside if r[2] == kind_rows), key=lambda r: r[3])
                    for n, r in enumerate(rows_of, 1):
                        extra = {k: r[8][k] for k in ("product", "protein_id") if k in r[8]}
                        out.append([*r[:8], {"ID": f"{tid}-{kind_rows.lower()}{n}", "Parent": tid, **extra}])
    for row in out:
        if isinstance(row, str):
            sys.stdout.write(row)
        else:
            print("\t".join(map(str, row[:8])), render(row[8]), sep="\t")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
