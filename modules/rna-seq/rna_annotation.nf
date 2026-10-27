include { dotenv } from 'plugin/nf-dotenv'

process RNA_ANNOTATION {
    cpus 1
    memory '4 GB'
    time '1h'
    conda "${moduleDir}/../../environments/GENESIS_RNA_COUNTS.yaml"
    container "${dotenv('GENESIS_RNA_COUNTS_IMAGE')}"
    tag "${ref.id}"
    input:
    tuple val(ref), path(annotation), path(chrom_sizes)
    output:
    tuple val(ref), path("${ref.id}.genes.gtf"), path("${ref.id}.rrna_genes.txt"),
        path("${ref.id}.max_intron.txt"), emit: annotation
    script:
    // Plant annotations come as GFF3 (Araport11, MaizeGDB, Phytozome) and sometimes as GTF, also
    // under a .gff3 name; gffread reads either and writes one GTF with gene_id on every exon.
    // The copy handling below only acts on GTF-style transcript_id attributes.
    """
    gzip -dcf '${annotation}' > annotation.txt
    # Lifted annotations (Liftoff, e.g. Col-CEN Araport11) reuse one transcript ID for copies at
    # other loci, such as rDNA and tRNA arrays, which gffread would join into one transcript with
    # huge introns. Copies are told apart by their transcript records' loci; a line outside all of
    # its transcript's records (e.g. a CDS-only copy) is dropped.
    awk -F '\\t' -v OFS='\\t' '
        function tid() {return match(\$9, /transcript_id "[^"]+"/) ? substr(\$9, RSTART + 15, RLENGTH - 16) : ""}
        NR == FNR {
            id = tid()
            if (id == "" || \$3 ~ /^(gene|exon|CDS|five_prime_UTR|three_prime_UTR|UTR|start_codon|stop_codon)\$/) next
            for (i = 1; i <= n[id]; i++) if (contig[id, i] == \$1 && \$4 <= last[id, i] && \$5 >= first[id, i]) {
                if (\$4 < first[id, i]) first[id, i] = \$4
                if (\$5 > last[id, i]) last[id, i] = \$5
                next
            }
            i = ++n[id]; contig[id, i] = \$1; first[id, i] = \$4; last[id, i] = \$5
            next
        }
        {
            id = tid()
            if (id == "" || !n[id] || \$3 == "gene") {print; next}
            for (i = 1; i <= n[id]; i++) if (contig[id, i] == \$1 && \$4 >= first[id, i] && \$5 <= last[id, i]) break
            if (i > n[id]) next
            if (i > 1) sub(/transcript_id "[^"]+"/, "transcript_id \\"" id "_copy" i "\\"", \$9)
            print
        }' annotation.txt annotation.txt > copies.txt
    gffread copies.txt -T -o genes.gtf
    # gffread keeps overlapping records that share an ID apart (e.g. a mature miRNA inside its
    # primary transcript); give each its own transcript ID too, or STAR would merge them.
    # Ensembl Plants GFF3 IDs carry gene: and transcript: prefixes; counts use the bare IDs.
    awk -F '\\t' -v OFS='\\t' '
        {sub(/gene_id "gene:/, "gene_id \\"", \$9); sub(/transcript_id "transcript:/, "transcript_id \\"", \$9)}
        \$3 == "transcript" {
            match(\$9, /transcript_id "[^"]+"/); id = substr(\$9, RSTART + 15, RLENGTH - 16)
            copy = ++seen[id]
        }
        copy > 1 {sub(/transcript_id "[^"]+"/, "transcript_id \\"" id "_dup" copy "\\"", \$9)}
        {print}' genes.gtf > '${ref.id}.genes.gtf'
    # Annotation and FASTA from different releases often name contigs differently (Chr1 vs 1).
    awk -F '\\t' 'NR == FNR {size[\$1] = 1; next} !(\$1 in size) {print \$1}' '${chrom_sizes}' \\
        '${ref.id}.genes.gtf' | sort -u > unknown.txt
    if [ -s unknown.txt ]; then
        printf 'Annotation contigs not in %s: %s\\n' '${ref.filename}' "\$(head -5 unknown.txt | paste -sd ' ')" >&2
        exit 1
    fi
    # rRNA genes: rRNA features (GFF3, Liftoff) or an rRNA biotype (Ensembl GTF); the gene is
    # gene_id in GTF-style attributes and Parent in GFF3 ones.
    awk -F '\\t' '\$3 == "rRNA" || \$9 ~ /(gene|transcript)_biotype "rRNA"/ {
        if (match(\$9, /gene_id "[^"]+"/)) id = substr(\$9, RSTART + 9, RLENGTH - 10)
        else if (match(\$9, /Parent=[^;]+/)) id = substr(\$9, RSTART + 7, RLENGTH - 7)
        else next
        sub(/^gene:/, "", id)
        print id
    }' annotation.txt | sort -u > '${ref.id}.rrna_genes.txt'
    # The longest annotated intron: STAR then aligns every annotated junction, and no novel one
    # longer than that, instead of its default of about 590 kb, which suits mammals.
    gffread '${ref.id}.genes.gtf' --table @exons 2> /dev/null | awk '{
        n = split(\$1, exons, ",")
        for (i = 2; i <= n; i++) {
            split(exons[i - 1], before, "-"); split(exons[i], after, "-")
            if (after[1] - before[2] - 1 > longest) longest = after[1] - before[2] - 1
        }
    } END {print longest + 0}' > '${ref.id}.max_intron.txt'
    """
    stub:
    """
    printf 'chr1\\tstub\\texon\\t1\\t100\\t.\\t+\\t.\\tgene_id "g1"; transcript_id "t1";\\n' > '${ref.id}.genes.gtf'
    touch '${ref.id}.rrna_genes.txt'
    printf '1000\\n' > '${ref.id}.max_intron.txt'
    """
}
