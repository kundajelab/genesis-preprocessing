include { dotenv } from 'plugin/nf-dotenv'

process SCATAC_FRAGMENTS {
    cpus 2
    memory { 2.GB * task.attempt }
    time { 1.h * (1 + fragments.size() / 1e9) * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption)
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/DAP_SEQ_ALIGNMENT.yaml"
    container "${dotenv('DAP_SEQ_ALIGNMENT_IMAGE')}"
    tag "${meta.id}"
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}", mode: 'copy',
        pattern: '*.fragments.tsv.gz*'
    publishDir "${params.run_dir}/output/${meta.species}/${meta.id}/qc", mode: 'copy',
        pattern: '{*.cells.tsv,*_mqc.tsv}'
    input:
    tuple val(meta), path(fragments), path(cells), path(chrom_sizes)
    output:
    tuple val(meta), path("${meta.id}.fragments.tsv.gz"), path("${meta.id}.fragments.tsv.gz.tbi"),
        emit: fragments
    tuple val(meta), path('groups/*.fragments.bed'), emit: groups
    tuple val(meta), path("${meta.id}.cells.tsv"), emit: cells
    tuple val(meta), path("${meta.id}.scatac_library_mqc.tsv"), emit: summary
    script:
    def organelles = meta.organelles.join(',')
    // Fragments are kept exactly as the producer wrote them (10x: BED0 half-open, Tn5 offsets
    // already applied); only rows for listed cells on nuclear contigs continue.
    """
    # A misspelled organelle contig would silently count organelle fragments as nuclear.
    for contig in \$(printf '%s' '${organelles}' | tr ',' ' '); do
        awk -v contig="\$contig" '\$1 == contig {found = 1} END {exit !found}' '${chrom_sizes}' || {
            printf 'Organelle contig %s is not in %s\\n' "\$contig" '${meta.reference_fasta}' >&2
            exit 1
        }
    done
    # The cell list: a barcode<TAB>group header, then one row per cell. Groups name output files.
    awk -F '\\t' '
        NR == 1 {if (\$0 != "barcode\\tgroup") {print "cells: expected header barcode<TAB>group" > "/dev/stderr"; exit 1}; print; next}
        \$2 !~ /^[A-Za-z0-9_][A-Za-z0-9_.-]*\$/ {
            print "cells line " NR ": groups use only letters, digits, _ . -" > "/dev/stderr"; exit 1
        }
        \$1 in first_line {
            print "cells line " NR ": barcode " \$1 " is listed again (first on line " first_line[\$1] \\
                "); each cell must appear once, in one group" > "/dev/stderr"; exit 1
        }
        {first_line[\$1] = NR; print}
    ' '${cells}' > cells.tsv
    # Fragments from another assembly stop here, then keep listed cells. Organelle fragments and
    # unlisted barcodes are counted, not kept.
    bgzip --decompress --stdout --threads ${task.cpus} '${fragments}' | awk -F '\\t' -v OFS='\\t' \\
        -v organelles='${organelles}' -v sizes='${chrom_sizes}' -v cells=cells.tsv '
        BEGIN {
            while ((getline line < sizes) > 0) {split(line, f, "\\t"); length_of[f[1]] = f[2]}
            n = split(organelles, names, ","); for (i = 1; i <= n; i++) organelle[names[i]] = 1
            getline line < cells
            while ((getline line < cells) > 0) {split(line, f, "\\t"); group[f[1]] = f[2]}
        }
        /^#/ {next}
        {
            if (!(\$1 in length_of) || \$3 + 0 > length_of[\$1] + 0) {
                print "fragments line " NR ": not a fragment on the reference: " \$0 > "/dev/stderr"
                exit 1
            }
            rows++
            if (!(\$4 in group)) {unlisted++; next}
            if (\$1 in organelle) {organelle_rows[\$4]++; organelle_total++; next}
            print
        }
        END {
            printf "%d\\t%d\\t%d\\n", rows, unlisted, organelle_total > "counts.tsv"
            for (barcode in organelle_rows) print barcode, organelle_rows[barcode] > "organelle.tsv"
        }' | LC_ALL=C sort -S 1G -T . -t \$'\\t' -k1,1 -k2,2n -k3,3n -k4,4 > kept.tsv
    touch organelle.tsv
    bgzip --threads ${task.cpus} -c kept.tsv > '${meta.id}.fragments.tsv.gz'
    tabix -p bed '${meta.id}.fragments.tsv.gz'
    # Per cell: nuclear fragments, their read-pair support, and organelle fragments.
    awk -F '\\t' -v OFS='\\t' -v organelle=organelle.tsv -v cells=cells.tsv '
        BEGIN {
            while ((getline line < organelle) > 0) {split(line, f, "\\t"); org[f[1]] = f[2]}
        }
        {fragments[\$4]++; support[\$4] += \$5}
        END {
            print "barcode", "group", "fragments", "read_pairs", "organelle_fragments", "organelle_fraction"
            getline line < cells
            while ((getline line < cells) > 0) {
                split(line, f, "\\t"); b = f[1]; n = fragments[b] + 0; o = org[b] + 0
                print b, f[2], n, support[b] + 0, o, (n + o ? o / (n + o) : "NA")
            }
        }' kept.tsv > '${meta.id}.cells.tsv'
    # One pseudobulk per group: its cells' fragments, still sorted.
    mkdir groups
    awk -F '\\t' -v OFS='\\t' -v cells=cells.tsv -v id='${meta.id}' '
        BEGIN {getline line < cells; while ((getline line < cells) > 0) {split(line, f, "\\t"); group[f[1]] = f[2]}}
        {print \$1, \$2, \$3, \$4, \$5 > ("groups/" id "." group[\$4] ".fragments.bed")}' kept.tsv
    # No listed barcode matching any fragment usually means a barcode suffix mismatch, e.g. -1.
    test -s kept.tsv || {
        printf 'No fragments for any listed barcode; check that the barcodes match the fragments\\n' >&2
        exit 1
    }
    {
        printf '# id: scatac_library\\n# section_name: scATAC libraries\\n# plot_type: table\\n'
        printf 'sample\\tsource_fragments\\tunlisted_barcode_fragments\\torganelle_fragments\\tcells\\tcell_fragments\\n'
        printf '%s\\t%s\\t%s\\n' '${meta.id}' "\$(cut -f1-3 counts.tsv)" \\
            "\$(awk 'NR > 1 {cells++; n += \$3} END {print cells "\\t" n + 0}' '${meta.id}.cells.tsv')"
    } > '${meta.id}.scatac_library_mqc.tsv'
    """
    stub:
    """
    mkdir groups
    printf 'chr1\\t100\\t300\\tAAA-1\\t1\\n' > 'groups/${meta.id}.all.fragments.bed'
    bgzip -c 'groups/${meta.id}.all.fragments.bed' > '${meta.id}.fragments.tsv.gz'
    touch '${meta.id}.fragments.tsv.gz.tbi'
    printf 'barcode\\tgroup\\n' > '${meta.id}.cells.tsv'
    printf 'sample\\n${meta.id}\\n' > '${meta.id}.scatac_library_mqc.tsv'
    """
}
