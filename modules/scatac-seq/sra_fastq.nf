include { dotenv } from 'plugin/nf-dotenv'

process SRA_FASTQ {
    // Each conversion needs temporary disk of about the FASTQs' size (about 100 GB for 116 million
    // 10x spots), so only a few run at once.
    maxForks 2
    cpus 4
    memory '4 GB'
    time { 4.h * task.attempt }
    // retry with more memory and time when killed (OOM, walltime, or preemption); any other
    // failure is a run without the expected reads, which converting again cannot fix
    errorStrategy { task.exitStatus in 130..145 ? 'retry' : 'terminate' }
    maxRetries 2
    conda "${moduleDir}/../../environments/GENESIS_SRA.yaml"
    container "${dotenv('GENESIS_SRA_IMAGE')}"
    tag "${meta.id}:${run}"
    input:
    tuple val(meta), val(run), path(sra)
    output:
    tuple val(meta), val(run), path("${run}.R1.fastq.gz"), path("${run}.R2.fastq.gz"),
        path("${run}.barcode.fastq.gz"), emit: reads
    script:
    // Reads are told apart by length, as submitters order them differently (e.g. 8/50/16/49,
    // 50/16/50, 51/51/16): the barcode read is meta.barcode_length long, the genomic mates are the
    // two longest, and an 8-base sample index read is dropped.
    """
    # sra-tools keeps its settings outside the home folder.
    export NCBI_SETTINGS="\$PWD/ncbi-settings.mkfg"
    # fasterq-dump reports how many spots it read and reads it wrote.
    fasterq-dump --split-files --include-technical --threads ${task.cpus} --temp tmp \\
        --outdir split '${sra}' 2>&1 | tee fasterq-dump.log >&2
    rm -r tmp
    for file in split/*.fastq; do
        printf '%s\\t%s\\n' "\$file" "\$(head -n 2 "\$file" | tail -n 1 | tr -d '\\n' | wc -c)"
    done | sort -t \$'\\t' -k2,2n > reads.tsv
    cat reads.tsv >&2
    barcode=\$(awk -F '\\t' '\$2 == ${meta.barcode_length} {print \$1}' reads.tsv)
    if [ "\$(printf '%s' "\$barcode" | grep -c .)" != 1 ]; then
        printf '%s has no single %s-base barcode read; is it a 10x scATAC run?\\n' '${run}' \\
            '${meta.barcode_length}' >&2
        exit 1
    fi
    mates=\$(awk -F '\\t' '\$2 >= 30' reads.tsv | tail -n 2 | cut -f1)
    if [ "\$(printf '%s\\n' "\$mates" | grep -c .)" != 2 ]; then
        printf '%s has no two genomic reads of 30 bases or more\\n' '${run}' >&2
        exit 1
    fi
    # Every file must hold every spot: reads written = spots x files.
    spots=\$(awk -F ':' '/spots read/ {gsub(/[^0-9]/, "", \$2); print \$2}' fasterq-dump.log)
    written=\$(awk -F ':' '/reads written/ {gsub(/[^0-9]/, "", \$2); print \$2}' fasterq-dump.log)
    if [ -z "\$spots" ] || [ "\$written" != "\$(( spots * \$(wc -l < reads.tsv) ))" ]; then
        printf '%s: %s reads written for %s spots in %s files\\n' '${run}' "\$written" "\$spots" \\
            "\$(wc -l < reads.tsv)" >&2
        exit 1
    fi
    first=\$(printf '%s\\n' "\$mates" | head -n 1)
    second=\$(printf '%s\\n' "\$mates" | tail -n 1)
    # The lower-numbered genomic read is R1.
    if [ "\$first" \\> "\$second" ]; then tmp=\$first; first=\$second; second=\$tmp; fi
    pigz -p ${task.cpus} -c "\$first" > '${run}.R1.fastq.gz'
    pigz -p ${task.cpus} -c "\$second" > '${run}.R2.fastq.gz'
    pigz -p ${task.cpus} -c "\$barcode" > '${run}.barcode.fastq.gz'
    rm -r split
    """
    stub:
    """
    printf '@r1\\nACGT\\n+\\nIIII\\n' | gzip > '${run}.R1.fastq.gz'
    printf '@r1\\nACGT\\n+\\nIIII\\n' | gzip > '${run}.R2.fastq.gz'
    printf '@r1\\nACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIII\\n' | gzip > '${run}.barcode.fastq.gz'
    """
}
