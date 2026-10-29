include { SRA_PREFETCH } from './sra_prefetch'
include { SRA_FASTQ } from './sra_fastq'

// One task per run: download it, then split out its two genomic mates and its barcode read.
workflow SCATAC_SRA {
    take:
    runs  // tuple(meta, run accession)

    main:
    SRA_PREFETCH(runs)
    SRA_FASTQ(SRA_PREFETCH.out.sra)

    emit:
    reads = SRA_FASTQ.out.reads  // tuple(meta, run, R1, R2, barcode FASTQs)
}
