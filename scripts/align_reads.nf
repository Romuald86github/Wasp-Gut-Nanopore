#!/usr/bin/env nextflow

params.reads = 'results/demultiplexing/'
params.reference = 'data/reference/reference_db.fasta'
params.outdir = 'results/alignment/'

process align_reads {
    input:
    path reads from params.reads
    path reference from params.reference

    output:
    path params.outdir

    script:
    """
    # Create the output directory if it doesn't exist
    mkdir -p $params.outdir

    # Align the reads to the reference database using minimap2
    minimap2 -ax map-ont $reference $reads > $params.outdir/aligned_reads.sam

    # Convert the SAM file to BAM format
    samtools view -b -o $params.outdir/aligned_reads.bam $params.outdir/aligned_reads.sam

    # Sort the BAM file
    samtools sort -o $params.outdir/aligned_reads_sorted.bam $params.outdir/aligned_reads.bam

    # Index the sorted BAM file
    samtools index $params.outdir/aligned_reads_sorted.bam
    """
}