process align_reads {
    input:
    path demultiplexed_reads
    path reference

    output:
    path 'results/alignment/', emit: aligned_reads

    script:
    """
    # Create the output directory if it doesn't exist
    mkdir -p $params.alignment_out_dir

    # Align the reads to the reference database using minimap2
    # The -ax map-ont parameter specifies the alignment mode for Oxford Nanopore data
    minimap2 -ax map-ont $reference $demultiplexed_reads > $params.alignment_out_dir/aligned_reads.sam

    # Convert the SAM file to BAM format
    samtools view -b -o $params.alignment_out_dir/aligned_reads.bam $params.alignment_out_dir/aligned_reads.sam

    # Sort the BAM file
    samtools sort -o $params.alignment_out_dir/aligned_reads_sorted.bam $params.alignment_out_dir/aligned_reads.bam

    # Index the sorted BAM file
    samtools index $params.alignment_out_dir/aligned_reads_sorted.bam
    """
}