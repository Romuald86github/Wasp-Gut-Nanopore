process preprocess_reads {
    input:
    path raw_reads: params.raw_reads

    output:
    path 'results/preprocessing/filtered_reads.fastq.gz', emit: filtered_reads

    script:
    """
    # Create output directory if it doesn't exist
    mkdir -p $params.preprocess_out_dir

    # Quality filtering using NanoFilt
    # The -q 7 parameter sets the minimum quality score to 7
    # The -l 300 parameter sets the minimum read length to 300 bp
    gunzip -c ${raw_reads} | NanoFilt -q 7 -l 300 > $params.preprocess_out_dir/filtered_reads.fastq

    # Adapter trimming using Porechop
    # The --threads 4 parameter uses 4 CPU threads for the trimming
    porechop -i $params.preprocess_out_dir/filtered_reads.fastq -o $params.preprocess_out_dir/trimmed_reads.fastq --threads 4

    # Further filtering to retain reads around 314 bp
    # The seqtk seq -L 300 command filters for reads longer than 300 bp
    # The seqtk seq -M 330 command filters for reads shorter than 330 bp
    seqtk seq -L 300 $params.preprocess_out_dir/trimmed_reads.fastq | seqtk seq -M 330 > $params.preprocess_out_dir/filtered_reads.fastq

    # Compress the final output
    gzip $params.preprocess_out_dir/filtered_reads.fastq
    """
}