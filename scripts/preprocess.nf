#!/usr/bin/env nextflow

params.raw_reads = 'data/raw/sequencing_data.fastq'

process preprocess {
    input:
    path raw_reads: params.raw_reads

    output:
    path 'results/preprocessing/filtered_reads.fastq.gz', emit: filtered_reads

    script:
    """
    # Create output directory if it doesn't exist
    mkdir -p results/preprocessing/

    # Quality filtering using NanoFilt
    gunzip -c ${raw_reads} | NanoFilt -q 7 -l 300 > results/preprocessing/filtered_reads.fastq

    # Adapter trimming using Porechop
    porechop -i results/preprocessing/filtered_reads.fastq -o results/preprocessing/trimmed_reads.fastq --threads 4

    # Further filtering to retain reads around 314bp
    seqtk seq -L 300 results/preprocessing/trimmed_reads.fastq | seqtk seq -M 330 > results/preprocessing/filtered_reads.fastq

    # Compress the final output
    gzip results/preprocessing/filtered_reads.fastq
    """
}