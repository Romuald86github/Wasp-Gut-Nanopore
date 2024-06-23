#!/bin/bash

filtered_reads=$1
output_dir=$2

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Basecalling using Guppy
# The --cpu_threads_per_caller 4 parameter uses 4 CPU threads per basecaller
# The --chunks_per_runner 512 parameter sets the number of chunks per runner to 512
guppy_basecaller \
    -i $filtered_reads \
    -s $output_dir \
    --config dna_r9.4.1_450bps_fast.cfg \
    --cpu_threads_per_caller 4 \
    --chunks_per_runner 512