#!/bin/bash

filtered_reads=$1
output_dir=$2

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Basecalling using Guppy
guppy_basecaller \
    -i $filtered_reads \
    -s $output_dir \
    --config dna_r9.4.1_450bps_fast.cfg \
    --num_callers 4 \
    --gpu_runners_per_device 2 \
    --chunks_per_runner 512