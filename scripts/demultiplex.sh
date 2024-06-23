#!/bin/bash

basecalled_reads=$1
output_dir=$2

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Demultiplexing using qcat
# The --threads 4 parameter uses 4 CPU threads for the demultiplexing
qcat -i $basecalled_reads -o $output_dir --threads 4