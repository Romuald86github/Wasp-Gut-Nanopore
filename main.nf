#!/usr/bin/env nextflow

params.raw_reads = 'data/raw/sequencing_data.fastq'
params.reference = 'data/reference/reference_db.fasta'

workflow {
    preprocess_reads(params.raw_reads)
    basecall_reads(preprocess_reads.out.filtered_reads)
    demultiplex_reads(basecall_reads.out.basecalled_reads)
    align_reads(demultiplex_reads.out.demultiplexed_reads, params.reference)
    analyze_results(align_reads.out.aligned_reads)
}

include { preprocess_reads } from './scripts/preprocess.nf'
include { basecall_reads } from './scripts/basecalling.sh'
include { demultiplex_reads } from './scripts/demultiplex.sh'
include { align_reads } from './scripts/align_reads.nf'

process analyze_results {
    input:
    path aligned_reads

    output:
    path 'results/analysis/'

    script:
    """
    # Create the output directory if it doesn't exist
    mkdir -p $params.analysis_out_dir

    # Run analysis scripts
    Rscript scripts/analysis/species_identification.R $aligned_reads $params.analysis_out_dir/species_identification/
    Rscript scripts/analysis/species_abundance.R $aligned_reads $params.analysis_out_dir/species_abundance/
    Rscript scripts/analysis/OTUs.R $aligned_reads $params.analysis_out_dir/OTUs/
    """
}