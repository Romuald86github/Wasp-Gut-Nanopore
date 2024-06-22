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

process preprocess_reads {
    input:
    path raw_reads

    output:
    path 'results/preprocessing/filtered_reads.fastq.gz', emit: filtered_reads

    script:
    """
    bash scripts/preprocess.sh $raw_reads results/preprocessing/filtered_reads.fastq.gz
    """
}

process basecall_reads {
    input:
    path filtered_reads

    output:
    path 'results/basecalling/', emit: basecalled_reads

    script:
    """
    bash scripts/basecalling.sh $filtered_reads results/basecalling/
    """
}

process demultiplex_reads {
    input:
    path basecalled_reads

    output:
    path 'results/demultiplexing/', emit: demultiplexed_reads

    script:
    """
    bash scripts/demultiplex.sh $basecalled_reads results/demultiplexing/
    """
}

process align_reads {
    input:
    path demultiplexed_reads
    path reference

    output:
    path 'results/alignment/', emit: aligned_reads

    script:
    """
    nextflow scripts/align_reads.nf --reads $demultiplexed_reads --reference $reference --outdir results/alignment/
    """
}

process analyze_results {
    input:
    path aligned_reads

    output:
    path 'results/analysis/'

    script:
    """
    # Run analysis scripts
    Rscript scripts/analysis/species_identification.R results/alignment/ results/analysis/species_identification/
    Rscript scripts/analysis/species_abundance.R results/alignment/ results/analysis/species_abundance/
    Rscript scripts/analysis/OTUs.R results/alignment/ results/analysis/OTUs/
    """
}