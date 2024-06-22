#!/bin/bash

# Create the directory structure
mkdir -p project/{data/{raw,processed,reference},results/{preprocessing,basecalling,demultiplexing,alignment,analysis/{species_identification,species_abundance,OTUs}},scripts/{analysis,utils}}

# Create necessary files
touch project/data/raw/sequencing_data.fastq
touch project/data/reference/reference_db.fasta

touch project/scripts/preprocess.nf
touch project/scripts/basecalling.sh
touch project/scripts/demultiplex.sh
touch project/scripts/align_reads.nf
touch project/scripts/analysis/species_identification.R
touch project/scripts/analysis/species_abundance.R
touch project/scripts/analysis/OTUs.R
touch project/scripts/utils/quality_control.R
touch project/scripts/utils/species_identification.R
touch project/scripts/utils/abundance_calculation.R
touch project/scripts/utils/plot_generation.R

touch project/nextflow.config
touch project/main.nf
touch project/requirements.txt

echo "Project structure created successfully."
