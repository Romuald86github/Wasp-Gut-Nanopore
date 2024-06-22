library(dplyr)
library(tidyr)
library(readr)
library(Rsamtools)
library(Biostrings)

# Load aligned reads
aligned_bam <- "results/alignment/aligned_reads_sorted.bam"

# Function to identify OTUs
identify_otus <- function(bam_file) {
  # Load BAM file
  aln <- scanBam(bam_file)
  
  # Extract read sequences
  read_seqs <- sapply(aln[[1]]$seq, toString)
  
  # Perform clustering to identify OTUs (simplified example)
  otus <- unique(read_seqs)
  
  # Create a data frame for results
  df <- data.frame(SampleID = seq_along(otus), OTU = otus)
  return(df)
}

# Save results
identify_otus(aligned_bam) %>%
  write_csv("results/analysis/OTUs/species_OTU_results.csv")