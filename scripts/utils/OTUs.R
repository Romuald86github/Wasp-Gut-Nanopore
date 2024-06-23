library(dplyr)
library(tidyr)
library(readr)
library(Rsamtools)
library(Biostrings)
library(DECIPHER)

# Load aligned reads
aligned_bam <- as.character(snapshotParam("aligned_reads"))

# Function to identify OTUs
identify_otus <- function(bam_file) {
  # Load BAM file
  aln <- scanBam(bam_file)

  # Extract read sequences and sample IDs
  read_seqs <- sapply(aln[[1]]$seq, toString)
  sample_ids <- aln[[1]]$qname

  # Perform OTU clustering using DECIPHER
  otus <- decipher(read_seqs, method = "unoise", processors = 4, cutoff = 0.98)

  # Create a data frame for results
  df <- data.frame(SampleID = sample_ids, OTU = otus$name)
  return(df)
}

# Identify OTUs
otu_results <- identify_otus(aligned_bam)

# Save results
write_csv(otu_results, "$params.analysis_out_dir/OTUs/species_OTU_results.csv")