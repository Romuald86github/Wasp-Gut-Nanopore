library(dplyr)
library(tidyr)
library(readr)
library(Rsamtools)
library(Biostrings)

# Load aligned reads
aligned_bam <- "results/alignment/aligned_reads_sorted.bam"
reference_db <- "data/reference/reference_db.fasta"

# Function to calculate species abundance
calculate_abundance <- function(bam_file, reference_db) {
  # Load BAM file
  aln <- scanBam(bam_file)
  
  # Extract read sequences
  read_seqs <- sapply(aln[[1]]$seq, toString)
  
  # Load reference database
  ref_seqs <- readDNAStringSet(reference_db)
  
  # Perform alignment and calculate abundance
  results <- sapply(read_seqs, function(seq) {
    hits <- vmatchPattern(seq, ref_seqs)
    if (length(hits) > 0) {
      return(names(hits)[1])
    } else {
      return(NA)
    }
  })
  
  # Calculate abundance
  abundance <- table(results)
  
  # Create a data frame for results
  df <- as.data.frame(abundance)
  colnames(df) <- c("Species", "Abundance")
  return(df)
}

# Save results
calculate_abundance(aligned_bam, reference_db) %>%
  write_csv("results/analysis/species_abundance/species_abundance_results.csv")