library(dplyr)
library(tidyr)
library(readr)
library(Rsamtools)
library(Biostrings)
library(DECIPHER)

# Load aligned reads
aligned_bam <- as.character(snapshotParam("aligned_reads"))
reference_db <- "data/reference/reference_db.fasta"

# Function to run species identification using DECIPHER
species_identification <- function(bam_file, reference_db) {
  # Load BAM file
  aln <- scanBam(bam_file)

  # Extract read sequences and sample IDs
  read_seqs <- sapply(aln[[1]]$seq, toString)
  sample_ids <- aln[[1]]$qname

  # Load reference database
  ref_seqs <- readDNAStringSet(reference_db)

  # Perform species identification using DECIPHER
  results <- decipher(read_seqs, ref_seqs, method = "blastn", processors = 4)

  # Create a data frame for results
  df <- data.frame(SampleID = sample_ids, Species = results$name)
  return(df)
}

# Run species identification
species_identification_results <- species_identification(aligned_bam, reference_db)

# Save results
write_csv(species_identification_results, "$params.analysis_out_dir/species_identification/species_identification_results.csv")