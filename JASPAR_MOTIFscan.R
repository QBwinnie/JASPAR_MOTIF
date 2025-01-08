# Install necessary packages if you haven't already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("JASPAR2020")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")

# Load the libraries
library(JASPAR2020)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)

# Define the path to your BED file (I had mine in CSV)
csv_file <- "~/your_input.csv"

# Read the CSV file
sequences_df <- read.csv(csv_file, stringsAsFactors = FALSE)

# Check the structure of the data frame
print(head(sequences_df))

# Create a DNAStringSet from the sequences
# Assuming the sequences are in the 'sequence' column
dna_seqs <- DNAStringSet(sequences_df$sequence)

# Optionally, set names if you have a sequence name column
if ("sequence_name" %in% colnames(sequences_df)) {
  names(dna_seqs) <- sequences_df$sequence_name
}

# Check the DNAStringSet
print(dna_seqs)

# Load JASPAR database and set options
# JASPAR requires a numeric taxonomy ID for species (9606 for Homo sapiens)
# collection="CORE" selects the core collection.
opts <- list()
opts[["species"]] <- 9606
opts[["collection"]] <- "CORE"

# pfm_list is a PFMatrixList
pfm_list <- TFBSTools::getMatrixSet(JASPAR2020, opts = opts)
# Convert each PFMatrix to a PWMatrix
pwm_list <- lapply(pfm_list, toPWM)  # Now each element is a PWMatrix

# Define a function to find motifs in a single DNA sequence
find_motifs <- function(dna_sequence, pwm_list, min_score="80%") {
  motif_hits <- list()
  
  # Loop over each PWMatrix and search in 'dna_sequence'
  for (motif_name in names(pwm_list)) {
    pwm <- pwm_list[[motif_name]]
    
    # searchSeq works with PWMatrix
    hits <- searchSeq(pwm, dna_sequence, min.score = min_score)
    
    if (length(hits) > 0) {
      motif_hits[[motif_name]] <- hits
    }
  }
  
  return(motif_hits)
}

# Find motifs in each DNA sequence
results <- lapply(dna_seqs, function(seq) {
  find_motifs(seq, pwm_list)
})

# Print or save the results
for (i in seq_along(results)) {
  cat(paste("Motifs for sequence", names(dna_seqs)[i], ":\n"))
  print(results[[i]])
  cat("\n")
}

# Optionally, save results as an RDS (more convenient for lists of lists)
saveRDS(results, file = "motif_results.rds")

#combine all motif scans of each seuqence in one csv file
combined_df <- data.frame()  # Create an empty data frame to store combined results
for (i in 1:40) {
  combined_df <- data.frame()  # Reset for each results[[i]]
  
  for (j in seq_along(results[[i]])) {
    site_set <- results[[i]][[j]]  # Extract each inner list
    
    if (class(site_set) == "SiteSet") {
      df <- as.data.frame(as.data.frame(site_set))
    } else {
      df <- as.data.frame(site_set)
    }
    
    combined_df <- rbind(combined_df, df)
  }
  
  # Create a filename for each combined CSV
  combined_filename <- paste0("combined_motif_results_", i, ".csv")
  
  # Write the combined data frame to a CSV file
  write.csv(combined_df, file = combined_filename, row.names = FALSE)
}
