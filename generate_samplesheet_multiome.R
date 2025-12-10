#!/usr/bin/env Rscript

library(optparse)
library(dplyr)

option_list <- list(
    make_option(c("-r", "--rna_dir"), type = "character", 
                help = "Directory containing RNA sample folders"),
    make_option(c("-a", "--atac_dir"), type = "character", 
                help = "Directory containing ATAC sample folders"),
    make_option(c("-o", "--output_file"), type = "character", default = "./libraries.csv",
                help = "Output libraries file name [default = %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$rna_dir) || is.null(opt$atac_dir)) {
    stop("Both --rna_dir and --atac_dir must be specified")
}

cat("RNA directory:", opt$rna_dir, "\n")
cat("ATAC directory:", opt$atac_dir, "\n")
cat("Output file:", opt$output_file, "\n")

# Get sample directories
rna_samples <- list.dirs(opt$rna_dir, full.names = TRUE, recursive = FALSE)
atac_samples <- list.dirs(opt$atac_dir, full.names = TRUE, recursive = FALSE)

# Initialize libraries dataframe
libraries <- data.frame(fastqs = character(), sample = character(), library_type = character(), stringsAsFactors = FALSE)

# Process RNA samples
for (sample_dir in rna_samples) {
    sample_id <- basename(sample_dir)
    # Extract meta ID (remove _merged suffix if present)
    meta_id <- sub("_merged$", "", sample_id)
    
    libraries <- libraries %>%
        add_row(fastqs = sample_dir, sample = meta_id, library_type = "Gene Expression")
}

# Process ATAC samples
for (sample_dir in atac_samples) {
    sample_id <- basename(sample_dir)
    meta_id <- sub("_merged$", "", sample_id)
    
    libraries <- libraries %>%
        add_row(fastqs = sample_dir, sample = meta_id, library_type = "Chromatin Accessibility")
}

# Write libraries.csv
write.csv(libraries, opt$output_file, row.names = FALSE, quote = FALSE)

cat("Libraries file generated:", opt$output_file, "\n")
cat("Total libraries:", nrow(libraries), "\n")