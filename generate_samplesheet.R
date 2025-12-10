#!/usr/bin/env Rscript

library(optparse)
library(dplyr)

# Set up command line argument parsing
option_list <- list(
    make_option(c("-d", "--data_dir"), type = "character", default = "data",
                help = "Directory containing sample folders [default = %default]"),
    make_option(c("-o", "--output_file"), type = "character", default = "./data/samplesheet.csv",
                help = "Output samplesheet file name [default = %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Extract options
data_dir <- opt$data_dir
output_file <- opt$output_file

# Print the settings
cat("Data directory:", data_dir, "\n")
cat("Output file:", output_file, "\n")

# Get a list of all sample directories
samples <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)

# Initialize an empty data frame to store sample information
samplesheet <- data.frame(id = character(), fastq_dir = character(), stringsAsFactors = FALSE)

# Iterate over each sample directory
for (sample_dir in samples) {
    # Extract the sample ID from the directory name
    sample_id <- basename(sample_dir)

    # Extract the meta ID before the first "_S[ID]"
    meta_id_temp <- sub("_S\\d+.*$", "", basename(list.files(sample_dir)[1]))

    meta_id <- ifelse(meta_id_temp == '', basename(sample_dir), meta_id_temp) 
    
    # Check if the directory contains at least R1 and R2 files
    r1_files <- list.files(sample_dir, pattern = "_R1_.*\\.fastq\\.gz$", full.names = TRUE)
    r2_files <- list.files(sample_dir, pattern = "_R2_.*\\.fastq\\.gz$", full.names = TRUE)
    index_files <- list.files(sample_dir, pattern = "_I[12]*\\.fastq\\.gz$", full.names = TRUE)
    
    # Make sure all required files are present
    if (length(r1_files) == 0 || length(r2_files) == 0) {
        cat("Warning: Missing R1 or R2 files for sample", sample_id, "\n")
        next
    }

    # Add the directory to the samplesheet
    samplesheet <- samplesheet %>%
        add_row(id = meta_id, fastq_dir = sample_dir)
}

# Write the samplesheet to a CSV file
write.csv(samplesheet, output_file, row.names = FALSE, quote = FALSE)

cat("Samplesheet generated:", output_file, "\n")
