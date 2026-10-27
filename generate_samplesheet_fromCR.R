#!/usr/bin/env Rscript

library(optparse)

# Build a samplesheet (id,cellranger_dir) from existing Cell Ranger count/arc output,
# for use with: nextflow run main.nf --skip_cellranger --modality <count|arc> ...
option_list <- list(
    make_option(c("-d", "--cr_dir"), type = "character",
                help = "Directory containing Cell Ranger sample folders (each with an outs/ subfolder)"),
    make_option(c("-m", "--modality"), type = "character", default = "count",
                help = "Cell Ranger output type: count or arc [default = %default]"),
    make_option(c("-o", "--output_file"), type = "character", default = "./data/samplesheet.csv",
                help = "Output samplesheet file name [default = %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$cr_dir)) {
    stop("--cr_dir must be specified")
}
if (!opt$modality %in% c("count", "arc")) {
    stop("--modality must be 'count' or 'arc'")
}

cr_dir <- opt$cr_dir
output_file <- opt$output_file
bam_name <- ifelse(opt$modality == "arc", "gex_possorted_bam.bam", "possorted_genome_bam.bam")

cat("Cell Ranger directory:", cr_dir, "\n")
cat("Modality:", opt$modality, "\n")
cat("Output file:", output_file, "\n")

# Files the workflow needs inside outs/
required <- c("raw_feature_bc_matrix.h5", "filtered_feature_bc_matrix.h5",
              "filtered_feature_bc_matrix", bam_name, paste0(bam_name, ".bai"))

samples <- list.dirs(cr_dir, full.names = TRUE, recursive = FALSE)

samplesheet <- data.frame(id = character(), cellranger_dir = character(), stringsAsFactors = FALSE)

for (sample_dir in samples) {
    sample_id <- basename(sample_dir)
    outs_dir <- file.path(sample_dir, "outs")

    if (!dir.exists(outs_dir)) {
        cat("Warning: no outs/ folder for sample", sample_id, "- skipping\n")
        next
    }

    missing <- required[!file.exists(file.path(outs_dir, required))]
    if (length(missing) > 0) {
        cat("Warning: sample", sample_id, "is missing:", paste(missing, collapse = ", "), "- skipping\n")
        next
    }

    samplesheet[nrow(samplesheet) + 1, ] <- c(sample_id, normalizePath(sample_dir))
}

if (nrow(samplesheet) == 0) {
    stop("No complete Cell Ranger ", opt$modality, " samples found in ", cr_dir)
}

write.csv(samplesheet, output_file, row.names = FALSE, quote = FALSE)

cat("Samplesheet generated:", output_file, "with", nrow(samplesheet), "samples\n")
