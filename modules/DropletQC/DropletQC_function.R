#!/usr/bin/env Rscript

library(DropletQC)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
thread <- as.integer(args[2])

# Check if the filtered matrix file exists
matrix_file <- file.path(input_file, "filtered_feature_bc_matrix.h5")
if (!file.exists(matrix_file)) {
    stop("Filtered matrix file not found: ", matrix_file)
}

message("STEP 1: Calculate nuclear fraction")
nf1 <- nuclear_fraction_tags(
    outs = input_file,
    tiles = 1,
    cores = thread,
    verbose = TRUE
)

nf1$cellnames <- rownames(nf1)
message(paste0("Calculated nuclear fraction for ", nrow(nf1), " cells."))


message("STEP 2: Identify empty droplets")

# Read sparse matrix directly from filtered_feature_bc_matrix.h5
h5 <- H5File$new(matrix_file, mode = "r")

# Load data
data <- h5[["matrix"]]
barcodes <- data[["barcodes"]][]
genes <- data[["features"]][["name"]][]
shape <- data[["shape"]][]

indptr <- data[["indptr"]][]
indices <- data[["indices"]][]
data_values <- data[["data"]][]

# Convert to sparse Matrix (dgCMatrix)
counts <- new("dgCMatrix",
              Dim = as.integer(shape),
              Dimnames = list(genes, barcodes),
              p = as.integer(indptr),
              i = as.integer(indices),
              x = as.numeric(data_values))

h5$close_all()

# Compute UMI counts per cell
umi_counts <- Matrix::colSums(counts)

nf_umi <- data.frame(
    cellnames = colnames(counts),
    nf = nf1$nuclear_fraction[match(colnames(counts), nf1$cellnames)],
    umi = umi_counts
)

nf_umi <- nf_umi[!is.na(nf_umi$nf), ]

ed <- identify_empty_drops(
    nf_umi = nf_umi[, c("nf", "umi")], 
    include_plot = TRUE, 
    pdf_png = "png", 
    plot_path ='.',
    plot_name = "empty_droplets_qc.png", 
    nf_rescue = 0.05, 
    umi_rescue = 1000
)
rownames(ed) <- nf_umi$cellnames

message("# # # Empty droplet output:")
print(head(ed))
print(table(ed$cell_status))
message(paste0("Identified empty cells for ", nrow(ed), " cells."))

write.csv(ed, "nf_ed_qc.csv", row.names = TRUE)
message("Output written to nf_ed_qc.csv")
