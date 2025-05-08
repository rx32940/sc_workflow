library(tidyverse)
library(data.table)

search_files <- function(indir, samples, pattern) {
  indir <- indir[[1]]
  patterns <- sprintf(pattern,samples)
  allfiles <- list.files(indir, recursive=TRUE, full.names=TRUE)
  matches <- lapply(patterns, function(p) allfiles[grepl(p,allfiles)])
  return(matches)
}

metadata <- fread("/home/vmartin/projects/snmcseq/multiome/northcott_cerebellumAtlas_metadata.csv")
allsamples <- paste(c(metadata$sample,str_replace(metadata$sample,"_","-")), collapse = "|")
datadir <- "/research/groups/northcgrp/home/common/HumanCBMB_singlecell/Cerebellum/scATAC/Raw/Nuclei/Droplet/Raw"
allfiles <- list.files(datadir, pattern = "\\.gz$", recursive=TRUE, full.names=TRUE)

df <- data.frame(fastq_file = allfiles) %>%
  mutate(
    sample =  sapply(strsplit(basename(fastq_file),"_S"), `[`, 1),
    type = case_when(
      grepl("_I1_", fastq_file) ~ "fastqi",  
      grepl("_R1_", fastq_file) ~ "fastq1", 
      grepl("_R2_", fastq_file) ~ "fastq2", 
      grepl("_R3_", fastq_file) ~ "fastq3"
    )
  ) %>% pivot_wider(names_from = type, values_from = fastq_file,values_fn = list) %>%
  unnest(c(fastqi,fastq1,fastq2,fastq3)) %>%
  dplyr::select(sample,everything())

fwrite(df,"/home/vmartin/projects/pipeline/cellranger/samplesheet_atac.csv")


