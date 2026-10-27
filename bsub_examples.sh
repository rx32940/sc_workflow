#!/bin/bash
# Examples for running the workflow on the st.jude cluster under different situations.
# Copy the block you need into bsub.sh (or save it as its own file) and submit with:  bsub < bsub.sh
# Every example uses the same LSF header; edit the paths marked <...>.

#BSUB -P cellranger_workflow
#BSUB -R "rusage[mem=50GB]"
#BSUB -R "span[hosts=1]"
#BSUB -q standard
#BSUB -n 5
#BSUB -J cellranger_workflow
#BSUB -o cellranger_workflow.%J.out
#BSUB -e cellranger_workflow.%J.out

module load nextflow/23.10.0

WORKFLOW="/home/rxu28/software/sc_workflow"
OUTDIR=/path/to/output            # <...>
SAMPLESHEET=$OUTDIR/samplesheet.csv
mkdir -p $OUTDIR
cd $WORKFLOW

# Rscript for the samplesheet generators
RSCRIPT=/research/groups/northcgrp/home/common/Rachel/envs/R_nf_4.3.1/bin/Rscript


# =====================================================================================
# 1) FASTQs, RNA (cellranger count -> CellBender, Velocyto, DropletQC -> merge)
# =====================================================================================
# $RSCRIPT $WORKFLOW/generate_samplesheet.R -d <dir with one folder of FASTQs per sample> -o $SAMPLESHEET
# nextflow run main.nf --species human --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 2) FASTQs, multiome (cellranger-arc count -> downstream)
# =====================================================================================
# samplesheet columns: fastqs,sample,library_type   (see generate_samplesheet_multiome.R)
# nextflow run main.nf --modality arc --species mouse --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 3) FASTQs, ATAC only (cellranger-atac count, mapping only)
# =====================================================================================
# nextflow run main.nf --modality atac --species mouse --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 4) Cell Ranger already run (count) -> start at CellBender / Velocyto / DropletQC
# =====================================================================================
# every sample folder needs outs/ with raw + filtered h5, filtered_feature_bc_matrix/, BAM + .bai
# $RSCRIPT $WORKFLOW/generate_samplesheet_fromCR.R -d <dir containing Cell Ranger sample folders> -m count -o $SAMPLESHEET
# nextflow run main.nf --modality count --skip_cellranger --species mouse \
#   --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 5) Cell Ranger ARC already run (multiome)
# =====================================================================================
# BAM is gex_possorted_bam.bam
# $RSCRIPT $WORKFLOW/generate_samplesheet_fromCR.R -d <dir containing Cell Ranger ARC sample folders> -m arc -o $SAMPLESHEET
# nextflow run main.nf --modality arc --skip_cellranger --species mouse \
#   --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 6) BD Rhapsody output (default: filter the unfiltered matrix at --min_umi 500)
# =====================================================================================
# samplesheet:  id,rhapsody_dir   (rhapsody_dir holds *_Bioproduct.bam and *_RSEC_MolsPerCell_Unfiltered_MEX.zip)
# --gtf should match your Rhapsody reference; without it the Cell Ranger GENCODE GTF for --species is used
# nextflow run main.nf --modality rhapsody --species human \
#   --gtf <genes.gtf matching the Rhapsody reference> --min_umi 500 \
#   --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 7) BD Rhapsody with a different UMI floor for the "filtered" barcode set
# =====================================================================================
# lower = more barcodes for DropletQC/Velocyto (slower); higher = closer to the called cells
# nextflow run main.nf --modality rhapsody --species human --min_umi 200 \
#   --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 8) BD Rhapsody with your own h5 files
# =====================================================================================
# samplesheet:  id,rhapsody_dir,raw_h5,filtered_h5   (either h5 column may be left empty)
# barcodes in the h5 files must be the Rhapsody cell labels found in the BAM CB tag
# nextflow run main.nf --modality rhapsody --species human \
#   --samplesheet $SAMPLESHEET --outdir $OUTDIR -resume


# =====================================================================================
# 9) Merge only a subset of samples in the final step
# =====================================================================================
# append to any command above:
#   --selected_samples sampleA,sampleB
#   --selected_batches_file /path/to/batches.txt     # one id per line


# =====================================================================================
# 10) Resume after a failure / change a setting and re-run
# =====================================================================================
# add -resume (already in the examples) and re-submit with the SAME -w work dir; to use a separate work dir:
#   -w /path/to/work
