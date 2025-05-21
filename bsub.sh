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
WORK="/home/rxu28/software/sc_workflow"

# ######### 1) step 1: generate samplesheet #########
# source ~/miniconda3/bin/activate /research/groups/northcgrp/home/common/Rachel/envs/R_nf_4.3.1

# DATADIR="/research/groups/northcgrp/projects/northcgrp_hartwell/common/illumina"
# PROJECT="northcgrp_850454_10XscRNAseq-1"

# Rscript $WORKFLOW/generate_samplesheet.R -d $DATADIR/$PROJECT -o $WORK/data/samplesheet.csv

# source  ~/miniconda3/bin/deactivate

# ######### 1) step 2: run workflow #########
cd $WORKFLOW

OUTDIR=/home/rxu28/software/sc_workflow/Processed

mkdir -p $OUTDIR
nextflow run main.nf --species human --samplesheet $WORK/data/samplesheet.csv --outdir $OUTDIR -resume


