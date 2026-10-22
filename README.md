# 10X Single-Cell RNA Preprocessing Pipeline

- **current implementation is only for Northcott lab to use on st.jude HPC cluster.**

___

This Nextflow pipeline automates pre-processing of 10X Genomics single-cell RNA-seq by running the following tools:

* `Cell Ranger`: read alignment and quantification
* `CellBender`: ambient RNA removal
* `Velocyto`: spliced/unspliced quantification for RNA velocity
* `DropletQC`: nuclear fraction and empty droplet detection

**use `downstream/downstream_filtering.ipynb` to further filter the data.**


## Pipeline Overview

```
samplesheet.csv → CellRanger → CellBender, Velocyto, DropletQC → Merge -> .h5ad
```

- will also take multiome and ATAC** file

  _** mapping only_

## Input Sample Sheet Format

`samplesheet.csv` (comma-separated, header required):


*IMPORTANT:

- data dir name cannot be same as id
- fastq files have to follow cellranger's default naming protocol: 
`<id>_S<ANY_NUM>_L<lane_number>_R<read_number>_001.fastq.gz`

---

```csv
id,fastq_dir
3300000_dummy_scRNA,./data/3300000
3300001_dummy_scRNA,./data/3300001
```

* for multiome
```csv
fastqs,sample,library_type
./data/3300000,3300000_dummy_scRNA,Gene Expression
./data/3300000,3300000_dummy_scRNA,Chromatin Accessibility
./data/3300001,3300001_dummy_scRNA,Gene Expression
./data/3300001,3300001_dummy_scRNA,Chromatin Accessibility
```
## Running the Pipeline

```bash
nextflow run main.nf \
  --samplesheet ./data/samplesheet.csv \
  --species human \
  --outdir /path/to/output \
  -resume
```

* for multiome
```bash
nextflow run main.nf --modality arc \
--species mouse \
--samplesheet ./data/samplesheet.csv \
--outdir /path/to/output
```

* for ATAC
```bash
nextflow run main.nf --modality atac \
--species mouse \
--samplesheet ./data/samplesheet.csv \
--outdir /path/to/output
```

* for existing Cell Ranger output (count or arc; skips mapping)
```csv
id,cellranger_dir
sampleA,/path/to/cellranger/sampleA
```
`cellranger_dir` is the sample folder containing `outs/` (or the `outs/` folder itself); it needs
`raw_feature_bc_matrix.h5`, `filtered_feature_bc_matrix.h5` and `filtered_feature_bc_matrix/`, and the BAM
with its `.bai` (`possorted_genome_bam.bam` for count, `gex_possorted_bam.bam` for arc). Files are symlinked
into the work directory, so the original run is never modified.
```bash
nextflow run main.nf --modality count --skip_cellranger \
--samplesheet ./data/samplesheet.csv --outdir /path/to/output   # use --modality arc for multiome
```

to submit a job to the st.jude cluster:

- edit `bsub.sh` file in this repo

```bash
bsub < bsub.sh
```



### Optional Parameters

| Parameter            | Description                                                  | Default                         |
| -------------------- | ------------------------------------------------------------ | ------------------------------- |
| `--genome_ref`       | Custom Cell Ranger reference path                            | Auto-selected based on species  |
| `--gtf`              | GTF file (for Velocyto); auto-detected if not provided       | `${genome_ref}/genes/genes.gtf` |
| `--species`          | `human` or `mouse`; used to auto-select genome reference     | `human`                       |
| `--outdir`           | Output directory root                                        | `./results`                     |
| `--selected_samples` | Optional comma-separated sample IDs to merge in final step   | all                             |

---

## Output Structure

```
outdir/
├── cellranger/
│   └── <sample_id>/
├── cellbender/
│   └── <sample_id>_cb/
├── Velocyto/
│   └── <sample_id>.loom
├── DropletQC/
│   └── <sample_id>/nf_ed_qc.csv
└── writes/
    └── merged_cts.h5ad
```

---

# Downstream filtering

- Customized filtering threshold is recommanded based on metrics produced by **cellbender** and **DropletQC**
- An example filtering workflow is available in [downstream/downstream_filtering.ipynb](https://github.com/rx32940/sc_workflow/blob/main/downstream/downstream_filtering.ipynb)

---

## Resuming

Use `-resume` to skip already completed steps:

```
nextflow run main.nf -resume
```

