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

The workflow can be started at three points (details in the sections below):

| Starting point | Samplesheet | Flags |
| --- | --- | --- |
| FASTQs (Cell Ranger runs first) | `id,fastq_dir` (multiome: `fastqs,sample,library_type`) | `--modality rna\|arc\|atac` |
| Existing Cell Ranger count / arc output | `id,cellranger_dir` (make it with `generate_samplesheet_fromCR.R`) | `--skip_cellranger --modality count\|arc` |
| BD Rhapsody output | `id,rhapsody_dir` | `--modality rhapsody [--min_umi 500]` |

Ready-to-copy `bsub` commands for each case are in [bsub_examples.sh](bsub_examples.sh).

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

Steps: `generate_samplesheet_fromCR.R` lists the sample folders that have everything the workflow needs (and warns about
those that don't) -> `STAGE_CELLRANGER` links the files -> CellBender, Velocyto and DropletQC -> merge. There is no
`cellranger/` folder in the output, because Cell Ranger is not re-run. See `bsub_examples.sh` for ready-to-copy commands.
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
| `--modality`         | `rna` (Cell Ranger count), `arc` (multiome) or `atac`        | `rna`                           |
| `--skip_cellranger`  | start from existing Cell Ranger count/arc output (samplesheet `id,cellranger_dir`) | `false`    |
| `--outdir`           | Output directory root                                        | `./results`                     |
| `--selected_samples` | Optional comma-separated sample IDs to merge in final step   | all                             |
| `--min_umi`          | rhapsody only: UMI floor defining the filtered barcode set    | `500`                           |

---

## BD Rhapsody input

For data already mapped by the BD Rhapsody pipeline (skips Cell Ranger):

```csv
id,rhapsody_dir
3502749,/path/to/rhapsody/output/3502749
```
```bash
nextflow run main.nf --modality rhapsody \
--species human \
--gtf /path/to/matching/genes.gtf \
--min_umi 500 \
--samplesheet ./data/samplesheet.csv \
--outdir /path/to/output
```
`rhapsody_dir` must contain exactly one `*_Bioproduct.bam` and one `*_RSEC_MolsPerCell_Unfiltered_MEX.zip`.
A conversion step (`modules/rhapsody/rhapsody_to_cellranger.py`) builds a Cell Ranger-like `outs/` folder:

1. `raw_feature_bc_matrix.h5` from the unfiltered MEX zip
2. `filtered_feature_bc_matrix.h5` = raw matrix restricted to barcodes with at least `--min_umi` UMIs (default 500)
3. `possorted_genome_bam.bam` = the Rhapsody BAM with `UB` (from `MA`) and `RE` (exon/intron/intergenic, from `XF`) tags added

CellBender (raw h5), DropletQC (BAM + filtered h5) and Velocyto (BAM + filtered barcodes) then run as usual.
Cell IDs are the Rhapsody integer cell labels, so merged cells are named `<label>_<id>`. Notes:

- **`--min_umi`** (default 500) sets how many barcodes DropletQC and Velocyto see. DropletQC needs empty droplets to
  compare against, so this is deliberately wider than the cells Rhapsody called. Lower it to include more barcodes
  (slower), raise it to get closer to the called cells. The merged `.h5ad` keeps every barcode and stores DropletQC's
  call in `obs['cell_status']` (`cell` / `empty_droplet`) - filter on it in the downstream notebook.
- **`--gtf`**: Velocyto defaults to the Cell Ranger GENCODE GTF for `--species`; pass the annotation matching your
  Rhapsody reference if you have it.
- Sample IDs may be numeric (e.g. `3502749`). Output is in `outdir/cellbender`, `Velocyto`, `DropletQC` and `writes/`
  as for Cell Ranger input; the converted `outs/` folder stays in the Nextflow work directory.
- The BAM is read once to add the Cell Ranger-style tags, so large BAMs (20 GB+) add some minutes to the first step.
- See `bsub_examples.sh` for ready-to-copy commands.

Optional `raw_h5` / `filtered_h5` columns let you supply your own
10x-style h5 files instead (barcodes must be the Rhapsody cell labels in the BAM `CB` tag):

```csv
id,rhapsody_dir,raw_h5,filtered_h5
3502749,/path/to/rhapsody/output/3502749,/path/to/raw.h5,/path/to/filtered.h5
```

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

