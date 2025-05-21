#!/usr/bin/env python

import os
import argparse
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
from cellbender.remove_background.downstream import anndata_from_h5
import anndata as ad
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser(description="Combine and filter batches")
    parser.add_argument("--workdir", required=True, help="Working directory")
    parser.add_argument("--samplesheet", required=True, help="Path to samplesheet.csv")
    parser.add_argument("--output", required=True, help="Output AnnData file path")
    parser.add_argument("--selected_batches", required=False, default=None,
                        help="Path to a file containing a list of selected batches (one per line) or a comma-separated list. Defaults to all batches in the samplesheet.")
    args = parser.parse_args()

    workdir = args.workdir
    samplesheet_file = args.samplesheet
    output_file = args.output
    selected_batches_file = args.selected_batches

    # Read sample IDs from samplesheet.csv
    samplesheet = pd.read_csv(samplesheet_file)
    all_batches = set(samplesheet['id'].tolist())

    # Handle selected batches
    if selected_batches_file:
        # Check if this is a file or a comma-separated list
        if os.path.isfile(selected_batches_file):
            with open(selected_batches_file, "r") as f:
                selected_batches = set(line.strip() for line in f.readlines())
        else:
            selected_batches = set([batch.strip() for batch in selected_batches_file.split(",")])

        # Filter to only include batches present in the samplesheet
        batches = list(all_batches.intersection(selected_batches))
        print(f"Selected {len(batches)} batches from the provided list.")
    else:
        batches = list(all_batches)
        print(f"Using all {len(batches)} batches from the samplesheet.")

    # Sort batches for consistent ordering
    batches.sort()

    batch_dict = {}
    for batch in tqdm(batches, desc="Processing batches"):
        print(f"Processing batch: {batch}")

        # Define file paths
        cellbender_file = os.path.join(workdir, "cellbender", batch + '_cb', "cellbender_output.h5")
        dropletqc_file = os.path.join(workdir,  "DropletQC", batch, "nf_ed_qc.csv")
        velocity_file = os.path.join(workdir,  "Velocyto", batch + ".loom")
        
        # Check if all required files exist
        if not (os.path.exists(cellbender_file) and os.path.exists(dropletqc_file) and os.path.exists(velocity_file)):
            print(f"Skipping batch {batch} due to missing files.")
            continue

        # Load CellBender data
        cellbender_filtered = anndata_from_h5(cellbender_file)

        # Load DropletQC data
        dropletQC_df = pd.read_csv(dropletqc_file, index_col=0)
        cellbender_filtered.obs = cellbender_filtered.obs.merge(dropletQC_df, left_index=True, right_index=True, how="left")

        # Filter cells
        cellbender_filtered = cellbender_filtered[~pd.isna(cellbender_filtered.obs['nf'])].copy()

        # Load RNA velocity data
        vadata = sc.read(velocity_file)
        cellbender_filtered = scv.utils.merge(cellbender_filtered, vadata)

        # Add batch information to cell names
        cellbender_filtered.obs.index = cellbender_filtered.obs.index + "_" + batch
        cellbender_filtered.var_names_make_unique()
        batch_dict[batch] = cellbender_filtered

    # Merge all batches
    merged_adata = ad.concat(batch_dict, join='outer', label="batch_names")
    merged_adata.obs_names_make_unique()
    merged_adata.var_names_make_unique()

    # Save the merged AnnData object
    merged_adata.write(output_file)
    print(f"Saved merged AnnData object to {output_file}")

if __name__ == "__main__":
    main()
