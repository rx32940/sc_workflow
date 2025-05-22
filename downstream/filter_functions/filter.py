import os
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import seaborn as sns
from scipy import stats
import numpy as np

def filter_outliers_up(number_array, up_threshold = 3):
    Zscores = stats.zscore(number_array)
    not_outliers = [(x < up_threshold) for x in Zscores.flatten()]
    upper_filter_index = np.argmin(np.abs(Zscores - up_threshold))
    print(f"Filtering out numbers { up_threshold } std up:\n Upper filter: {number_array[upper_filter_index]}")
    return not_outliers

def filter_outliers_down(number_array, down_threshold=3, hard_filter=200):
    # Calculate Z-scores
    Zscores = stats.zscore(number_array)

    # Determine the value at 3 standard deviations
    mean = np.mean(number_array)
    std_dev = np.std(number_array)
    threshold_value = mean - down_threshold * std_dev

    # Use the hard filter if the threshold value is below the hard filter
    filter_value = max(threshold_value, hard_filter)
    # Create a boolean mask to keep numbers above the filter value
    keep = [x > filter_value for x in number_array]

    # Print the filter value used
    print(f"Filtering out numbers below {filter_value:.2f}.")

    # Return the boolean mask
    return keep

def get_dropletQC_outlier(adata, sample, threshold_nf=1.5, threshold_umi=1.5, nf_cut=0.1, umi_cut=3.5, workdir='.'):
    cur_batch = adata.obs.loc[adata.obs["batch_names"] == sample, :]
    print(sample)
    nf = cur_batch["nf"]
    filtered_nf = nf
    nf_Zscore = stats.zscore(filtered_nf)
    
    if min(nf_Zscore) < (-threshold_nf):  # if the min zscore < zscore threshold
        nf_z_score_index = np.argmin(np.abs(nf_Zscore + threshold_nf))  # set a soft threshold with zscore threshold
        z_score_neg_nf_soft = filtered_nf[nf_z_score_index] # get the exact metrics of the soft threshold
        if z_score_neg_nf_soft > nf_cut: # if the metrics of the soft threshold is > hard threshold, set filter to soft threshold
            z_score_neg_nf = z_score_neg_nf_soft
        else: # else filter is hard threshold
            z_score_neg_nf = nf_cut 
    else:  # soft threshold equals to hard threshold
        z_score_neg_nf = nf_cut

    logumi = np.log10(cur_batch["umi"])
    filtered_umi = logumi
    logumi_Zscore = stats.zscore(filtered_umi)
    
    if min(logumi_Zscore) < (-threshold_umi):  # if the min zscore < zscore threshold
        logumi_z_score_index = np.argmin(np.abs(logumi_Zscore + threshold_umi))  # set a soft threshold with zscore threshold
        z_score_neg_logumi_soft = filtered_umi[logumi_z_score_index]
        if z_score_neg_logumi_soft > umi_cut: # if the metrics of the soft threshold is > hard threshold, set filter to soft threshold
            z_score_neg_logumi = z_score_neg_logumi_soft
        else: # else filter is hard threshold
            z_score_neg_logumi = umi_cut
    else:  # soft threshold equals to hard threshold
        z_score_neg_logumi = umi_cut

    # Identifying outliers
    filter_a = set(filtered_nf[filtered_nf < z_score_neg_nf].index).union(set(logumi[logumi < z_score_neg_logumi].index))
    outliers = filter_a

    cur_batch["isFiltered"] = ['red' if idx in outliers else 'blue' for idx in cur_batch.index]
    
    # Plotting
    plt.figure(figsize=(8, 6))
    plt.grid(visible=False)
    plt.scatter(cur_batch["nf"], np.log10(cur_batch["umi"]), c=cur_batch['isFiltered'], s =0.8, alpha = 0.7)
    
    # Adding dashed lines
    plt.axvline(x=z_score_neg_nf, color='gray', linestyle='--')
    plt.axhline(y=z_score_neg_logumi, color='gray', linestyle='--')
    
    # Labels and Title
    plt.xlabel('Nuclear Fraction')
    plt.ylabel('Log10(UMI)')
    plt.title(sample)
    

    # Save plot
    output_dir = os.path.join(workdir, "figures", "DropletQC_filter")
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"{sample}_DropletQC_filter.png"), bbox_inches='tight')
    plt.show()
    
    return outliers