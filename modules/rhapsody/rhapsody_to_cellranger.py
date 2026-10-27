#!/usr/bin/env python
"""
Convert BD Rhapsody pipeline output into a Cell Ranger-like `outs/` directory so the
existing CellBender / DropletQC / Velocyto steps can run unchanged.

Steps:
  1. raw h5       <- Rhapsody *_RSEC_MolsPerCell_Unfiltered_MEX.zip (or a 10x-style .h5 given to --unfiltered)
  2. filtered h5  <- raw h5 restricted to barcodes with >= --min_umi total UMIs
                     (or a 10x-style .h5 given to --filtered, used as is)
  3. BAM          <- Rhapsody *_Bioproduct.bam with Cell Ranger-style tags added (see below)

Produces under <outdir>/outs:
    raw_feature_bc_matrix.h5
    filtered_feature_bc_matrix.h5
    filtered_feature_bc_matrix/{barcodes,features}.tsv.gz
    possorted_genome_bam.bam (+ .bai)

BAM tags added (all existing tags are kept; reads are never dropped):
    UB  <- MA   adjusted UMI (Cell Ranger name; needed by velocyto). Reads without MA get no UB.
    RE  <- XF   region type (Cell Ranger name; needed by DropletQC), uniquely mapped reads only:
                  gene hit (no "__" prefix)  -> E (exon)
                  __introns_                 -> N (intron)
                  __intergenic               -> I (intergenic)
                other XF classes (ambiguous / antisense) get no RE tag.
Cell barcodes (CB) are the Rhapsody integer cell labels, identical in the BAM and the MEX barcodes.
"""
import argparse
import gzip
import os
import shutil
import sys
import tempfile
import zipfile

import h5py
import numpy as np
import pandas as pd
import pysam
from scipy import io as spio
from scipy import sparse


def read_mex(zip_path):
    """Read a Rhapsody MEX zip; return (csc matrix features x barcodes, barcodes, features df)."""
    tmp = tempfile.mkdtemp(prefix="mex_")
    try:
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(tmp)
        mat = spio.mmread(os.path.join(tmp, "matrix.mtx.gz")).tocsc()
        barcodes = pd.read_csv(os.path.join(tmp, "barcodes.tsv.gz"), header=None, sep="\t", dtype=str)[0].tolist()
        features = pd.read_csv(os.path.join(tmp, "features.tsv.gz"), header=None, sep="\t", dtype=str)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if mat.shape != (len(features), len(barcodes)):
        sys.exit(f"MEX shape {mat.shape} does not match features ({len(features)}) x barcodes ({len(barcodes)}) in {zip_path}")
    return mat, barcodes, features


def read_10x_h5(path):
    """Read a 10x-style h5; return (csc matrix features x barcodes, barcodes, features df)."""
    with h5py.File(path, "r") as f:
        if "matrix" not in f or "barcodes" not in f["matrix"]:
            sys.exit(f"{path} is not a 10x-style h5 (no matrix/barcodes)")
        g = f["matrix"]
        barcodes = [b.decode() for b in g["barcodes"][:]]
        shape = tuple(int(x) for x in g["shape"][:])
        mat = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape)
        ft = g["features"]
        cols = [[x.decode() for x in ft[k][:]] for k in ("id", "name", "feature_type")]
        features = pd.DataFrame(dict(enumerate(cols)))
    return mat, barcodes, features


def load_matrix(path):
    return read_mex(path) if path.endswith(".zip") else read_10x_h5(path)


def write_10x_h5(path, mat, barcodes, features, genome):
    """Write a 10x v3-style h5 (cells as columns of a CSC matrix)."""
    mat = sparse.csc_matrix(mat)
    mat.sort_indices()
    n = len(features)
    ids = features[0].astype(str).values
    names = features[1].astype(str).values if features.shape[1] > 1 else ids
    types = features[2].astype(str).values if features.shape[1] > 2 else np.array(["Gene Expression"] * n)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        g.create_dataset("data", data=mat.data.astype(np.int32))
        g.create_dataset("indices", data=mat.indices.astype(np.int64))
        g.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array(mat.shape, dtype=np.int32))
        ft = g.create_group("features")
        ft.create_dataset("id", data=np.array(ids, dtype="S"))
        ft.create_dataset("name", data=np.array(names, dtype="S"))
        ft.create_dataset("feature_type", data=np.array(types, dtype="S"))
        ft.create_dataset("genome", data=np.array([genome] * n, dtype="S"))
        ft.create_dataset("_all_tag_keys", data=np.array(["genome"], dtype="S"))


def write_filtered_tables(filt_dir, barcodes, features):
    """barcodes.tsv.gz is read by DropletQC and velocyto run10x."""
    os.makedirs(filt_dir, exist_ok=True)
    with gzip.open(os.path.join(filt_dir, "barcodes.tsv.gz"), "wt") as f:
        f.write("\n".join(barcodes) + "\n")
    features.to_csv(os.path.join(filt_dir, "features.tsv.gz"), sep="\t", header=False, index=False)


def region_tag(xf):
    if not xf.startswith("__"):
        return "E"
    if xf.startswith("__introns_"):
        return "N"
    if xf == "__intergenic":
        return "I"
    return None


def retag_bam(bam_in, bam_out, threads):
    n = n_ub = n_re = 0
    with pysam.AlignmentFile(bam_in, "rb", threads=threads) as fin, \
         pysam.AlignmentFile(bam_out, "wb", template=fin, threads=threads) as fout:
        for read in fin:
            n += 1
            if read.has_tag("MA"):
                read.set_tag("UB", read.get_tag("MA"), "Z")
                n_ub += 1
            if (not read.is_unmapped and not read.is_secondary and read.has_tag("XF")
                    and read.has_tag("NH") and read.get_tag("NH") == 1):
                re_val = region_tag(read.get_tag("XF"))
                if re_val:
                    read.set_tag("RE", re_val, "A")
                    n_re += 1
            fout.write(read)
    pysam.index("-@", str(threads), bam_out)
    print(f"BAM reads: {n}; UB added: {n_ub}; RE added: {n_re}", file=sys.stderr)
    if n_ub == 0:
        sys.exit("No reads had an MA tag; cannot create UB tags. Is this a Rhapsody *_Bioproduct.bam?")
    if n_re == 0:
        sys.exit("No reads had a usable XF tag; cannot create RE tags. Is this a Rhapsody *_Bioproduct.bam?")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="Rhapsody *_Bioproduct.bam")
    p.add_argument("--unfiltered", required=True,
                   help="all barcodes: *_RSEC_MolsPerCell_Unfiltered_MEX.zip or a 10x-style .h5")
    p.add_argument("--filtered", default=None,
                   help="optional 10x-style .h5 of called cells, used instead of filtering by --min_umi")
    p.add_argument("--min_umi", type=int, default=500,
                   help="keep barcodes with at least this many UMIs in the filtered matrix [default = %(default)s]")
    p.add_argument("--outdir", required=True, help="Output sample directory; files go in <outdir>/outs")
    p.add_argument("--genome", default="unknown", help="Genome label written into the h5 feature table")
    p.add_argument("--threads", type=int, default=4)
    a = p.parse_args()

    outs = os.path.join(a.outdir, "outs")
    filt_dir = os.path.join(outs, "filtered_feature_bc_matrix")
    os.makedirs(filt_dir, exist_ok=True)

    raw, raw_bcs, feats = load_matrix(a.unfiltered)
    write_10x_h5(os.path.join(outs, "raw_feature_bc_matrix.h5"), raw, raw_bcs, feats, a.genome)
    print(f"raw: {raw.shape[0]} features x {raw.shape[1]} barcodes", file=sys.stderr)

    if a.filtered:
        filt, filt_bcs, filt_feats = read_10x_h5(a.filtered)
        print(f"filtered: {filt.shape[1]} barcodes (user-supplied h5)", file=sys.stderr)
    else:
        umi = np.asarray(raw.sum(axis=0)).ravel()
        keep = np.where(umi >= a.min_umi)[0]
        if len(keep) == 0:
            sys.exit(f"No barcodes have >= {a.min_umi} UMIs; lower --min_umi")
        filt = raw[:, keep]
        filt_bcs = [raw_bcs[i] for i in keep]
        filt_feats = feats
        print(f"filtered: {len(keep)} of {raw.shape[1]} barcodes with >= {a.min_umi} UMIs", file=sys.stderr)
    write_10x_h5(os.path.join(outs, "filtered_feature_bc_matrix.h5"), filt, filt_bcs, filt_feats, a.genome)
    write_filtered_tables(filt_dir, filt_bcs, filt_feats)

    retag_bam(a.bam, os.path.join(outs, "possorted_genome_bam.bam"), a.threads)


if __name__ == "__main__":
    main()
