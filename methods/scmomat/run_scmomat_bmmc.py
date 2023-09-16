#!/usr/bin/python

import sys 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))

#==== method specific ==== 
import scmomat
import torch

#from matplotlib import rcParams
from anndata import AnnData
import anndata as ad
import scipy
import numpy as np
import pandas as pd
import scipy.io as sio
import os
import scanpy as sc
from copy import deepcopy
from utils_eval import read_mtx_folder, write_adata
import timeit

# ------ Helper Functions -----
# copied from https://github.com/PeterZZQ/scMoMaT/blob/main/scmomat/utils.py
# slight adjustment to preprocess() line 91 of this document. An error occured due to third dimension introduced in the previous calculations. 
def quantile_norm(X):
    from scipy import stats
    """Normalize the columns of X to each have the same distribution.
    Given an expression matrix (microarray data, read counts, etc) of M genes
    by N samples, quantile normalization ensures all samples have the same
    spread of data (by construction).
    The data across each row are averaged to obtain an average column. Each
    column quantile is replaced with the corresponding quantile of the average
    column.
    Parameters
    ----------
    X : 2D array of float, shape (M, N)
        The input data, with M rows (genes/features) and N columns (samples).
    Returns
    -------
    Xn : 2D array of float, shape (M, N)
        The normalized data.
    """
    # compute the quantiles
    quantiles = np.mean(np.sort(X, axis=0), axis=1)

    # compute the column-wise ranks. Each observation is replaced with its
    # rank in that column: the smallest observation is replaced by 1, the
    # second-smallest by 2, ..., and the largest by M, the number of rows.
    ranks = np.apply_along_axis(stats.rankdata, 0, X)

    # convert ranks to integer indices from 0 to M-1
    rank_indices = ranks.astype(int) - 1

    # index the quantiles for each rank with the ranks matrix
    Xn = quantiles[rank_indices]

    return(Xn)

def quantile_norm_log(X, log = True):
    if log:
        logX = np.log1p(X)
    else:
        logX = X
    logXn = quantile_norm(logX)
    return logXn


def preprocess(counts, modality = "RNA", log = True):
    """\
    Description:
    ------------
        Preprocess the dataset, for count, interaction matrices
    
    Parameters:
    -------------
        counts: count matrix
        modality: modality that the matrix belong to, can be ``ATAC'', ``RNA'', ``protein''
        log: log transform the data or not
    """
    if modality == "ATAC":
        # make binary, maximum is 1
        counts = (counts > 0).astype(np.float) 

    else:
        # other cases, e.g. Protein, RNA, etc
        counts = quantile_norm_log(counts, log = log)[:,:,0]
        counts = counts/np.max(counts)

    return counts
# ------ Helper Functions End -----
def run_scmomat_fn(in_dir,out_dir):
    start = timeit.default_timer()
    
    adata_prna = read_mtx_folder(os.path.join(in_dir,"paired_RNA/"),
                                       "Gene Expression",
                                       ["gene"],
                                       ["barcodes"])

    adata_patac = read_mtx_folder(os.path.join(in_dir,"paired_ATAC/"),
                                       "Peaks",
                                       ["peak"],
                                       ["barcodes"])

    adata_urna = read_mtx_folder(os.path.join(in_dir,"unpaired_RNA/"),
                                       "Gene Expression",
                                       ["gene"],
                                       ["barcodes"])

    adata_uatac = read_mtx_folder(os.path.join(in_dir,"unpaired_ATAC/"),
                                       "Peaks",
                                       ["peak"],
                                       ["barcodes"])
    
    adata_prna.obs['dataset'] = 'multiome'
    adata_patac.obs['dataset'] = 'multiome'
    adata_urna.obs['dataset'] = 'scRNA'
    adata_uatac.obs['dataset'] = 'snATAC'
    
    # find the hvg for RNA modality
    adata_prna.layers["counts"] = adata_prna.X.copy()
    sc.pp.normalize_total(adata_prna, target_sum=1e4)
    sc.pp.log1p(adata_prna)
    sc.pp.highly_variable_genes(adata_prna,n_top_genes=5000)
    hvg_prna = adata_prna.var['highly_variable'].index[adata_prna.var['highly_variable']].tolist()

    adata_urna.layers["counts"] = adata_urna.X.copy()
    sc.pp.normalize_total(adata_urna, target_sum=1e4)
    sc.pp.log1p(adata_urna)
    sc.pp.highly_variable_genes(adata_urna,n_top_genes=5000)
    hvg_urna = adata_urna.var['highly_variable'].index[adata_urna.var['highly_variable']].tolist()

    hvg_rna_common = np.intersect1d(hvg_prna,hvg_urna)
    
    ## calculate the conversion between ATAC to pseudoRNA counts
    A_long = pd.read_csv('/home/myylee/scmint/methods_eval/dataset/bmmc/bmmc_peak_GxR.csv')
    A_long_sel = A_long[A_long.loc[:, 'gene.name'].isin(hvg_rna_common.tolist())]
    A_sel = pd.crosstab(A_long_sel.loc[:, 'peak'], A_long_sel.loc[:, 'gene.name'])

    # filter for peaks that are within 2000bp of TSS or along the gene body, as well as genes that have at least one peak
    gene_final = A_sel.columns.values.squeeze()
    peak_final = A_sel.index.values.squeeze()

    counts_prna = adata_prna[:,gene_final].layers['counts'].todense()
    counts_urna = adata_urna[:,gene_final].layers['counts'].todense()

    counts_prna = preprocess(counts_prna, modality = "RNA", log = False) 
    counts_urna = preprocess(counts_urna, modality = "RNA", log = False) 


    counts_patac = adata_patac[:,peak_final].X.todense()
    counts_uatac = adata_uatac[:,peak_final].X.todense()

    counts_patac = preprocess(counts_patac, modality = "ATAC")
    counts_uatac = preprocess(counts_uatac, modality = "ATAC")

    counts_uatac_pseudoRNA = counts_uatac @ A_sel
    #BINARIZE, still is able to see the cluster pattern, much denser than scRNA-Seq (cluster pattern clearer)
    counts_uatac_pseudoRNA = (counts_uatac_pseudoRNA!=0).astype(int)
    counts_uatac_pseudoRNA = np.asmatrix(counts_uatac_pseudoRNA.to_numpy())
    
    counts = {"rna":[counts_prna,counts_urna,counts_uatac_pseudoRNA], "atac": [counts_patac,None,counts_uatac] }
    
    feats_name = {"rna": gene_final, "atac": peak_final}
    counts["feats_name"] = feats_name
    # might change this is there are 2 batches or more in either data types.
    counts["nbatches"] = 3

    # Running scMoMaT model
    device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
    lamb = 0.001
    batchsize = 0.1
    # running seed
    seed = 0
    # number of latent dimensions
    K = 30
    interval = 1000
    T = 4000
    lr = 1e-2
    # 1st stage training, learning cell factors
    model = scmomat.scmomat_model(counts = counts, K = K, device = device)
    losses = model.train_func(T = T)

    # extract cell factors/latent representations
    
    zs = model.extract_cell_factors()
    
    cell_bc = np.concatenate([adata_prna.obs.index.values.squeeze(),
                          adata_urna.obs.index.values.squeeze(),
                          adata_uatac.obs.index.values.squeeze()])
    
    res_df = pd.DataFrame(np.concatenate(zs, axis=0),
                          index=cell_bc)
    # set column names as latent_x 
    res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")
    res_df['dataset'] = np.concatenate([adata_prna.obs['dataset'].values.squeeze(),
                                        adata_urna.obs['dataset'].values.squeeze(),
                                        adata_uatac.obs['dataset'].values.squeeze()])
    res_df['dataset'] = res_df['dataset'].astype("string")
    
    # save latent representation and model
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir,"scmomat"), exist_ok=True)
    os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    
    csv_out = os.path.join(out_dir, "scmomat","scmomat_result.csv")
    res_df.to_csv(csv_out)
    stop = timeit.default_timer()
    
    print('Time(s): ', stop - start)  
    # record time 
    runtime_out = os.path.join(out_dir,"runtime","scmomat_runtime.txt")
    print(stop - start,  file=open(runtime_out, 'w'))
    print("------ Done ------")
    
    print("------ No Prediction ------")

print("argument 1:",sys.argv[1])
print("argument 2:",sys.argv[2])

run_scmomat_fn(sys.argv[1],sys.argv[2])

