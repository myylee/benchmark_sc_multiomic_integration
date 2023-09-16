#!/usr/bin/python

import sys 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))

#==== method specific ==== 
import networkx as nx
import scglue
from itertools import chain
import seaborn as sns
from matplotlib import rcParams

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

def run_glue_fn(in_dir,out_dir):
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
    
    adata_prna.obs['dataset'] = 'multiomeRNA'
    adata_patac.obs['dataset'] = 'multiomeATAC'
    adata_urna.obs['dataset'] = 'scRNA'
    adata_uatac.obs['dataset'] = 'snATAC'
    
    rna = adata_urna
    atac = adata_uatac
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir,"glue"), exist_ok=True)
    os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    
    # preprocessing of scRNA
    rna.layers["counts"] = rna.X.copy()
    sc.pp.filter_genes(rna, min_cells=3)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    sc.pp.scale(rna)
    sc.tl.pca(rna, n_comps=100, svd_solver="auto")
    
    # preprocessing of snATAC
    sc.pp.filter_genes(atac,min_counts=1)
    scglue.data.lsi(atac, n_components=100, n_iter=15)

    # build graph
    scglue.data.get_gene_annotation(
        # this works for mouse mm10 genome-build (downloaded from https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz)
        rna, gtf="/home/myylee/scmint/methods_eval/mm10_genes.gtf.gz",
        gtf_by="gene_name"
    )
    rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
    
    split = atac.var_names.str.split(r"[--]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    atac_chrs = atac.var['chrom'].value_counts().index.tolist()
    row_keep = rna.var_names[rna.var['chrom'].isin(atac_chrs).tolist()]
    rna = rna[:,row_keep].copy()
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    scglue.graph.check_graph(guidance, [rna, atac])
    
    # prepare for training 
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="counts", use_rep="X_pca"
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi"
    )

    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()
    
    # GLUE training 
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        fit_kws={"directory": os.path.join(out_dir,"glue")}
    )
    
    dx = scglue.models.integration_consistency(
        glue, {"rna": rna, "atac": atac}, guidance_hvf
    )
    print(dx)
    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    combined = ad.concat([rna, atac])

    # extract latent representation
    res_df = pd.DataFrame(combined.obsm['X_glue'],index=combined.obs.index)
    # set column names as latent_x 
    res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")
    res_df['dataset'] = combined.obs['dataset']
    res_df['dataset'] = res_df['dataset'].astype("string")
    
    # save latent representation and model
    
    csv_out = os.path.join(out_dir, "glue","glue_result.csv")
    res_df.to_csv(csv_out)
    model_out = os.path.join(out_dir,"glue","glue.dill")
    glue.save(model_out)
    stop = timeit.default_timer()
    
    print('Time(s): ', stop - start)  
    # record time 
    runtime_out = os.path.join(out_dir,"runtime","glue_runtime.txt")
    print(stop - start,  file=open(runtime_out, 'w'))
    print("------ Done ------")
    print("------ No prediction ------")

print("argument 1:",sys.argv[1])
print("argument 2:",sys.argv[2])

run_glue_fn(sys.argv[1],sys.argv[2])

