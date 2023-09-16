#!/usr/bin/python
import sys 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))


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

#==== method specific ==== 
import scvi

scvi.settings.seed = 420
def run_multivi_fn(in_dir,out_dir):
    ##Test
    #os.makedirs(os.path.join(out_dir,"multivi"), exist_ok=True)
    #os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    #csv_out = os.path.join(out_dir, "multivi","multivi_result.csv")
    #pd.DataFrame([1,2,3,4]).to_csv(csv_out)
    
    # save latent representation and model 
    start = timeit.default_timer()
    scvi.settings.seed = 420
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
    # Horizontally stack two modalities of paired dataset 
    adata_paired = AnnData(scipy.sparse.hstack((deepcopy(adata_prna.X), 
                                            deepcopy(adata_patac.X)), 
                                           format='csr'),
                           obs = deepcopy(adata_prna.obs),
                           var = pd.concat([deepcopy(adata_prna.var[["modality"]]),deepcopy(adata_patac.var[["modality"]])]))
    # organize_mulitome_anndatats: concatenate paired and two unpaired anndata
    adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired, adata_urna, adata_uatac)
    # gene features need to be before chromatin peaks (algorithm assumption)
    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
    sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
    # setup batch annotation
    scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key='modality')
    # setup model 
    mvi = scvi.model.MULTIVI(
        adata_mvi,
        n_genes=(adata_mvi.var['modality']=='Gene Expression').sum(),
        n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
    )
    # train 
    mvi.train()
    os.makedirs(out_dir, exist_ok=True)
    # get latent representation 
    adata_mvi.obsm["MultiVI_latent"] = mvi.get_latent_representation()
   
    adata_mvi.obs = adata_mvi.obs.set_axis([s. split("_", 1)[0] for s in adata_mvi.obs.index], axis='index')

    # extract latent representation
    res_df = pd.DataFrame(adata_mvi.obsm['MultiVI_latent'],index=adata_mvi.obs.index)
    # set column names as latent_x 
    res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")
    # save modality information as dataset 
    res_df['dataset'] = adata_mvi.obs['modality']
    # convert categories to the ["snATAC","scRNA","Multiome"]
    res_df['dataset'] = res_df['dataset'].astype("category")
    res_df['dataset'].cat.categories = ["snATAC","scRNA","Multiome"]
    res_df['dataset'] = res_df['dataset'].astype("string")
    
    # save latent representation and model 
    os.makedirs(os.path.join(out_dir,"multivi"), exist_ok=True)
    os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    
    csv_out = os.path.join(out_dir, "multivi","multivi_result.csv")
    res_df.to_csv(csv_out)
    model_out = os.path.join(out_dir,"multivi","trained_multivi")
    mvi.save(model_out, overwrite=True)
    stop = timeit.default_timer()
    print('Time(s): ', stop - start)  
    # record time 
    runtime_out = os.path.join(out_dir,"runtime","multivi_runtime.txt")
    print(stop - start,  file=open(runtime_out, 'w'))
    print("------ Done ------")
    print("------ No prediction ------")
    return(mvi)

print("argument 1:",sys.argv[1])
print("argument 2:",sys.argv[2])

mvi = run_multivi_fn(sys.argv[1],sys.argv[2])
