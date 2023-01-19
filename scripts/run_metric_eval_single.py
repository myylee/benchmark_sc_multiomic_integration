#!/usr/bin/python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import sys 
import os
import pandas as pd
import numpy as np
import pickle
import utils_eval
import matplotlib.pyplot as plt
from anndata import AnnData
import scanpy as sc
import pickle
import anndata as ad 
os.environ['R_HOME'] = '/home/myylee/anaconda3/envs/scib2/lib/R/'
import scib
from copy import deepcopy 
import tempfile
import shutil

def metrics(adata_plot, pred, truth, batch): 
    # Clustering and annotation  quality 
    # ARI
    ari = scib.metrics.ari(adata_plot,pred,truth)
    # NMI
    nmi = scib.metrics.nmi(adata_plot,pred,truth)
    
    # Cell type ASW batch (cell types are separated)
    ct_asw = scib.metrics.silhouette(adata_plot,truth,"embed")
    
    # LISI
    lisi_res = scib.metrics.lisi.lisi_graph(adata_plot,batch,truth,type_="embed")
    # cell type - cLSi (scaled to [0,1])
    #clisi = scib.metrics.lisi.lisi_graph(adata_plot,batch,truth)[1]
    clisi = lisi_res[1]
    
    # Batch removal 
    # batch iLSi (scaled to [0,1]) [did not consider the size imbalance between batches]
    ilisi = lisi_res[0]

    # ASW - batch: batches are well mixed for each cell type
    b_asw = scib.metrics.silhouette_batch(adata_plot, batch, truth,"embed")

    #kBET (type_="knn"[for graph-based integration], type_=None [default, for joint_embed, or feature matrix])
    kbet = scib.metrics.kBET(adata_plot,batch,truth,embed="embed")

    # graph connectivity
    gconn = scib.metrics.graph_connectivity(adata_plot, truth)
    
    # isolated cluster
    # isolated labels - F1
    iso_f1 = scib.metrics.isolated_labels(adata_plot, truth,batch,"embed",iso_threshold=2)
    
    # isolated labels - ASW
    iso_asw = scib.metrics.isolated_labels(adata_plot, truth,batch,"embed",cluster=False,iso_threshold=2)

    metrics = [ari, nmi, ct_asw, clisi, b_asw, ilisi, kbet, gconn, iso_f1, iso_asw]
    
    #metrics_tbl = tabulate(metrics, headers=["ARI", "NMI", "ct_aws","clisi","b_saw","ilisi","kbet","gconn","iso_f1","iso_asw"],floatfmt=".4f")
    df = pd.DataFrame(metrics).T
    df = df.set_axis(["ARI", "NMI", "ct_aws","clisi","b_saw","ilisi","kbet","gconn","iso_f1","iso_asw"],axis=1)
    
    return(df)
    
def latent_method_eval(out_dir,file_path,ct_ref,nclust=None):
    from os.path import exists
    # having issue with file being deleted when it is accessed with read_csv function. Thus, copy the result file into working directory with a temporary filename, read it, and then delete the temporary file. 
    path = os.path.join(out_dir, file_path)
    res_df = pd.read_csv(path,index_col=0)
    
    col_sel = [col for col in res_df.columns if 'latent' in col]
    sub_df = res_df.loc[: , col_sel]

    cell_type = pd.read_csv(ct_ref,index_col=0)
    cell_type.columns = ['cell_type']
    # Assuming everything before '_' is not useful 
    cell_bc = [sub.split('_')[len(sub.split('_'))-1] for sub in list(res_df.index)]
    res_df['cell_type'] = list(cell_type.loc[cell_bc,'cell_type'])

    adata_plot = AnnData(res_df.loc[: , col_sel].to_numpy(),
                         obs=res_df.loc[:,['dataset','cell_type']],
                         dtype=np.float32)
    sc.pp.neighbors(adata_plot, n_neighbors=15, use_rep="X")
    umap_col_sel = [col for col in res_df.columns if 'umap' in col]
    if len(umap_col_sel) > 0: 
        adata_plot.obsm['umap'] = res_df.loc[: , umap_col_sel].to_numpy()
    else:
        sc.tl.umap(adata_plot, min_dist=0.2)

    if nclust is not None:
        resolution = utils_eval.find_resolution_louvain(adata_plot,nclust)
        sc.tl.louvain(adata_plot, resolution = resolution, random_state = 0)
        adata_plot.obs['predicted_ct'] = adata_plot.obs['louvain']
    else:
        adata_plot.obs['predicted_ct'] = res_df.loc[:,['predicted_ct']]
        adata_plot.obs['predicted_ct'] = adata_plot.obs['predicted_ct'].astype('category')
    
    adata_plot.obsm['embed'] = adata_plot.X
    adata_plot.obsm['X_emb'] = adata_plot.X
    
    plot_name = out_dir+file_path
    plot_name = "_"+plot_name.replace("/", "-").replace(".csv", ".svg")
    
    plot_name2 = out_dir+file_path
    plot_name2 = "_"+plot_name.replace("/", "-").replace(".csv", "_with_multiome.svg")
    
    idx = list(adata_plot.obs['dataset'].isin(["scRNA","snATAC"]))
    adata_plot_sel = deepcopy(adata_plot)[idx,:]
    adata_plot_sel.obs = deepcopy(adata_plot.obs).loc[idx,:]
    
    if len(umap_col_sel) > 0: 
        sc.pl.umap(adata_plot_sel, 
                   color=['cell_type','dataset','predicted_ct'],
                   use_raw =False,
                   layer='umap',
                   wspace=0.4,
                   save=plot_name,
                   legend_fontsize="x-small")
    else:
        sc.pl.umap(adata_plot_sel, 
                   color=['cell_type','dataset','predicted_ct'],
                   wspace=0.4,
                   save=plot_name,
                   legend_fontsize="x-small")

    # if multiome dataset is used, plot all cells as well in a different UMAP plot
    if adata_plot_sel.shape[0] != adata_plot.shape[0]:
        if len(umap_col_sel) > 0: 
            sc.pl.umap(adata_plot, 
                       color=['cell_type','dataset','predicted_ct'],
                       use_raw =False,
                       layer='umap',
                       wspace=0.4,
                       save=plot_name2,
                   legend_fontsize="x-small")
        else:
            sc.pl.umap(adata_plot, 
                       color=['cell_type','dataset','predicted_ct'],
                       wspace=0.4,
                       save=plot_name2,
                   legend_fontsize="x-small")

    return(adata_plot_sel)

def run_metric(out_dir,file_path,ct_ref,nclust=None):
    if nclust is not None:
        nclust = int(nclust)
    adata_plot = latent_method_eval(out_dir,file_path,ct_ref,nclust)
    
    pred = 'predicted_ct'
    truth = 'cell_type'
    batch = 'dataset'
    
    print("Metric calculated using this number of cells: ",str(adata_plot.shape[0]))
    res = metrics(adata_plot,pred,truth,batch)
    #res['run'] = keyword # should contain cond_key, split, repeat, method 

    #os.makedirs(out_dir+"results/", exist_ok=True)
    res_dir, filename  = os.path.split(os.path.join(out_dir,file_path))
    key_splits = ct_ref.split('_')
    key = key_splits[len(key_splits)-1]
    key = str(key.split('.')[0])
    print(key)
    metric_path = os.path.join(res_dir, "{}_metric.csv".format(key)) 
    print(metric_path)
    res.to_csv(metric_path)

#===== Running codes ===== 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))

print("argument 1:",sys.argv[1])
print("argument 2:",sys.argv[2])
print("argument 3:",sys.argv[3])


if (len(sys.argv) == 5):
    print("argument 4:",sys.argv[4])
    cell_type = pd.read_csv(sys.argv[3],index_col=0)
    print(cell_type)
    run_metric(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
else:
    run_metric(sys.argv[1],sys.argv[2],sys.argv[3])

print("-------Finished evaluation-------")

