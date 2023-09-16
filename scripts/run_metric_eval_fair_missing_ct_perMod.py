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

from numpy.random import choice
subsample_size = lambda vector, size: choice(vector, size = size, replace = False) 
    
# find the max f1 for a given cell type
# try all possible cluster corresponding to this ct 
# record TP, FP, FN for each possible cluster id
def f1_all(obs, pred, truth, ct_target):
    """cluster optimizing over largest F1 score of isolated label"""
    possible_cluster = obs.loc[obs[truth] == ct_target][pred].unique()
    col_names = ['Cluster','N','TP','FP','FN','F1']
    score = pd.DataFrame(index=range(len(possible_cluster)),columns=range(len(col_names)))
    score.columns = col_names
    for i in range(len(possible_cluster)):
        cluster = possible_cluster[i]
        print(cluster)
        obs_sel = obs.loc[(obs[truth] == ct_target) | (obs[pred]==cluster)]
        TP = sum((obs_sel[truth] == ct_target) & (obs_sel[pred] == cluster))
        FP = sum((obs_sel[truth] != ct_target) & (obs_sel[pred] == cluster)) 
        FN = sum((obs_sel[truth] == ct_target) & (obs_sel[pred] != cluster)) 
        n = obs_sel.shape[0]
        F1 = TP/(TP+0.5*(FP+FN))
        score.loc[i,:] = [cluster,n,TP,FP,FN,F1]
    score = score.sort_values('F1', ascending=False)
    return score

def latent_method_eval(out_dir,file_path,ct_ref,nclust=None):
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
    plot_name = "_"+plot_name.replace("/", "-").replace(".csv", ".png")
    
    plot_name2 = out_dir+file_path
    plot_name2 = "_"+plot_name.replace("/", "-").replace(".csv", "_with_multiome.png")
    
    idx = list(adata_plot.obs['dataset'].isin(["Multiome-RNA","Multiome-ATAC","multiomeRNA","multiomeATAC"]))
    adata_plot.obs['dataset2'] = adata_plot.obs['dataset'].tolist()
    adata_plot.obs['dataset2'][idx] = 'Multiome'
    if len(idx) > 0:
        duplicated_idx = [i for i, x in enumerate(idx) if x]
        idx_kept = list(~adata_plot.obs['dataset'].isin(["Multiome-RNA","Multiome-ATAC","multiomeRNA","multiomeATAC"]))
        idx_kept = [i for i, x in enumerate(idx_kept) if x]
        idx_f = idx_kept+ subsample_size(duplicated_idx,np.int32(np.floor(len(duplicated_idx)/2))).tolist()
        adata_plot_sel = deepcopy(adata_plot)[idx_f,:]
        adata_plot_sel.obs = deepcopy(adata_plot.obs).iloc[idx_f,:]
    else:
        adata_plot_sel = adata_plot

    if len(umap_col_sel) > 0: 
        sc.pl.umap(adata_plot_sel, 
                   color=['cell_type','dataset','predicted_ct'],
                   use_raw =False,
                   layer='umap',
                   wspace=0.4,
                   save=plot_name,
                   legend_fontsize="xx-small")
    else:
        sc.pl.umap(adata_plot_sel, 
                   color=['cell_type','dataset','predicted_ct'],
                   wspace=0.4,
                   save=plot_name,
                   legend_fontsize="xx-small")

    # if multiome dataset is used, plot all cells as well in a different UMAP plot
    if adata_plot_sel.shape[0] != adata_plot.shape[0]:
        if len(umap_col_sel) > 0: 
            sc.pl.umap(adata_plot, 
                       color=['cell_type','dataset','predicted_ct'],
                       use_raw =False,
                       layer='umap',
                       wspace=0.4,
                       save=plot_name2,
                       legend_fontsize="xx-small")
        else:
            sc.pl.umap(adata_plot, 
                       color=['cell_type','dataset','predicted_ct'],
                       wspace=0.4,
                       save=plot_name2,
                       legend_fontsize="xx-small")

    return(adata_plot_sel)


def run_metric(out_dir,file_path,ct_ref,rare_ct_path,nclust=None):
    if nclust is not None:
        nclust = int(nclust)
    adata_plot = latent_method_eval(out_dir,file_path,ct_ref,nclust)

    rare_ct_list = pd.read_csv(rare_ct_path,index_col=0).iloc[:,0].tolist()
    print("Metric calculated using this number of cells: ",str(adata_plot.shape[0]))


    pred = 'predicted_ct'
    truth = 'cell_type'
    batch = 'dataset2'
    print("Breakdown by dataset: ")
    print(adata_plot.obs[batch].value_counts())
    # calculate all metric using the same number of categories for the dataset column - 3 categories. 
    res_dir, filename  = os.path.split(os.path.join(out_dir,file_path))
    key_splits = ct_ref.split('_')
    key = key_splits[len(key_splits)-1]
    key = str(key.split('.')[0])

    obs = adata_plot.obs
    for ct_i in rare_ct_list:
        for mod in obs[batch].unique().tolist():
            obs_i = obs[obs[batch] == mod]
            ct_i_score = f1_all(obs_i,pred,truth,ct_i)
            metric_path_i = os.path.join(res_dir, "{}_{}_{}_metric.csv".format(key,ct_i,mod)) 
            ct_i_score.to_csv(metric_path_i)
    
#===== Running codes ===== 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))

# zero-based indexing 
# out_dir
print("out_dir:",sys.argv[1])
# file_path
print("file_path:",sys.argv[2])
# ct_ref
print("ct_ref:",sys.argv[3])
# rare_ct
print("rare_ct_path:",sys.argv[4])


if (len(sys.argv) == 6):
    # nclust
    print("nclust:",sys.argv[5])
    run_metric(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
else:
    run_metric(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

print("-------Finished evaluation-------")
