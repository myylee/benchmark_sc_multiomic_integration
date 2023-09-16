import numpy as np 
import pandas as pd 
from anndata import AnnData
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import svds
from sklearn.metrics import mean_squared_error
import os 
from copy import deepcopy
import sklearn.metrics as metrics
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi
import scipy.io as sio


def read_mtx_folder(dir_path,mod_key,var_list,obs_list):
    mtx_path = []
    tsv_path = []
    for file in os.listdir(dir_path):
        if file.endswith(".mtx"):
            mtx_path.append(os.path.join(dir_path, file))
    if(len(mtx_path)==0):
        print("no .mtx file found, exiting function")
        return
    counts = sio.mmread(mtx_path[0])
    if type(counts) != "scipy.sparse.csr.csr_matrix": counts = counts.tocsr()
    
    var_v = []
    for i in var_list:
        var_v.append(pd.read_csv(os.path.join(dir_path,i+".tsv"),header=None))
    var_df = pd.concat(var_v, axis=1)
    var_df = var_df.set_axis(var_list, axis='columns')
    obs_v = []
    for i in obs_list:
        obs_v.append(pd.read_csv(os.path.join(dir_path,i+".tsv"),header=None))
    obs_df = pd.concat(obs_v, axis=1)
    obs_df = obs_df.set_axis(obs_list, axis='columns')
    if(counts.shape[0]!=obs_df.shape[0]): counts=deepcopy(counts.transpose())
    adata = AnnData(counts,obs=obs_df,var=var_df)
    adata.obs = adata.obs.set_axis(list(adata.obs["barcodes"]),axis="index")
    adata.var= adata.var.set_axis(list(adata.var[var_list[0]]),axis="index")
    adata.var["modality"] = mod_key
    return(adata)
        

def one_sample_exclusive_split(adata_rna,adata_atac,hold_ct,md_col,u1=True,r=1/3,transpose=False):
    from numpy.random import choice
    from numpy import setdiff1d
    from math import floor
    
    # randomly sample from 'vector' with the length defined by 'frac'
    subsample = lambda vector, frac: choice(vector, size = round(len(vector) * frac), replace = False) 

    if transpose:
        adata_rna = adata_rna.transpose()
        adata_atac = adata_atac.transpose()
    assert(adata_rna.shape[0] == adata_atac.shape[0])
    assert(md_col in adata_rna.obs.columns)
    ncells = adata_rna.shape[0]

    cells_excluded = np.nonzero(np.array(~adata_rna.obs[md_col].isin([hold_ct])))[0]
    r_adj_p = ncells*r / len(cells_excluded)
    paired_idx  = subsample(cells_excluded,r_adj_p)
    
    unpaired_idx = setdiff1d(range(ncells), paired_idx)
    unpaired_idx_excluded = setdiff1d(cells_excluded, paired_idx)
    
    r_adj_u = 0.5*len(unpaired_idx) / len(unpaired_idx_excluded)
    if u1:
        u2_idx = subsample(unpaired_idx_excluded,r_adj_u)
        u1_idx = setdiff1d(unpaired_idx, u2_idx)
    else: 
        u1_idx = subsample(unpaired_idx_excluded,r_adj_u)
        u2_idx = setdiff1d(unpaired_idx, u1_idx)

    adata_p1 = adata_rna[paired_idx,:]
    adata_p2 = adata_atac[paired_idx,:]
    adata_u1 = adata_rna[u1_idx,:]
    adata_u2 = adata_atac[u2_idx,:]
    
    if transpose:
        adata_p1 = adata_p1.copy().transpose()
        adata_p2 = adata_p2.copy().transpose()
        adata_u1 = adata_u1.copy().transpose()
        adata_u2 = adata_u2.copy().transpose()
    
    return((adata_p1,adata_p2,adata_u1,adata_u2))



def unpair_exclusive_split(adata_rna,adata_atac,hold_ct,md_col,r=1/3,transpose=False):
    from numpy.random import choice
    from numpy import setdiff1d
    from math import floor
    
    # randomly sample from 'vector' with the length defined by 'frac'
    subsample = lambda vector, frac: choice(vector, size = round(len(vector) * frac), replace = False) 

    if transpose:
        adata_rna = adata_rna.transpose()
        adata_atac = adata_atac.transpose()
    assert(adata_rna.shape[0] == adata_atac.shape[0])
    assert(md_col in adata_rna.obs.columns)
    ncells = adata_rna.shape[0]

    cells_excluded = np.nonzero(np.array(~adata_rna.obs[md_col].isin([hold_ct])))[0]
    r_adj = ncells*r / len(cells_excluded)
    paired_idx  = subsample(cells_excluded,r_adj)
    unpaired_idx = setdiff1d(range(ncells), paired_idx)
    u1_idx = subsample(unpaired_idx,1/2)
    u2_idx = setdiff1d(unpaired_idx, u1_idx)

    adata_p1 = adata_rna[paired_idx,:]
    adata_p2 = adata_atac[paired_idx,:]
    adata_u1 = adata_rna[u1_idx,:]
    adata_u2 = adata_atac[u2_idx,:]
    
    if transpose:
        adata_p1 = adata_p1.copy().transpose()
        adata_p2 = adata_p2.copy().transpose()
        adata_u1 = adata_u1.copy().transpose()
        adata_u2 = adata_u2.copy().transpose()
    
    return((adata_p1,adata_p2,adata_u1,adata_u2))



def pair_unpair_split(adata_rna,adata_atac,r=1/3,transpose=False):
    from numpy.random import choice
    from numpy import setdiff1d
    from math import floor
    
    # randomly sample from 'vector' with the length defined by 'frac'
    subsample = lambda vector, frac: choice(vector, size = round(len(vector) * frac), replace = False) 

    if transpose:
        adata_rna = adata_rna.transpose()
        adata_atac = adata_atac.transpose()
    assert(adata_rna.shape[0] == adata_atac.shape[0])
    ncells = adata_rna.shape[0]
    paired_idx  = subsample(range(ncells),r)
    unpaired_idx = setdiff1d(range(ncells), paired_idx)
    u1_idx = subsample(unpaired_idx,1/2)
    u2_idx = setdiff1d(unpaired_idx, u1_idx)

    adata_p1 = adata_rna[paired_idx,:]
    adata_p2 = adata_atac[paired_idx,:]
    adata_u1 = adata_rna[u1_idx,:]
    adata_u2 = adata_atac[u2_idx,:]
    
    if transpose:
        adata_p1 = adata_p1.copy().transpose()
        adata_p2 = adata_p2.copy().transpose()
        adata_u1 = adata_u1.copy().transpose()
        adata_u2 = adata_u2.copy().transpose()
    
    return((adata_p1,adata_p2,adata_u1,adata_u2))

def write_adata(adata, folder_path,modality,feature,bc="rna.bc",feature_name=None,transpose=False,additional_obs=None):
    import os
    os.makedirs(folder_path, exist_ok=True)
    if transpose:
        adata = adata.copy().transpose()
    from scipy.sparse import coo_matrix

    sio.mmwrite(os.path.join(folder_path,modality+"_counts.mtx"),coo_matrix(adata.X))
    #adata.var[celltype].to_csv(folder_path+"cell_type.tsv",header=False,index=False)
    adata.var[bc].to_csv(os.path.join(folder_path,"barcodes.tsv"),header=False,index=False)
    if feature_name is not None:
        adata.obs[feature_name].to_csv(os.path.join(folder_path,feature+".tsv"),header=False,index=False)
    else:
        adata.obs.iloc[:,0].to_csv(os.path.join(folder_path,feature+".tsv"),header=False,index=False)
    if additional_obs is not None:
        for obs_i in additional_obs:
            adata.var[obs_i].to_csv("{}/{}.tsv".format(folder_path,obs_i),header=False,index=False)
        

def find_resolution_leiden(adata_, n_clusters, random = 0): 
    import scanpy as sc
    adata = adata_.copy()
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1000.]
    
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        sc.tl.leiden(adata, resolution = current_res, random_state = random)
        labels = adata.obs['leiden']
        obtained_clusters = len(np.unique(labels))
        
        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res
        
        iteration = iteration + 1
        
    return current_res

def find_resolution_louvain(adata_, n_clusters, random = 0): 
    import scanpy as sc
    adata = adata_.copy()
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1000.]
    
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        sc.tl.louvain(adata, resolution = current_res, random_state = random)
        labels = adata.obs['louvain']
        obtained_clusters = len(np.unique(labels))
        
        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res
        
        iteration = iteration + 1
        
    return current_res

def calculate_lisi(umap_df,metadata_df,metadata_keys): 
    """A function to compute lisi, batch-mixing metric"""
    # uses r-lisi package. Install all necessary packages if already. Assumes R is installed.
    
    # ===== load rpy2 =====#
    try:
        __import__("rpy2")
    except ImportError:
        import os
        os.system('pip install rpy2')      
        
    # function: ensure R package is installed
        # default - download from CRAN 
        # if devtools_path present, then use devtools to download from a github
        # if contriburl present, download from that url
    def importr_tryhard(packname, devtools_path=None,contriburl=None):
        #import rpy2.rinterface as rpy2_int
        from rpy2.robjects.packages import PackageNotInstalledError
        try:
            rpack = importr(packname)
        except PackageNotInstalledError:
            utils = importr('utils')
            if devtools_path is not None:
                devtools_r = importr_tryhard("devtools")
                devtools_r.install_github(devtools_path)
            elif contriburl is not None:
                utils.install_packages(packname, contriburl = contriburl)
            else: 
                utils.chooseCRANmirror(ind=1)
                utils.install_packages(packname)
            # try loading again after it should've begin installed. Print message if still cannot be loaded correctly.
            try:
                rpack = importr(packname)
            except PackageNotInstalledError:
                print(packname+" cannot be installed from correctly. Either the name is not found in CRAN, or if no input for devtools::install_github or url option")
                return
        return rpack


    # ===== load necessary functions=====##
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    from rpy2.robjects.vectors import StrVector
    from rpy2.robjects.packages import importr
    utils = importr('utils')

    # loads lisi library, install using devtools from the github path if lisi is not installed
    lisi = importr_tryhard("lisi",devtools_path="immunogenomics/lisi")
    
    # convert df to r-object
    r_df1 = ro.conversion.py2rpy(umap_df)
    r_df2 = ro.conversion.py2rpy(metadata_df)
    # run lisi calculation
    lisi_score = lisi.compute_lisi(r_df1,r_df2,StrVector(metadata_keys))
    return lisi_score

def importr_tryhard(packname, devtools_path=None,contriburl=None):
    from rpy2.robjects.packages import importr
    from rpy2.robjects.packages import PackageNotInstalledError
    try:
        rpack = importr(packname)
    except PackageNotInstalledError:
        print("WANRING: please install " + packname)
        return
    return rpack


def downsample_rna(mat, rna_ratio, bycol=False): 
    """A function to compute lisi, batch-mixing metric"""
    # uses r-lisi package. Install all necessary packages if already. Assumes R is installed.
    
    # ===== load rpy2 =====#
    try:
        __import__("rpy2")
    except ImportError:
        import os
        os.system('pip install rpy2')      
        
    # function: ensure R package is installed
        # default - download from CRAN 
        # if devtools_path present, then use devtools to download from a github
        # if contriburl present, download from that url
    # ===== load necessary functions=====##
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    from rpy2.robjects.vectors import StrVector
    from rpy2.robjects.packages import importr
    #utils = importr('utils')

    # loads DropletUtils library, print warning message if not installed
    scuttle = importr_tryhard("scuttle")
    base = importr('base')
    
    # convert df to r-object
    r_mat = ro.conversion.py2rpy(mat)
    r_rna_ratio = ro.conversion.py2rpy(rna_ratio)
    r_bycol = ro.conversion.py2rpy(bycol)
    
    # run downsampling
    r_mat_down = base.as_matrix(
        scuttle.downsampleMatrix(
            r_mat, r_rna_ratio, 
            bycol = r_bycol)
    )
    return r_mat_down

def downsample_atac(bc, atac_ratio, frag_out, 
                   frag_path =  "",
                    pkg_source="/home/myylee/anaconda3/envs/scib2/bin/",
                feature_path="",
                    ncore=1, cal_mat=True): 
    """A function to compute lisi, batch-mixing metric"""
    # uses r-lisi package. Install all necessary packages if already. Assumes R is installed.
    import os
    # ===== load rpy2 =====#
    try:
        __import__("rpy2")
    except ImportError:
        os.system('pip install rpy2')      

    # ===== load necessary functions=====##
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    from rpy2.robjects.vectors import StrVector
    from rpy2.robjects.packages import importr
    #utils = importr('utils')

    # load r_utils.R under the same directory 
    import rpy2.robjects as robjects
    r = robjects.r
    print("current directory: "+ os. getcwd()+"; loading r_utils.R")
    r.source('r_utils.R')
    base = importr('base')
    
    # convert df to r-object
    r_bc = ro.conversion.py2rpy(bc)
    r_atac_ratio = ro.conversion.py2rpy(atac_ratio)
    os.makedirs(frag_out, exist_ok=True)
    r_frag_out = ro.conversion.py2rpy(frag_out)   
    
    # ==== default parameters ===== 
    r_pkg_source = ro.conversion.py2rpy(pkg_source)
    r_feature_path = ro.conversion.py2rpy(feature_path)
    r_ncore = ro.conversion.py2rpy(ncore)
    r_cal_mat = ro.conversion.py2rpy(cal_mat)
    r_frag_path = ro.conversion.py2rpy(frag_path)

    atac_mat_down = base.as_matrix(r.downsampleFragments(r_bc,r_atac_ratio,r_frag_out,frag_path=r_frag_path,
                                                         pkg_source=r_pkg_source,feature_path=r_feature_path,
                                                         ncore=r_ncore,cal_mat=r_cal_mat))
    print("done preparing atac_mat_down")
    return atac_mat_down


def purity_score(y_true, y_pred):
    """A function to compute cluster purity"""
    # compute contingency matrix (also called confusion matrix)
    contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)

    return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix)


def clustering_metrics(pred,truth):
    metrics_ = [ari, nmi, purity_score]
    ARI, NMI, Purity = [metric(pred, truth) for metric in metrics_]
    return((ARI, NMI, Purity))

    
def export_cobolt_res(model,resolution=0.5,out_dir=None):
    from cobolt.utils import SingleData, MultiomicDataset
    from cobolt.model import Cobolt

    latent = model.get_all_latent()
    res_df = pd.DataFrame(latent[0])
    res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")
    res_df["barcodes"] = latent[1]
    clusters = model.get_clusters(algo="leiden", resolution=resolution)
    res_df["leiden"] = clusters
    res_df["UMAP_1"] = model.reduction["UMAP2"]["embedding"][:,0]
    res_df["UMAP_2"] = model.reduction["UMAP2"]["embedding"][:,1]
    res_df["dataset"] = np.array([model.dataset.dataset[b] for b in res_df["barcodes"]])
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        res_df.to_csv(out_dir+"cobolt_result.csv")
    return(res_df)

def find_resolution_cobolt(model_, n_clusters, algo="leiden"): 
    model = deepcopy(model_)
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1000.]
    
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        model.clustering(algo=algo, resolution=current_res)
        labels = model.get_clusters(algo=algo, resolution=current_res)
        obtained_clusters = len(np.unique(labels))
        
        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res
        
        iteration = iteration + 1
        
    return current_res


# add s to every elements in v, add s before v if before = True, after v otherwise
def append_string(v,s, before=True):
    if before:
        return [s + vi  for vi in v]
    else:
        return [vi + s  for vi in v]



def tf_idf(adata):
    """A function to perfrom TF-IDF calculation on adata.X and save as a new layer "normalized".
    TD-IDF calculation is following the Signac paper (https://www.biorxiv.org/content/10.1101/2020.11.09.373613v1)
    
    
    Arguments:
    ------------------------------------------------------------------
    - adata: `AnnData`
    
    Returns:
    ------------------------------------------------------------------
    - adata: `AnnData`, with extra layer called "normalized" stored
    """
    X = adata.X
    # convert X to sparse matrix 
    X_csr = csr_matrix(X)
    rsums = X_csr.sum(axis=0)
    rsums = np.asarray(rsums).reshape(-1)
    if len(rsums==0) ==0:
        print("some peaks have 0 counts, please filter before running this function")
        #return adata

    npeaks = X_csr.sum(axis=1)
    npeaks = np.asarray(npeaks).reshape(-1)
    tf = X_csr.transpose().dot(csr_matrix(np.diag(1/npeaks)))
    print(tf.shape)
    idf = X_csr.shape[0]/rsums
    idf = np.asarray(idf).reshape(-1)
    idf_csr = csr_matrix(np.diag(idf))
    print(idf_csr.shape)
    norm_data_csr = idf_csr.dot(tf)
    norm_data_csr_log = norm_data_csr*10000
    norm_data_csr_log = norm_data_csr_log.log1p()
    
    adata.layers["normalized"] = norm_data_csr_log.transpose()
    return adata


def svd(adata,k,scale_embed,scale_max=10):
    """A function to perfrom SVD calculation on adata.layers["normalized"](if exist). This is the second part of LSI dimensional reduction.
    TD-IDF calculation is following the Signac paper (https://www.biorxiv.org/content/10.1101/2020.11.09.373613v1) and Signac's implmentation
    of RunSVD(https://github.com/timoast/signac/blob/master/R/dimension_reduction.R)
    
    
    Arguments:
    ------------------------------------------------------------------
    - adata: `AnnData`
    - k: `int`, the number of dimension to reduce to
    - scale_embed: `boolean`, indicate whether to standardize the cell embeddings
    - scale_max: `int or None`, represent the max value of cell embeddings
    
    Returns:
    ------------------------------------------------------------------
    - adata: `AnnData`, with extra obsm, varm and uns added
    """
    if "normalized" not in list(adata.layers):
        X = adata.X
    else: 
        X = adata.layers["normalized"]
    u, s, vt = svds(X, k=k)
    idx = np.argsort(-s)
    s = s[idx]
    u = u[:,idx]
    vt = vt[idx,:]
    
    sdev = s / np.sqrt(max(1, X.shape[1] - 1))
    cell_embeddings = u
    if scale_embed:
        print("Scaling cell embeddings")
        embed_mean = np.mean(cell_embeddings,axis=0)
        embed_sd = np.std(cell_embeddings,axis=0)
        norm_embeddings = (cell_embeddings - embed_mean)/ (embed_sd)
        if scale_max is not None:
            norm_embeddings[norm_embeddings > scale_max] <- scale_max
            norm_embeddings[norm_embeddings < -scale_max] <- -scale_max
    else:
        norm_embeddings = cell_embeddings
        
    
    adata.obsm['cell_embeddings'] = norm_embeddings
    adata.varm['feature_loadings'] = vt.transpose()
    adata.uns['lsi_sdev'] = sdev
    
    return adata

 

# mean l2 norm between target and pred 
def modal_msd(pred,target,pred_id,target_id,label=None):
    """A function to calculates average L2 distance between predicted and target, which should be multidimensional latent space coodinates.
    
    Arguments:
    ------------------------------------------------------------------
    - pred: `Numpy ndarray`, predicted coordinates in the latent space (predicted from another modality)
    - target: `Numpy ndarray`, target coordinates in the latent space 
    - pred_id: `list/Panda dataframe`, cell id used to find its corresponding cells in target data. All numbers must exist in target_id.
    - target_id: `list/Panda dataframe`, cell id. Could be longer than pred_id
    - label: 'Panda dataframe or None'. If not None, same length as pred_id and a label specific msd is calculated for each label value 
    
    Returns:
    ------------------------------------------------------------------
    - msd: `panda dataframe`, average l2 distance. If label is present, row is indexed by unique label values, else, the index is 'All'
    """
    
    pred_id = list(pred_id)
    target_id = list(target_id)
    # reorder target_id to match pred_id, get relative idx 
    idx = [target_id.index(i) for i in pred_id]
    # reorder target rows 
    target = target[idx,:]
    if label is not None:
        label = np.array(label)
        label_uni = np.unique(label)
        msd = np.zeros(len(label_uni))
        for i, label_i in enumerate(label_uni):
            idx_i = (label==label_i)
            pred_i = pred[idx_i,:]
            target_i = target[idx_i,:]
            msd[i] = np.mean(np.linalg.norm(np.array(target_i)- np.array(pred_i),axis=1)**2)
        msd_pd = pd.DataFrame(data=msd, index=label_uni)
    else:
        msd = np.mean(np.linalg.norm(np.array(target)- np.array(pred),axis=1)**2)
        msd_pd = pd.DataFrame(data=[msd],index=['All'])
        
        
    return msd_pd