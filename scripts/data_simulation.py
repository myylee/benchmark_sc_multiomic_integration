# If symlink is broken or it exists, return a boolean indicating whether symlink can be deleted
def delete_symlink(filename):
    import os 
    is_symlink = os.path.islink(filename) 
    file_exists = os.path.exists(filename)
    bad_symlink = is_symlink and not file_exists
    return (bad_symlink or file_exists)
    
# dataset split function that requires number of cells to be specified for each of the three data types
# specify where the fragment file is and will create symlink 
# patent_dir is used to check if there is already fragment files 
def pair_unpair_split_size(adata_rna,adata_atac,n_multi,n_urna=None,n_uatac=None,
                           transpose=False,fragment_path=None,parent_dir=None):
    import numpy as np
    from numpy.random import choice
    from numpy import setdiff1d
    from math import floor
    import os 
    from copy import deepcopy
    from anndata import AnnData
    from scipy.sparse import csr_matrix
    # randomly sample from 'vector' with the length defined by size
    subsample_size = lambda vector, size: choice(vector, size = size, replace = False) 
    
    if transpose:
        adata_rna = adata_rna.copy().transpose()
        adata_atac = adata_atac.copy().transpose()
    assert(adata_rna.shape[0] == adata_atac.shape[0])
    adata_rna = deepcopy(adata_rna)
    adata_atac = deepcopy(adata_atac)
    ncells = adata_rna.shape[0]
    paired_idx  = subsample_size(range(ncells),n_multi)
    unpaired_idx = setdiff1d(range(ncells), paired_idx)
    if n_urna is None and n_uatac is None:
        n_urna = round(len(unpaired_idx)/2)
        n_uatac = len(unpaired_idx) - n_urna
    elif n_urna is not None and n_uatac is None:
        n_uatac = len(unpaired_idx) - n_urna
    elif n_urna is None and n_uatac is not None:
        n_urna = len(unpaired_idx) - n_uatac
    elif (n_urna+n_uatac+n_multi) > adata_rna.shape[0]:
        print("sum of all cells greater than available, rerun the function")
        return
    u1_idx = subsample_size(unpaired_idx,n_urna)
    remain_idx = setdiff1d(unpaired_idx, u1_idx)
    u2_idx = subsample_size(remain_idx,n_uatac)

    adata_p1 = AnnData(csr_matrix(deepcopy(adata_rna.X[paired_idx,:].todense())),
                       obs=adata_rna.obs.iloc[paired_idx,:],
                       var=adata_rna.var,dtype=np.float32)
    adata_p2 = AnnData(csr_matrix(deepcopy(adata_atac.X[paired_idx,:].todense())),
                       obs=adata_atac.obs.iloc[paired_idx,:],
                       var=adata_atac.var,dtype=np.float32)
    adata_u1 = AnnData(csr_matrix(deepcopy(adata_rna.X[u1_idx,:].todense())),
                       obs=adata_rna.obs.iloc[u1_idx,:],
                       var=adata_rna.var,dtype=np.float32)
    adata_u2 = AnnData(csr_matrix(deepcopy(adata_atac.X[u2_idx,:].todense())),
                       obs=adata_atac.obs.iloc[u2_idx,:],
                       var=adata_atac.var,dtype=np.float32)
                    
    if transpose:
        adata_p1 = adata_p1.copy().transpose()
        adata_p2 = adata_p2.copy().transpose()
        adata_u1 = adata_u1.copy().transpose()
        adata_u2 = adata_u2.copy().transpose()
    
    if (fragment_path is not None) and (parent_dir is not None):
        # copy for paired_atac
        # copy for unpaired_atac
        from os.path import exists
        tbi_file = fragment_path+".tbi"
        if exists(fragment_path) and exists(tbi_file):
            pfolder = parent_dir+"/paired_ATAC/"
            ufolder = parent_dir+"/unpaired_ATAC/"
            os.makedirs(pfolder, exist_ok=True)
            os.makedirs(ufolder, exist_ok=True)
            
            if delete_symlink(pfolder+"fragments.tsv.gz"):
                os.remove(pfolder+"fragments.tsv.gz")
            if delete_symlink(pfolder+"fragments.tsv.gz.tbi"):
                os.remove(pfolder+"fragments.tsv.gz.tbi")
            if delete_symlink(ufolder+"fragments.tsv.gz"):
                os.remove(ufolder+"fragments.tsv.gz")
            if delete_symlink(ufolder+"fragments.tsv.gz.tbi"):
                os.remove(ufolder+"fragments.tsv.gz.tbi")
                
            os.symlink(fragment_path, pfolder+"fragments.tsv.gz")
            os.symlink(tbi_file, pfolder+"fragments.tsv.gz.tbi")
            #os.getcwd()+"/"+
            os.symlink(fragment_path, ufolder+"fragments.tsv.gz")
            os.symlink(tbi_file, ufolder+"fragments.tsv.gz.tbi")
    
    return((adata_p1,adata_p2,adata_u1,adata_u2))


# simulate with batch structure 
#=== if same_unpair_origin == True, assert that adata_urna.obs_names == adata_uatac.obs_names. And that n_urna + n_uatac <= adata_urna.shape[0]
#=== if same_unpair_origin == False, then adata_urna and adata_uatac are from different adata 
def pair_unpair_split_size_batch(adata_prna,adata_patac,adata_urna,adata_uatac,n_multi,n_urna=None,n_uatac=None,
                                 transpose=False,fragment_path=None,parent_dir=None,same_unpair_origin=False):
    import numpy as np
    from numpy.random import choice
    from numpy import setdiff1d
    from math import floor
    import os 
    from copy import deepcopy
    from anndata import AnnData
    from scipy.sparse import csr_matrix
    # randomly sample from 'vector' with the length defined by size
    subsample_size = lambda vector, size: choice(vector, size = size, replace = False) 
    
    if transpose:
        adata_prna = adata_prna.copy().transpose()
        adata_patac = adata_patac.copy().transpose()
        adata_urna = adata_urna.copy().transpose()
        adata_uatac = adata_uatac.copy().transpose()
        
    assert adata_prna.shape[0] == adata_patac.shape[0]
    
    ncells = adata_prna.shape[0]
    ncells_urna = adata_urna.shape[0]
    ncells_uatac = adata_uatac.shape[0]
    

    adata_prna = deepcopy(adata_prna)
    adata_patac = deepcopy(adata_patac)
    adata_urna = deepcopy(adata_urna)
    adata_uatac = deepcopy(adata_uatac)
    
    assert ncells>=n_multi,"no enough number of cells to generate paired sample"
    paired_idx  = subsample_size(range(ncells),n_multi)
    
    # if unpaired_RNA and unpaired_ATAC are from the same sample, cells with the same barcode will either be in RNA or ATAC, but not both
    if same_unpair_origin: 
        assert adata_urna.obs_names.tolist() == adata_uatac.obs_names.tolist(), "the input data for unpaired RNA and ATAC samples have different cell barcodes or wrong order"
        assert ncells_urna >= (n_urna + n_uatac), "no enough number of cells to generate unpaired sample"
        u1_idx = subsample_size(range(ncells_urna),n_urna)
        remain_idx = setdiff1d(range(ncells_urna), u1_idx)
        u2_idx = subsample_size(remain_idx,n_uatac)

    else:
    # if unpaired_RNA and unpaired_ATAC are not from the same sample, different reference adata is used to simulate unpaired_RNA and unpaired_ATAC
        assert ncells_urna >= n_urna, "no enough number of cells to generate unpaired RNA"
        assert ncells_uatac >= n_uatac,"no enough number of cells to generate unpaired ATAC"
        u1_idx  = subsample_size(range(ncells_urna),n_urna)
        u2_idx  = subsample_size(range(ncells_uatac),n_uatac)
    
    adata_p1 = AnnData(csr_matrix(deepcopy(adata_prna.X[paired_idx,:].todense())),
                       obs=adata_prna.obs.iloc[paired_idx,:],
                       var=adata_prna.var,dtype=np.float32)
    adata_p2 = AnnData(csr_matrix(deepcopy(adata_patac.X[paired_idx,:].todense())),
                       obs=adata_patac.obs.iloc[paired_idx,:],
                       var=adata_patac.var,dtype=np.float32)
    adata_u1 = AnnData(csr_matrix(deepcopy(adata_urna.X[u1_idx,:].todense())),
                       obs=adata_urna.obs.iloc[u1_idx,:],
                       var=adata_urna.var,dtype=np.float32)
    adata_u2 = AnnData(csr_matrix(deepcopy(adata_uatac.X[u2_idx,:].todense())),
                       obs=adata_uatac.obs.iloc[u2_idx,:],
                       var=adata_uatac.var,dtype=np.float32)
                    
    if transpose:
        adata_p1 = adata_p1.copy().transpose()
        adata_p2 = adata_p2.copy().transpose()
        adata_u1 = adata_u1.copy().transpose()
        adata_u2 = adata_u2.copy().transpose()
    
    if (fragment_path is not None) and (parent_dir is not None):
        # copy for paired_atac
        # copy for unpaired_atac
        from os.path import exists
        tbi_file = fragment_path+".tbi"
        if exists(fragment_path) and exists(tbi_file):
            pfolder = parent_dir+"/paired_ATAC/"
            ufolder = parent_dir+"/unpaired_ATAC/"
            os.makedirs(pfolder, exist_ok=True)
            os.makedirs(ufolder, exist_ok=True)
            
            if delete_symlink(pfolder+"fragments.tsv.gz"):
                os.remove(pfolder+"fragments.tsv.gz")
            if delete_symlink(pfolder+"fragments.tsv.gz.tbi"):
                os.remove(pfolder+"fragments.tsv.gz.tbi")
            if delete_symlink(ufolder+"fragments.tsv.gz"):
                os.remove(ufolder+"fragments.tsv.gz")
            if delete_symlink(ufolder+"fragments.tsv.gz.tbi"):
                os.remove(ufolder+"fragments.tsv.gz.tbi")
                
            os.symlink(fragment_path, pfolder+"fragments.tsv.gz")
            os.symlink(tbi_file, pfolder+"fragments.tsv.gz.tbi")
            #os.getcwd()+"/"+
            os.symlink(fragment_path, ufolder+"fragments.tsv.gz")
            os.symlink(tbi_file, ufolder+"fragments.tsv.gz.tbi")
    
    return((adata_p1,adata_p2,adata_u1,adata_u2))

# requires r-scuttle
# downsamples fragment file by default. 
# if depth_x < 1, fragment downsampling will happen for ATAC datasets. Fragment symlink file created previously will be removed and a new downsampled fragment file will be created.
# if depth_x < 1, cell-peak counts matrix will be counted again.
def downsample_samples(adata_p1,adata_p2,adata_u1,adata_u2,depth_multi, depth_urna, depth_uatac,downsample_fragment=True,parent_dir=None,frag_path=""):
    import utils_eval
    import pandas as pd
    import os 
    
    assert depth_multi > 0 and depth_multi <= 1, "please input depth_multi in the correct range (0,1]"
    assert depth_urna > 0 and depth_urna <= 1, "please input depth_urna in the correct range (0,1]"
    assert depth_uatac > 0 and depth_uatac <= 1, "please input depth_uatac in the correct range (0,1]"

    # ==== RNA portion (run function if ratio != 1) ====
    if depth_multi != 1:
        print("downsampling paired RNA")
        adata_p1.X = utils_eval.downsample_rna(adata_p1.X.todense().T,depth_multi,bycol=True).T
    if depth_urna != 1:
        print("downsampling unpaired RNA")
        adata_u1.X = utils_eval.downsample_rna(adata_u1.X.todense().T,depth_urna).T

    # ==== ATAC portion ====
    if downsample_fragment:
        if parent_dir is None:
            print("please specify where the fragment files should be saved")
            return
        else: 
            pd.DataFrame(adata_p2.var_names).to_csv("tmp.csv",header=False,index=False)
            if depth_multi < 1:
                print("downsampling paired ATAC")
                print("deleting previously written symlink file(fragments.tsv.gz)")
                pfolder = parent_dir+"/paired_ATAC/"
                if delete_symlink( os.path.join(pfolder,"fragments.tsv.gz")):
                    os.remove( os.path.join(pfolder,"fragments.tsv.gz"))
                if delete_symlink( os.path.join(pfolder,"fragments.tsv.gz.tbi")):
                    os.remove( os.path.join(pfolder,"fragments.tsv.gz.tbi"))
                    
                adata_p2.X = utils_eval.downsample_atac(list(adata_p2.obs.index),
                                                        depth_multi,
                                                        os.path.join(parent_dir,"paired_ATAC"),
                                                        frag_path=frag_path,
                                                        feature_path="tmp.csv").T
            if depth_uatac < 1:
                print("downsampling unpaired ATAC")
                print("deleting previously written symlink file(fragments.tsv.gz)")
                pfolder = os.path.join(parent_dir,"unpaired_ATAC")
                if delete_symlink( os.path.join(pfolder,"fragments.tsv.gz")):
                    os.remove( os.path.join(pfolder,"fragments.tsv.gz"))
                if delete_symlink( os.path.join(pfolder,"fragments.tsv.gz.tbi")):
                    os.remove( os.path.join(pfolder,"fragments.tsv.gz.tbi"))
                adata_u2.X = utils_eval.downsample_atac(list(adata_u2.obs.index),
                                                        depth_uatac,
                                                         os.path.join(parent_dir,"unpaired_ATAC"),
                                                        frag_path=frag_path,
                                                        feature_path="tmp.csv").T
            os.remove("tmp.csv")
    else:
        if depth_multi != 1:
            print("downsampling paired ATAC")
            adata_p2.X = utils_eval.downsample_rna(adata_p2.X.todense().T,depth_multi,bycol=True).T
        if depth_uatac != 1:
            print("downsampling unpaired ATAC")
            adata_u2.X = utils_eval.downsample_rna(adata_u2.X.todense().T,depth_uatac).T

    return(adata_p1,adata_p2,adata_u1,adata_u2)

    
# batch and wait_time only matters if wait is true
def eval_test(folder_dir,script,eval_script,conda_envs,scripts,py_langs,file_paths,
              ct_ref,nclust,dir_path,cond_key,iter_list,repeats,repeat_start=1,wait=True,
              wait_time=10*60,batch=1,output_folder="results",ncore=8):
    import os
    from os.path import exists
    import time
    import subprocess
    
    counter=0
    file_check_list= []
    for j in range(repeat_start,repeats+1):
        for i in range(len(iter_list)):
            in_dir_i = "{}{}{}_{}/".format(dir_path,cond_key,iter_list[i],j)
            out_dir_i = "{}{}{}_{}/{}/".format(dir_path,cond_key,iter_list[i],j,output_folder)

            for m in range(len(conda_envs)):
                print("processing {} repeat {} using {}".format(in_dir_i,j,scripts[m]))
                lang_string = "-p" if py_langs[m] else "-r"

                job_arg = ["bsub",
                           "-q","mingyao_normal",
                           "-n", str(ncore),
                          "-o", "{}/job_outs/{}{}_{}_results_{}_%J.txt".format(folder_dir,cond_key,iter_list[i],j,conda_envs[m]),
                          "-R","span[hosts=1]"
                          ]

                script_arg =["\'sh" ,script,
                             "-i", in_dir_i,
                             "-w", out_dir_i,
                             "-c", conda_envs[m],
                             "-s", scripts[m],
                             lang_string,
                             "-e", eval_script,
                             "-f", file_paths[m],
                             "-t", ct_ref,
                             "-l", str(nclust),
                             "\'"
                            ]
                bsub_line  = ' '.join(job_arg) + " " + ' '.join(script_arg)
                subprocess.check_call(bsub_line,shell=True)
                if wait:
                    # wait until this job is done to submit the next
                    counter +=1 
                    file_check_list.append(out_dir_i+file_paths[m])
                    if (counter == batch):
                        # if not all exists, then keep checking 
                        while(not all(list(map(exists, file_check_list)))):
                            time.sleep(wait_time)
                        counter = 0
                        file_check_list = []
                        

#====== Compatible with evaluate_vary_situations_all.ipynb and submit_job_per_condition_n_eval2.sh ======#
#=== folder_dir: working directory 
#=== script: job submission script 
#=== conda_envs: a list of conda environments 
#=== scripts: a list of method_scripts
#=== py_langs: a list of indicators for whether the script is based on python or R 
#=== file_paths: a list of output file location for evaluation script to be run on 
#=== method_keys: a list of strings indicating what method the code is based on 
#=== ct_ref: path to the cell type reference file 
#=== nclust: integer indicating number of clusters
#=== dir_path: input directory path 
#=== cond_key: condition key that structures input and output folder 
#=== iter_list: a list of variables that is varying 
#=== repeats: integer, repeats to run each method
#=== repeat_start: integer, start index for replicate
#*** replications to be run will be [repeat_start, repeats] 
#=== wait: boolean, whether to wait for the script to finish
#=== wait_time: seconds to wait before checking for each iteration 
#=== batch: integer, number of runs to submit for each iteration 
#*** batch and wait_time only matters if wait is true
#=== output_folder: string, output folder for results to be stored 
#=== ncore: integer, number of cores to use 
#=== gp_script: string|None, path to gene-peak pair evaluation script 
#=== gp_truth: string|None, path to gene-peak pair ground truth csv file 

def eval_test_all(folder_dir,script,eval_script,conda_envs,scripts,py_langs,file_paths,
                  method_keys,ct_ref,nclust,dir_path,cond_key,iter_list,repeats,repeat_start=1,wait=True,
                  wait_time=10*60,batch=1,output_folder="results",ncore=8,gp_script=None,gp_truth=None):
    import os
    from os.path import exists
    import time
    import subprocess
    
    counter=0
    file_check_list= []
    for j in range(repeat_start,repeats+1):
        for i in range(len(iter_list)):
            in_dir_i = os.path.join(dir_path,"{}{}_{}/".format(cond_key,iter_list[i],j))
            out_dir_i = os.path.join(dir_path,"{}{}_{}/".format(cond_key,iter_list[i],j),output_folder)

            for m in range(len(conda_envs)):
                print("processing {} repeat {} using {}".format(in_dir_i,j,scripts[m]))
                lang_string = "-p" if py_langs[m] else "-r"

                job_arg = ["bsub",
                           "-q","mingyao_normal",
                           "-n", str(ncore),
                          "-o", "{}/job_outs/{}{}_{}_results_{}_%J.txt".format(folder_dir,cond_key,iter_list[i],j,conda_envs[m]),
                          "-R","span[hosts=1]",
                           "-M","32GB"
                          ]
                if gp_script == None:
                    gp_script = "false"
                    
                if gp_truth == None:
                    gp_truth = "false"

                script_arg =["\'sh" ,script,
                             "-i", in_dir_i,
                             "-w", out_dir_i,
                             "-c", conda_envs[m],
                             "-s", scripts[m],
                             lang_string,
                             "-e", eval_script,
                             "-f", file_paths[m],
                             "-m", method_keys[m],
                             "-t", ct_ref,
                             "-l", str(nclust),
                             "-a", gp_script,
                             "-b", gp_truth,
                             "\'"
                            ]
                bsub_line  = ' '.join(job_arg) + " " + ' '.join(script_arg)
                subprocess.check_call(bsub_line,shell=True)
                if wait:
                    # wait until this job is done to submit the next
                    counter +=1 
                    file_check_list.append(os.path.join(out_dir_i,file_paths[m]))
                    if (counter == batch):
                        # if not all exists, then keep checking 
                        while(not all(list(map(exists, file_check_list)))):
                            time.sleep(wait_time)
                        counter = 0
                        file_check_list = []


def data_simulation(in_dir,adata_rna,adata_atac,iter_list,depth_multiome_list,depth_scrna_list,depth_snatac_list,
                   n_multiome_list,n_scrna_list,n_snatac_list,repeats,fragment_path,cond_key,
                   fn, downsample=False):
    import utils_eval
    for i in range(len(iter_list)): 
        s_i = iter_list[i]
        for j in range(1,repeats+1):
            folder_path = "{}{}{}_{}/".format(in_dir,cond_key,fn(s_i),str(j))
            print("spliting: paired dataset has {} cells, repeat {}".format(s_i,j))
            #print("spliting: paired dataset is {} depth, repeat {}".format(s_i,j))
            print("folder name :"+folder_path)
            #if not downsample:
            (adata_p1,adata_p2,adata_u1,adata_u2) = pair_unpair_split_size(adata_rna,
                                                                           adata_atac,
                                                                           n_multi = n_multiome_list[i],
                                                                           n_urna = n_scrna_list[i],
                                                                           n_uatac = n_snatac_list[i],
                                                                           transpose=False,
                                                                           fragment_path=fragment_path,
                                                                           parent_dir=folder_path)
#             else:
#                 (adata_p1,adata_p2,adata_u1,adata_u2) = pair_unpair_split_size(adata_rna,
#                                                                                adata_atac,
#                                                                                n_multi = n_multiome_list[i],
#                                                                                n_urna = n_scrna_list[i],
#                                                                                n_uatac = n_snatac_list[i],
#                                                                                transpose=False,
#                                                                                fragment_path=None,
#                                                                                parent_dir=None)
            if downsample:
                # downsampling ATAC from fragment files if depth ratio < 1
                (adata_p1,adata_p2,adata_u1,adata_u2) = downsample_samples(adata_p1,adata_p2,adata_u1,adata_u2,
                                                                           depth_multi = depth_multiome_list[i],
                                                                           depth_urna = depth_scrna_list[i],
                                                                           depth_uatac = depth_snatac_list[i],
                                                                           downsample_fragment=True,
                                                                           parent_dir=folder_path,
                                                                           frag_path = fragment_path)


            utils_eval.write_adata(adata_p1, folder_path+"/paired_RNA/","RNA","gene",feature_name='feature',transpose=True)
            utils_eval.write_adata(adata_p2, folder_path+"/paired_ATAC/","ATAC","peak",feature_name='feature',transpose=True)
            utils_eval.write_adata(adata_u1, folder_path+"/unpaired_RNA/","RNA","gene",feature_name='feature',transpose=True)
            utils_eval.write_adata(adata_u2, folder_path+"/unpaired_ATAC/","ATAC","peak",feature_name='feature',transpose=True)

            
def data_simulation_batch(in_dir,adata_prna,adata_patac,adata_urna,adata_uatac,
                          iter_list,depth_multiome_list,depth_scrna_list,depth_snatac_list,
                          n_multiome_list,n_scrna_list,n_snatac_list,repeats,fragment_path,
                          cond_key,fn, downsample=False, same_unpair_origin=False):
    import utils_eval
    for i in range(len(iter_list)): 
        s_i = iter_list[i]
        for j in range(1,repeats+1):
            folder_path = "{}{}{}_{}/".format(in_dir,cond_key,fn(s_i),str(j))
            print("spliting: paired dataset has {} cells, repeat {}".format(s_i,j))
            #print("spliting: paired dataset is {} depth, repeat {}".format(s_i,j))
            print("folder name :"+folder_path)
            #if not downsample:
            (adata_p1,adata_p2,adata_u1,adata_u2) = pair_unpair_split_size_batch(adata_prna,
                                                                               adata_patac,
                                                                               adata_urna,
                                                                               adata_uatac,
                                                                               n_multi = n_multiome_list[i],
                                                                               n_urna = n_scrna_list[i],
                                                                               n_uatac = n_snatac_list[i],
                                                                               transpose=False,
                                                                               fragment_path=fragment_path,
                                                                               parent_dir=folder_path,
                                                                                 same_unpair_origin=same_unpair_origin)
            if downsample:
                # downsampling ATAC from fragment files if depth ratio < 1
                (adata_p1,adata_p2,adata_u1,adata_u2) = downsample_samples(adata_p1,adata_p2,adata_u1,adata_u2,
                                                                           depth_multi = depth_multiome_list[i],
                                                                           depth_urna = depth_scrna_list[i],
                                                                           depth_uatac = depth_snatac_list[i],
                                                                           downsample_fragment=True,
                                                                           parent_dir=folder_path,
                                                                           frag_path = fragment_path)


            utils_eval.write_adata(adata_p1, folder_path+"/paired_RNA/","RNA","gene",
                                   feature_name='feature',transpose=True,additional_obs=['batch'])
            utils_eval.write_adata(adata_p2, folder_path+"/paired_ATAC/","ATAC","peak",
                                   feature_name='feature',transpose=True,additional_obs=['batch'])
            utils_eval.write_adata(adata_u1, folder_path+"/unpaired_RNA/","RNA","gene",
                                   feature_name='feature',transpose=True,additional_obs=['batch'])
            utils_eval.write_adata(adata_u2, folder_path+"/unpaired_ATAC/","ATAC","peak",
                                   feature_name='feature',transpose=True,additional_obs=['batch'])
