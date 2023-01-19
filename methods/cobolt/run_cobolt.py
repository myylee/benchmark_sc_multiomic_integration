#!/usr/bin/python
import sys 
print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:'+ str(sys.argv))

from cobolt.utils import SingleData, MultiomicDataset
from cobolt.model import Cobolt
import os
import pandas as pd
import numpy as np
import pickle
import timeit

def run_cobolt_fn(in_dir,out_dir):
    start = timeit.default_timer()
    # Read the SNARE-seq gene expression data.
    paired_rna = SingleData.from_file(path=os.path.join(in_dir, "paired_RNA"),
                                      dataset_name="Multiome",
                                      feature_name="GeneExpr",
                                      count_file="RNA_counts.mtx",
                                      feature_file="gene.tsv")

    # Read the SNARE-seq chromatin accessibility data.
    paired_atac = SingleData.from_file(path=os.path.join(in_dir, "paired_ATAC"),
                                      dataset_name="Multiome",
                                      feature_name="ChromAccess",
                                      count_file="ATAC_counts.mtx",
                                      feature_file="peak.tsv")

    unpaired_rna = SingleData.from_file(path=os.path.join(in_dir, "unpaired_RNA"),
                                    dataset_name="scRNA",
                                    feature_name="GeneExpr",
                                    feature_file="gene.tsv",
                                    count_file="RNA_counts.mtx")
    unpaired_atac = SingleData.from_file(path=os.path.join(in_dir, "unpaired_ATAC"),
                                dataset_name="snATAC",
                                feature_name="ChromAccess",
                                feature_file="peak.tsv",
                                count_file="ATAC_counts.mtx")


    # Quality filtering on features.
    paired_rna.filter_features(upper_quantile=0.99, lower_quantile=0.7)
    paired_atac.filter_features(upper_quantile=0.99, lower_quantile=0.7)
    unpaired_rna.filter_features(upper_quantile=0.99, lower_quantile=0.7)
    unpaired_atac.filter_features(upper_quantile=0.99, lower_quantile=0.7)
    
    multi_dt = MultiomicDataset.from_singledata(unpaired_rna, unpaired_atac, paired_rna, paired_atac)

    model = Cobolt(dataset=multi_dt, lr=0.001, n_latent=16)
    model.train(num_epochs=20)
    os.makedirs(os.path.join(out_dir,"runtime"), exist_ok=True)
    os.makedirs(os.path.join(out_dir,"cobolt"), exist_ok=True)
    model_out = os.path.join(out_dir,"cobolt","cobolt_model.pickle")
    pickle.dump(model, open(model_out, 'wb'))
    
    # save latent embedding as csv 
    latent = model.get_all_latent()
    res_df = pd.DataFrame(latent[0],index=latent[1])
    res_df = res_df.set_axis(["latent_" + s  for s in res_df.columns.astype("str").tolist()],axis="columns")
    res_df['dataset'] = np.array([model.dataset.dataset[b] for b in res_df.index])
    csv_out = os.path.join(out_dir, "cobolt","cobolt_result.csv")
    res_df.to_csv(csv_out)
    
    print("------ Done ------")
    stop = timeit.default_timer()
    print('Time(s): ', stop - start)
    runtime_out = os.path.join(out_dir, "runtime","cobolt_runtime.txt")
    print(stop - start,  file=open(runtime_out, 'w'))

    return(model)


print("argument 1:",sys.argv[1])
print("argument 2:",sys.argv[2])
mvi = run_cobolt_fn(sys.argv[1],sys.argv[2])


