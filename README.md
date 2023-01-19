# Benchmarking algorithms for joint integration of unpaired and paired single-cell RNA-seq and ATAC-seq data
**Michelle Y. Y. Lee, Klaus H. Kaestner, Mingyao Li**

This repository contains codes for data simulation and method evaluation described in this study. We evaluated seven methods at three different sceanrios simulated using two publically available datasets. 

## Benchmarking framework 

For each sceanrio, there are three major steps:
1. Simulate the data [Python kernel]
2. Define integration task (run each method with the necessary arguments and setup) and submit as an non-interactive job to computer clusters (e.g. LPC/HPC) [Python kernel] 
3. Plot integration result [R kernel]

We have one conda environment for data simulation and performance evaluation. Codes associated with data simulation are in jupyter notebooks under python environment, with calls to R functions. Codes associated with running each method is also in jupyter notebook in which is submits 

Default: we used HPC to submit each integration task as one job, with 8 cores and 32 GB of RAM. 

## Method-specific environment setup and execution
The list of methods evaluated are under the 'methods' folder. Each subfolder contains files related to one method. 

Installation 
- Option 1: run the '.*_env.txt' file line-by-line in linux and R to install the method and its dependencies. One conda environment is created for each method.
- Option 2: install using the .yml file to create the conda environment. E.g. run the code below to install the conda environment to run Seurat v4 or Seurat v3. 
    ```
    conda env create -f seurat.yml
    ```


Execution 
- Function(s) to run the method. The function can be called by a bash script or from linux terminal. The function accepts requires 2 or 3 arguments (depending on if the number of clusters need to be specified).
    - Inputs:
        - in_dir
        - out_dir
        - ncluster
    - Outputs: 
        - 

## Evaluation
- Eval single
- Eval fair 


## Plotting 
- Summary plot 
- UMAP plot

## Dataset preparation 
- PBMC
- BMMC

---
# Extras! [To be removed]
## Sceanrios
3 sceanrios, breif description

### Sceanrio 1
- Simulation 
- Run method 
- Evaluation 

### Sceanrio 2
- Simulation 
- Run method 
- Evaluation 

### Sceanrio 3
- Simulation 
- Run method 
- Evaluation 



Notes:
- One giant jupyter notebook script for data simulation and method execution 
- one giant jupyter notebook with R kernel for plotting summary plot and (UMAP (maybe?))

- A folder for each source data, with jupyter notebook detailing how the data is preprocessed, pmat tabulated, and peak-gene pair called. 
