[![DOI](https://zenodo.org/badge/580110879.svg)](https://zenodo.org/badge/latestdoi/580110879)

# Benchmarking algorithms for joint integration of unpaired and paired single-cell RNA-seq and ATAC-seq data
**Michelle Y. Y. Lee, Klaus H. Kaestner, Mingyao Li**

This repository contains codes for data simulation and method evaluation described in this study. We evaluated nine methods at five types of simulated scenarios using three publicly available datasets, plus one real data situation. 

## Benchmarking framework 

Codes for data simulation and method evaluation are in evaluate_vary_situations_public.ipynb. 

Codes for generating summary plots are in summary_metric_plot_public.ipynb. 

An example for generating UMAP plots are in umap_plot_generation_public.ipynb.

Plots related to the HPAP integration are generated using codes in hpap_result_plot_public.ipynb. 

Each file is separated into five major scenario plus one real data situation and a total of 17 challenges. For each challenge, first two steps are in evaluate_vary_situations_public.ipynb. First, data is simulated. Secondly, a list of methods are run and the performance is calculated. To generate a summary plot for each challenge, follow the steps in  summary_metric_plot_public.ipynb, which is again separated into the same six situations and 17 challenges. 

## Setup 
To run the methods as specified in the step above, one needs to create a conda environment for each method and install all the necessary packages. 

### Method-specific environment setup and execution
The list of methods evaluated are under the **methods** folder. Each subfolder contains files related to one method. 

Installation (one conda environment is created for each method) 
- Option 1: run the '.*_env.txt' file line-by-line in linux and R to install the method and its dependencies. 
- Option 2: install using the .yml file to create the conda environment. E.g. run the code below to install the conda environment to run Seurat v4 or Seurat v3. 
    ```
    conda env create -f seurat.yml
    ```
### Evaluation environment setup 
Data simulation, evaluation, and plot generation are run after activating the scib2 environment. Details about environment setup can be found in the **methods/eval_scib2** folder. 

## Data 
Processed source data and reference files used during data simulation and evaluations, as well as scripts to generate these files from raw dataset can be downloaded [here](https://upenn.box.com/s/jtua3rmmvzempjq55z4kj9xiij9dqoez).

