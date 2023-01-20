# Benchmarking algorithms for joint integration of unpaired and paired single-cell RNA-seq and ATAC-seq data
**Michelle Y. Y. Lee, Klaus H. Kaestner, Mingyao Li**

This repository contains codes for data simulation and method evaluation described in this study. We evaluated seven methods at three types of sceanrios using two publically available datasets. 

## Benchmarking framework 

Codes for data simulation and method evaluation are in evaluate_vary_situations_public.ipynb. 

Codes for generating summary plots are in summary_metric_plot_public.ipynb. 

Each file is separated into three major sceanrios and a total of eleven challenges. For each challenge, first two steps are in evaluate_vary_situations_public.ipynb. First, data is simulated. Secondly, a list of methods are run and the performance is calculated. To generate a summary plot for each challenge, follow the steps in  summary_metric_plot_public.ipynb, which is again separated into the same three sceanrios and eleven challenges. 

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
Data simulation, evaluation, and plot generation are run after activating the scib2 environment. Details about environment setup can be ofund in the **methods/eval_scib2** folder. 

## Data 
Proceesed source data and reference files used during data simulation and evaluations, as well as scripts to generate these files from raw dataset are in the **dataset** folder [coming soon]. 
 

