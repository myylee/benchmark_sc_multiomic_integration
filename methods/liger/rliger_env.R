# command line

conda create -n liger r-essentials r-base

conda activate liger
conda install -y igraph hdf5

R
# R
install.packages('Seurat')
install.packages('IRkernel')
IRkernel::installspec(name = 'liger', displayname = 'rliger')

install.packages('devtools')
system("conda install -c conda-forge r-devtools")
#system("conda install -c conda-forge r-gert")
library(devtools)
install_github('welch-lab/liger')

install.packages("BiocManager")
BiocManager::install(c("GenomeInfoDb","IRanges", "Rsamtools", "S4Vectors", "BiocGenerics"))
remotes::install_version("RSQLite", version = "2.2.5")
BiocManager::install(c("EnsDb.Hsapiens.v86","biovizBase"))
install.packages("Signac") 


library(rliger)
library(Seurat)
library(Signac)
library(stringr)
