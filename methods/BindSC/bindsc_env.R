# command line

conda create -n bindsc r-essentials r-base

conda activate bindsc
conda install -y igraph hdf5

R
# R
install.packages('Seurat')
install.packages('IRkernel')
IRkernel::installspec(name = 'bindsc', displayname = 'rbindsc')


system("conda install -y -c conda-forge r-devtools")
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(devtools)
install_github('KChen-lab/bindSC')


BiocManager::install(c("GenomeInfoDb","IRanges", "Rsamtools", "S4Vectors", "BiocGenerics"))
remotes::install_version("RSQLite", version = "2.2.5")
BiocManager::install(c("EnsDb.Hsapiens.v86","biovizBase"))
install.packages("Signac") 
remotes::install_github("mojaveazure/seurat-disk")


library(bindSC)
library(Seurat)
library(Signac)
library(stringr)
