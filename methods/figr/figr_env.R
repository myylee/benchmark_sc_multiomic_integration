# command line

conda create -n figr r-essentials r-base

conda activate figr
conda install -y igraph hdf5

R
# R
install.packages('Seurat')
install.packages('IRkernel')
IRkernel::installspec(name = 'figr', displayname = 'rfigr')


system("conda install -y -c conda-forge r-devtools")
install.packages("BiocManager")
library(devtools)

BiocManager::install(c("GenomeInfoDb","IRanges", "Rsamtools", "S4Vectors", "BiocGenerics"))
remotes::install_version("RSQLite", version = "2.2.5")
BiocManager::install(c("EnsDb.Hsapiens.v86","biovizBase"))
install.packages("Signac") 
remotes::install_github("mojaveazure/seurat-disk")

install.packages("optmatch")
BiocManager::install("chromVAR")
install.packages("pbmcapply")

system("conda install -c conda-forge r-ggrastr")
devtools::install_github("caleblareau/BuenColors")
BiocManager::install("ComplexHeatmap")
install.packages(c("networkD3","network","GGally","network"))

library(Seurat)
library(Signac)
library(stringr)

# figr download from github repo 
# git clone https://github.com/buenrostrolab/stimATAC_analyses_code.git

