#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(Seurat)
require(stringr)
require(Signac)
require(Matrix)
require(future)
source("r_utils.R")

orig_wd = getwd()
figr_path = "/home/myylee/scmint/methods_eval/stimATAC_analyses_code/R/"
setwd(figr_path)
source("optMatching_functions.R")
source("DORC_functions.R")
source("FigR_functions.R")

setwd(orig_wd)


run_rfigr_fn <- function(in_dir, out_dir){
    # starting time
    t1 <- Sys.time()
#     plan("multisession")
#     options(future.rng.onMisue = "ignore")
#     print(paste0("workers used:",nbrOfWorkers()))
#     options(future.globals.maxSize = 8000 * 1024^2)
    
    # load dataset 
    datasets = load_datasets(in_dir,obs=c("barcodes","batch"))
    paired_rna = datasets$paired_rna
    paired_atac = datasets$paired_atac
    unpaired_rna = datasets$unpaired_rna
    unpaired_atac = datasets$unpaired_atac
    
    # verify number of cells in each condition 
    dataset_vec <- rep(c("scRNA","snATAC","Multiome-RNA","Multiome-ATAC"),
                   c(ncol(unpaired_rna),
                     ncol(unpaired_atac),
                     ncol(paired_rna),
                    ncol(paired_atac)))
    names(dataset_vec) <- c(paste0(colnames(unpaired_rna)),
                            paste0(colnames(unpaired_atac)),
                            paste0("prna_",colnames(paired_rna)),
                            paste0("patac_",colnames(paired_atac)))
    print(table(dataset_vec))
    
    paired_rna <- RenameCells(paired_rna,add.cell.id = "prna",for.merge = FALSE)
    paired_atac <- RenameCells(paired_atac,add.cell.id = "patac",for.merge = FALSE)
    
    unpaired_atac@meta.data$technology <- "snATAC"
    unpaired_rna@meta.data$technology <- "scRNA"
    paired_atac@meta.data$technology <- "Multiome-ATAC"
    paired_rna@meta.data$technology <- "Multiome-RNA"

    # merging
    unpaired_rna <- merge(unpaired_rna,paired_rna)
    unpaired_atac <- merge(unpaired_atac,paired_atac)
    
    unpaired_rna@meta.data$group <- paste0(unpaired_rna$batch,"_",unpaired_rna$technology)
    unpaired_atac@meta.data$group <- paste0(unpaired_atac$batch,"_",unpaired_atac$technology)

    # Normalize gene expression and obtain highly variable genes 
    DefaultAssay(unpaired_rna) <- "RNA"
    unpaired_rna <- NormalizeData(unpaired_rna)
    rna.list <- SplitObject(unpaired_rna, split.by = "group")
        #select high variable features across samples 
        features = SelectIntegrationFeatures(
          rna.list,
          nfeatures = 5000,
          verbose = TRUE,
          fvf.nfeatures = 10000,
        )
    unpaired_rna@assays$RNA@var.features = features
    unpaired_rna <- ScaleData(unpaired_rna,features=features,split.by='group')
    
    # quantify gene activity
    gene.activities <- GeneActivity(unpaired_atac)
    unpaired_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
    # normalize gene activities
    DefaultAssay(unpaired_atac) <- "ACTIVITY"
    unpaired_atac <- NormalizeData(unpaired_atac)
    gene.activity.list <- SplitObject(unpaired_atac, split.by = "group")
        #select high variable features across samples 
    features2 = SelectIntegrationFeatures(
          gene.activity.list,
          nfeatures = 5000,
          verbose = TRUE,
          fvf.nfeatures = 10000,
    )
    unpaired_atac@assays$ACTIVITY@var.features = features2
    unpaired_atac <- ScaleData(unpaired_atac,split.by='group')

    gene_sel <- intersect(unpaired_rna@assays$RNA@var.features,unpaired_atac@assays$ACTIVITY@var.features)

    cca_res <- RunCCA(object1 = unpaired_rna, 
                      object2 = unpaired_atac,
                      assay1 = "RNA",
                      assay2 = "ACTIVITY",
                      num.cc = 30,
                      features = gene_sel,
                      renormalize = FALSE,
                      rescale = FALSE)
    cca_res <- RunUMAP(cca_res, dims = 1:30,reduction="cca")
    cell_loading <- cca_res@reductions$cca@cell.embeddings
    # cell_pair <- cell_pairing(cell_loading[colnames(unpaired_rna),],
    #                       cell_loading[colnames(unpaired_atac),])
    # RNA column is the list of RNA cell barcodes, while the ATAC column stores the ATAC profiles. The two lists are in the same order
    #colnames(cell_pair) <- c("RNA","ATAC")
    print("------ Saving integration result ------")
    # cca loadings
    df_umap = as.data.frame(cell_loading)
    # ===== added =====
    colnames(df_umap) = paste0("latent_",1:ncol(df_umap))
    df = cbind(df_umap,dataset=cca_res$technology)
    print(table(df$dataset))
    dir.create(file.path(out_dir,"rfigr"),recursive=TRUE)
    print("------ Saving integration and cell pairing result ------")
    write.csv(df,file.path(out_dir,"rfigr","rfigr_result.csv"))
    # # ASSUMING ONE-TO-ONE cell pairing 
    # cell_pair = cell_pair[which(!duplicated(cell_pair$RNA)),]
    # cell_pair = cell_pair[which(!duplicated(cell_pair$ATAC)),]
    # # saving cell pairs 
    # write.csv(cell_pair,file.path(out_dir,"rfigr","cell_pair.csv"))
    t2 <- Sys.time()
    dir.create(file.path(out_dir,"runtime"),recursive=TRUE)
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","rfigr_runtime.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Integration Done ------")
    print("------ No prediction ------")

    
}

print(paste0("argument 1: ",args[1]))
print(paste0("argument 2: ",args[2]))

if (length(args)<2) {
  stop("Insufficient number of arguments are supplied (input file).n", call.=FALSE)
}else if(length(args)==2) {
    run_rfigr_fn(args[1], args[2])
}else{
    print("More arguments than function needed are supplied, running function with the first two arugments")
    run_rfigr_fn(args[1], args[2])
}

