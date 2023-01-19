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

options(future.globals.maxSize = 3000*1024^3)

run_rfigr_fn <- function(in_dir, out_dir){
    # starting time
    t1 <- Sys.time()
    #plan("multisession")
    #options(future.rng.onMisue = "ignore")
    print(paste0("number of cores available:",availableCores()))
    # load dataset 
    datasets = load_datasets(in_dir)
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
    
    unpaired_atac@meta.data$dataset <- "snATAC"
    unpaired_rna@meta.data$dataset <- "scRNA"
    paired_atac@meta.data$dataset <- "Multiome-ATAC"
    paired_rna@meta.data$dataset <- "Multiome-RNA"

    # merging
    unpaired_rna <- merge(unpaired_rna,paired_rna)
    unpaired_atac <- merge(unpaired_atac,paired_atac)

    unpaired_rna <- NormalizeData(unpaired_rna)
    unpaired_rna <- FindVariableFeatures(unpaired_rna,nfeatures = 5000)
    # We exclude the first dimension as this is typically correlated with sequencing depth
    unpaired_atac <- RunTFIDF(unpaired_atac)
    unpaired_atac <- FindTopFeatures(unpaired_atac, min.cutoff = "q0")
    unpaired_atac <- RunSVD(unpaired_atac)
    unpaired_atac <- RunUMAP(unpaired_atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    # quantify gene activity
    gene.activities <- GeneActivity(unpaired_atac)
    unpaired_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
    # normalize gene activities
    DefaultAssay(unpaired_atac) <- "ACTIVITY"
    unpaired_atac <- NormalizeData(unpaired_atac)
    unpaired_atac <- FindVariableFeatures(unpaired_atac, nfeatures = 5000)
    unpaired_atac <- ScaleData(unpaired_atac)

    gene_sel <- intersect(unpaired_rna@assays$RNA@var.features,unpaired_atac@assays$ACTIVITY@var.features)

    cca_res <- RunCCA(object1 = unpaired_rna, 
                      object2 = unpaired_atac,
                      assay1 = "RNA",
                      assay2 = "ACTIVITY",
                      num.cc = 30,
                      features = gene_sel,
                      renormalize = TRUE,
                      rescale = TRUE)
    cca_res <- RunUMAP(cca_res, dims = 1:30,reduction="cca")
    cell_loading <- cca_res@reductions$cca@cell.embeddings
    cell_pair <- cell_pairing(cell_loading[colnames(unpaired_rna),],
                          cell_loading[colnames(unpaired_atac),])
    # RNA column is the list of RNA cell barcodes, while the ATAC column stores the ATAC profiles. The two lists are in the same order
    colnames(cell_pair) <- c("RNA","ATAC")
    print("------ Saving integration result ------")
    # cca loadings
    df_umap = as.data.frame(cell_loading)
    # ===== added =====
    colnames(df_umap) = paste0("latent_",1:ncol(df_umap))
    df = cbind(df_umap,dataset=cca_res$dataset)
    print(table(df$dataset))
    dir.create(file.path(out_dir,"rfigr"),recursive=TRUE)
    print("------ Saving integration and cell pairing result ------")
    write.csv(df,file.path(out_dir,"rfigr","rfigr_result.csv"))
    # ASSUMING ONE-TO-ONE cell pairing 
    cell_pair = cell_pair[which(!duplicated(cell_pair$RNA)),]
    cell_pair = cell_pair[which(!duplicated(cell_pair$ATAC)),]
    # saving cell pairs 
    write.csv(cell_pair,file.path(out_dir,"rfigr","cell_pair.csv"))
    t2 <- Sys.time()
    dir.create(file.path(out_dir,"runtime"),recursive=TRUE)
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","rfigr_runtime.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Integration Done ------")
    print("------ Prediction ------")
    # starting time
    t1 <- Sys.time()
    # save prediction result 
    paired_rna <- subset(unpaired_rna, cells=cell_pair$RNA)
    paired_atac <- subset(unpaired_atac, cells=cell_pair$ATAC)

    # keep genes expressed by at least 3 cells. Set RNA cells to have ATAC cell barcodes
    rna_counts_mtx <- GetAssayData(object = paired_rna, slot = "counts",assay="RNA")
    colnames(rna_counts_mtx) <- colnames(paired_atac)
    md<-paired_atac@meta.data
    
    paired_comb <- CreateSeuratObject(counts = rna_counts_mtx,min.cells = 3,meta.data = md, assay = "RNA")
    
    atac_counts_mtx <- GetAssayData(object = paired_atac, slot = "counts",assay="ATAC")
    
    chrom_assay <- CreateChromatinAssay(
      counts = atac_counts_mtx,
      sep = c("-", "-"),
      min.cells = 1,
      min.features = 0
    )

    paired_comb[["ATAC"]] <- chrom_assay

    DefaultAssay(paired_comb)<-'RNA'
    paired_comb <- NormalizeData(paired_comb)
    
    write_mtx_folder(file.path(out_dir,"rfigr","predicted","ATAC"),paired_comb,assay_key="ATAC",slot_key="counts","peak")
    write_mtx_folder(file.path(out_dir,"rfigr","predicted","RNA"),paired_comb,assay_key="RNA",slot_key="data","gene")
    
    t2 <- Sys.time()
    ## use '[[1]]' for clean output
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","rfigr_prediction_time.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Prediction Done ------")
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

