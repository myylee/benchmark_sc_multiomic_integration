#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(stringr)
require(Seurat)
require(Signac)
require(Matrix)
require(bindSC)
source("r_utils.R")
require(future)

run_rbindsc_fn <- function(in_dir, out_dir){
    # starting time
    t1 <- Sys.time()
    n_lat = 30
    plan("multisession")
    options(future.rng.onMisue = "ignore")
    print(paste0("workers used:",nbrOfWorkers()))
    
    datasets = load_datasets(in_dir)
    paired_rna=datasets$paired_rna
    paired_atac=datasets$paired_atac
    unpaired_rna=datasets$unpaired_rna
    unpaired_atac=datasets$unpaired_atac

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

    # merging
    unpaired_rna <- merge(unpaired_rna,paired_rna)
    unpaired_atac <- merge(unpaired_atac,paired_atac)

    DefaultAssay(unpaired_rna) <- "RNA"
    unpaired_rna <- NormalizeData(unpaired_rna)
    unpaired_rna <- FindVariableFeatures(unpaired_rna, nfeatures = 5000)
    unpaired_rna <- ScaleData(unpaired_rna)
    unpaired_rna <- RunPCA(unpaired_rna)
    unpaired_rna <- FindNeighbors(unpaired_rna, dims = 1:20, reduction = "pca")
    unpaired_rna <- FindClusters(unpaired_rna, resolution = 0.5)
    unpaired_rna <- RunUMAP(unpaired_rna, reduction = "pca", dims = 1:15)

    # quantify gene activity
    gene.activities <- GeneActivity(unpaired_atac, features = VariableFeatures(unpaired_rna))
    # add gene activities as a new assay
    unpaired_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
    DefaultAssay(unpaired_atac) <- "ACTIVITY"
    unpaired_atac <- FindVariableFeatures(unpaired_atac, nfeatures = 5000)
    
    DefaultAssay(unpaired_atac) <- "ATAC"
    unpaired_atac <- RunTFIDF(unpaired_atac)
    unpaired_atac <- FindTopFeatures(unpaired_atac, min.cutoff = 50)
    unpaired_atac <- RunSVD(unpaired_atac)
    unpaired_atac <- FindNeighbors(unpaired_atac, dims = 1:20, reduction = "lsi")
    unpaired_atac <- FindClusters(unpaired_atac, resolution = 0.5)
    unpaired_atac <- RunUMAP(unpaired_atac, reduction = "lsi", dims = 1:15)

    DefaultAssay(unpaired_atac) <- "ACTIVITY"
    gene.use <- intersect(VariableFeatures(unpaired_rna), 
                          VariableFeatures(unpaired_atac))
    
    X <- unpaired_rna[["RNA"]][gene.use,]
    Y <- unpaired_atac[["ATAC"]][]
    Z0 <- unpaired_atac[["ACTIVITY"]][gene.use]
    type <- c(rep("RNA", ncol(X)), rep("ATAC", ncol(X)))

    a <- rowSums(Y)
    Y <- Y[a>50,]
    
    out <- dimReduce(dt1 =  X, dt2 = Z0,  K = 30)
    x <- out$dt1
    z0 <- out$dt2
    y  <- unpaired_atac@reductions$lsi@cell.embeddings

    res <- BiCCA( X = t(x) ,
                 Y = t(y), 
                 Z0 =t(z0), 
                 X.clst = unpaired_rna$seurat_clusters,
                 Y.clst = unpaired_atac$seurat_clusters,
                 alpha = 0.5, 
                 lambda = 0.5,
                 K = 15,
                 temp.path  = "out",
                 num.iteration = 50,
                 tolerance = 0.01,
                 save = TRUE,
                 parameter.optimize = FALSE,
                 block.size = 0)

    df_umap <- as.data.frame(rbind(res$u, res$r))
    colnames(df_umap) <- paste0("latent_",1:dim(df_umap)[2])
    df_umap$dataset = "scRNA"
    df_umap[names(dataset_vec),"dataset"] = dataset_vec
    print("------ Saving integration result ------")
    dir.create(file.path(out_dir,"rbindsc"),recursive=TRUE)
    write.csv(df_umap,file.path(out_dir,"rbindsc","rbindsc_result.csv"))
    t2 <- Sys.time()
    dir.create(file.path(out_dir,"runtime"),recursive=TRUE)
    
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","rbindsc_runtime.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Done ------")
    
}


if (length(args) < 2) {
  stop("Insufficient number of arguments are supplied (input file).n", call.=FALSE)
}

print(paste0("argument 1: ",args[1]))
print(paste0("argument 2: ",args[2]))

run_rbindsc_fn(args[1], args[2])
