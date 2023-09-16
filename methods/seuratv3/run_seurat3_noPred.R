#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(Seurat)
require(Matrix)
source("seurat_wnn_project.R")
source("r_utils.R")
require(future)

run_seurat3_fn <- function(in_dir, out_dir){
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
    
    # Perform standard analysis of each modality independently RNA analysis
    unpaired_rna <- NormalizeData(unpaired_rna)
    unpaired_rna <- FindVariableFeatures(unpaired_rna)
    unpaired_rna <- ScaleData(unpaired_rna)
    unpaired_rna <- RunPCA(unpaired_rna)
    unpaired_rna <- RunUMAP(unpaired_rna, dims = 1:n_lat)
    
    # We exclude the first dimension as this is typically correlated with sequencing depth
    unpaired_atac <- RunTFIDF(unpaired_atac)
    unpaired_atac <- FindTopFeatures(unpaired_atac, min.cutoff = "q0")
    unpaired_atac <- RunSVD(unpaired_atac)
    unpaired_atac <- RunUMAP(unpaired_atac, reduction = "lsi", dims = 2:n_lat, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    
    # quantify gene activity
    gene.activities <- GeneActivity(unpaired_atac, features = VariableFeatures(unpaired_rna))

    # add gene activities as a new assay
    unpaired_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

    # normalize gene activities
    DefaultAssay(unpaired_atac) <- "ACTIVITY"
    unpaired_atac <- NormalizeData(unpaired_atac)
    unpaired_atac <- ScaleData(unpaired_atac, features = rownames(unpaired_atac))

    # Identify anchors
    transfer.anchors <- FindTransferAnchors(reference = unpaired_rna, 
                                            query = unpaired_atac, 
                                            features = VariableFeatures(object = unpaired_rna),
                                            reference.assay = "RNA", 
                                            query.assay = "ACTIVITY", 
                                            reduction = "cca")
    # note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
    # full transcriptome if we wanted to
    genes.use <- VariableFeatures(unpaired_rna)
    refdata <- GetAssayData(unpaired_rna, assay = "RNA", slot = "data")[genes.use, ]

    # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
    # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
    imputation <- TransferData(anchorset = transfer.anchors, 
                               refdata = refdata, 
                               weight.reduction = unpaired_atac[["lsi"]],
                               dims = 2:n_lat)
    unpaired_atac[["RNA"]] <- imputation

    coembed <- merge(x = unpaired_rna, y = unpaired_atac)

    # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
    # datasets
    coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
    coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
    coembed <- RunUMAP(coembed, dims = 1:n_lat)

    print("------ Saving integration result ------")
    df_umap = as.data.frame(coembed@reductions$pca@cell.embeddings[,1:n_lat])

    # ===== added =====
    colnames(df_umap) = paste0("latent_",1:ncol(df_umap))
    
    df = cbind(df_umap,dataset=coembed@meta.data$dataset)
    print(table(df_umap$dataset))
    dir.create(file.path(out_dir,"seurat3"),recursive=TRUE)
    write.csv(df,file.path(out_dir,"seurat3/","seurat3_result.csv"))
    t2 <- Sys.time()
    dir.create(file.path(out_dir,"runtime"),recursive=TRUE)
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","seurat3_runtime.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Done ------")
 
}


if (length(args)<2) {
  stop("Insufficient number of arguments are supplied (input file).n", call.=FALSE)
}else{ 
    print(paste0("Input directory: ",args[1]))
    print(paste0("Output directory: ",args[2]))
    run_seurat3_fn(args[1], args[2])
}
