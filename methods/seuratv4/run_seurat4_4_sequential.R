#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(Seurat)
require(Matrix)
source("seurat_wnn_project.R")
source("r_utils.R")
require(future)
require(patchwork)


# run_seurat4_4, perform sctransform per individual/donor, assuming 8GB of RAM
run_seurat4_fn <- function(in_dir, out_dir,nclust = 7){
    # starting time
    t1 <- Sys.time()
    
    #plan("multisession")
    #options(future.rng.onMisue = "ignore")
    #print(paste0("workers used:",nbrOfWorkers()))
    #options(future.globals.maxSize = 8000 * 1024^2)
    # load dataset 
    datasets = load_datasets(in_dir,obs=c("barcodes","batch"))
    pdata = datasets$paired_rna
    pdata[["ATAC"]] = datasets$paired_atac[["ATAC"]]
    DefaultAssay(pdata) <- "RNA"

    unpaired_rna = datasets$unpaired_rna
    unpaired_atac = datasets$unpaired_atac

    dir.create(file.path(out_dir,"seurat4int","figures"),recursive=TRUE)
    print("------ Processing Multiome ------")
    pdata <- process_paired_multidonor(pdata,
                            file_path = file.path(out_dir,"seurat4int","figures","paired_data_umap_seurat.pdf"),
                            nclust=nclust)
    unpaired_rna <- project_rna(pdata,unpaired_rna, 
                                file_path = file.path(out_dir,"seurat4int","figures","rna_data_umap_seurat.pdf"),
                                reference_assay = "rnaintegrated")

    print("------ Projecting snATAC ------")
    res <- project_atac_slsi(pdata,unpaired_atac, 
                             file_path = file.path(out_dir,"seurat4int","figures","atac_data_umap_seurat.pdf"),
                             return_anchor=T)
    # testing - set transfer of atac to based on their tutorial
#     res <- project_atac_ref_lsi(pdata,unpaired_atac, 
#                              file_path = file.path(out_dir,"seurat4int","figures","atac_data_umap_seurat.pdf"),
#                              return_anchor=T,reference_reduction = "integrated_lsi")
    unpaired_atac <- res[[1]]
    anchors_u2 <- res[[2]]

    dir.create(file.path(out_dir,"runtime/"),recursive=TRUE)
    print("------ Saving integration result ------")
    df_umap = do.call(rbind,list(pdata@reductions$wnn.umap@cell.embeddings,
                    unpaired_rna@reductions$ref.umap@cell.embeddings,
                    unpaired_atac@reductions$ref.umap@cell.embeddings))

    # ===== added =====
    colnames(df_umap) = paste0("latent_",1:ncol(df_umap))

    pdata@meta.data$predicted.ct = pdata$ct
    col_sel = "predicted.ct"
    df = data.frame(predicted_ct=c(pdata@meta.data[,col_sel],
           unpaired_rna@meta.data[,col_sel],
           unpaired_atac@meta.data[,col_sel]))
    # first df_umap is latent embedding, second df_umap is umap embedding, in this case, they are the same. 
    # During evaluation, if umap embedding is detected, no umap projection was run again
    df = cbind(df_umap,df_umap,df,
               dataset=c(rep("Multiome",ncol(pdata)),
                 rep("scRNA",ncol(unpaired_rna)),
                 rep("snATAC",ncol(unpaired_atac))))

    umap_col_range = c((ncol(df_umap)+1):(ncol(df_umap)*2))
    colnames(df)[umap_col_range] = paste0("umap_",1:length(umap_col_range))
    rownames(df) <- gsub("\\.","-",rownames(df))
    write.csv(df,file.path(out_dir,"seurat4int","seurat4int_result.csv"))
    t2 <- Sys.time()
    ## use '[[1]]' for clean output
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","seurat4int_runtime.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Done ------")
    
    print("------ Prediction ------")
    # starting time
    t1 <- Sys.time()

    # predict gene expression values
    rna <- TransferData(
      anchorset = anchors_u2,
      refdata = GetAssayData(pdata, assay = "SCT", slot = "data"),
      weight.reduction = "lsiproject"
    )

    # add predicted values as a new assay
    unpaired_atac[["SCT"]] <- rna
    
    #data_combined <- merge(pdata,unpaired_atac)
    
    write_mtx_folder(file.path(out_dir,"seurat4int","predicted","ATAC"),unpaired_atac,assay_key="ATAC",slot_key="counts","peak")
    write_mtx_folder(file.path(out_dir,"seurat4int","predicted","RNA"),unpaired_atac,assay_key="SCT",slot_key="data","gene")
    
    t2 <- Sys.time()
    ## use '[[1]]' for clean output
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","seurat4int_prediction_time.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Prediction Done ------")
    
}



if (length(args)<2) {
  stop("Insufficient number of arguments are supplied (input file).n", call.=FALSE)
}else if(length(args)==2) {
    print(paste0("Input directory: ",args[1]))
    print(paste0("Output directory: ",args[2]))
    run_seurat4_fn(args[1], args[2])
}else if(length(args) ==3){
    print(paste0("Input directory: ",args[1]))
    print(paste0("Output directory: ",args[2]))
    print(paste0("Number of clusters: ",args[3]))
    run_seurat4_fn(args[1], args[2],as.integer(args[3]))
}else{
    stop(paste0(length(args)," arguments are supplied (input file).n, please double check!"), call.=FALSE)
}

    
    