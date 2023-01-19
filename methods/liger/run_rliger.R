#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(rliger)
require(Seurat)
require(stringr)
require(Signac)
require(Matrix)
source("r_utils.R")
#require(future)

# must have a fragment file named fragments.tsv.gz 
run_rliger_fn <- function(in_dir, out_dir, nclust=7){
    # starting time
    t1 <- Sys.time()
    n_lat = 30
    #plan("multisession")
    #options(future.rng.onMisue = "ignore")
    #print(paste0("workers used:",nbrOfWorkers()))
    # run_rliger: Liger object construction, processing, and integration 
    run_rliger <- function(unpaired_rna,unpaired_atac,nclust,k=20){
        gene.activities <- GeneActivity(unpaired_atac)
        unpaired_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

        data_comb <- list(atac = unpaired_atac[["ACTIVITY"]]@counts, rna = unpaired_rna[["RNA"]]@counts)
        data_comb <- createLiger(data_comb)

        data_comb <- rliger::normalize(data_comb)
        data_comb <- selectGenes(data_comb, datasets.use = 2)
        data_comb <- scaleNotCenter(data_comb)

        data_comb <- optimizeALS(data_comb, k = k)

        data_comb <- quantile_norm(data_comb)
        res <- find_resolution_liger(data_comb, nclust)
        data_comb <- louvainCluster(data_comb, resolution = res)
        
        return(data_comb)
    }

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
    
    # merging
    unpaired_rna <- merge(unpaired_rna,paired_rna)
    unpaired_atac <- merge(unpaired_atac,paired_atac)
    
    # run LIGER
    data_integrated <- run_rliger(unpaired_rna,unpaired_atac,nclust)
    # run UMAP
    data_integrated <- runUMAP(data_integrated, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
    
    df_umap <- data.frame(cbind(data_integrated@H.norm,predicted_ct=data_integrated@clusters))
    colnames(df_umap) [1:20] <- paste0("latent_",1:20)

    df_umap$dataset = "scRNA"
    df_umap[names(dataset_vec),"dataset"] = dataset_vec
    #df_umap = cbind(df_umap,data_integrated@tsne.coords)
    #colnames(df_umap)[c(ncol(df_umap)-1):ncol(df_umap)] = paste0("umap_",1:2)
    
    dir.create(file.path(out_dir,"rliger/"),recursive=TRUE)
    print("------ Saving integration result ------")
    write.csv(df_umap,file.path(out_dir,"rliger","rliger_result.csv"))
    t2 <- Sys.time()
    dir.create(file.path(out_dir,"runtime/"),recursive=TRUE)
    ## use '[[1]]' for clean output
    write.table(difftime(t2, t1, units = "secs")[[1]], 
                file = file.path(out_dir,"runtime","rliger_runtime.txt"), 
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    print("------ Done ------")
}


if (length(args)<2) {
  stop("Insufficient number of arguments are supplied (input file).n", call.=FALSE)
}else if(length(args)==2) {
    print(paste0("Input directory: ",args[1]))
    print(paste0("Output directory: ",args[2]))
    run_rliger_fn(args[1], args[2])
}else{
    print(paste0("Number of arguments: ",length(args)))
    print(paste0("Input directory: ",args[1]))
    print(paste0("Output directory: ",args[2]))
    print(paste0("Number of clusters: ",args[3]))
    run_rliger_fn(args[1], args[2],as.integer(args[3]))
}

    