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
    run_rliger <- function(data_comb,nclust,k=20,datasets.use=2){

        data_comb <- createLiger(data_comb)

        data_comb <- rliger::normalize(data_comb)
        # change this to use the first n urna datasets
        data_comb <- selectGenes(data_comb, datasets.use = datasets.use)
        data_comb <- scaleNotCenter(data_comb)

        data_comb <- optimizeALS(data_comb, k = k)

        data_comb <- quantile_norm(data_comb)
        res <- find_resolution_liger(data_comb, nclust)
        data_comb <- louvainCluster(data_comb, resolution = res)

        return(data_comb)
    }
    datasets = load_datasets(in_dir,obs=c("barcodes","batch"))
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
    
    unpaired_rna@meta.data$technology = "scRNA"
    unpaired_atac@meta.data$technology = "snATAC"
    paired_rna@meta.data$technology = "Multiome-RNA"
    paired_atac@meta.data$technology = "Multiome-ATAC"

    unpaired_rna@meta.data$group <- paste0(unpaired_rna$batch,"_",unpaired_rna$technology)
    paired_rna@meta.data$group <- paste0(paired_rna$batch,"_",paired_rna$technology)

    unpaired_atac@meta.data$group <- paste0(unpaired_atac$batch,"_",unpaired_atac$technology)
    paired_atac@meta.data$group <- paste0(paired_atac$batch,"_",paired_atac$technology)

    # merging
    unpaired_rna <- merge(unpaired_rna,paired_rna)
    unpaired_atac <- merge(unpaired_atac,paired_atac)

    urna_list = SplitObject(unpaired_rna, split.by = "group")
    
    gene.activities <- GeneActivity(unpaired_atac)
    unpaired_atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
    uatac_list = SplitObject(unpaired_atac, split.by = "group")

    urna_counts_list = lapply(urna_list,function(x){x[["RNA"]]@counts})
    uatac_counts_list = lapply(uatac_list,function(x){x[["ACTIVITY"]]@counts})
    data_comb<-c(urna_counts_list,uatac_counts_list)
    
    datasets.use = c(1:length(urna_counts_list))
    # TO-BE-EDITED, PASS a list of objects to liger.
    data_integrated <- run_rliger(data_comb,nclust,datasets.use=datasets.use)
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

    