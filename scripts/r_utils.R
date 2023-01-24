
require(EnsDb.Hsapiens.v86)

# dir_path: directory path
# seurat_data: seurat object 
# assay_key: an element in the seurat_data@assays
# slot_key: which data (e.g. counts, data, scaled.data) to save 
# bc: column name for the barcode column in meta.data 
# feature: string used to save feature output (e.g. gene, peak)
# mtx_name: mtx name for the counts matrix file
write_mtx_folder <- function(dir_path,seurat_data,assay_key,slot_key,feature,bc="barcodes",mtx_name="counts.mtx"){
    dir.create(dir_path,recursive=T)
    counts_mtx <- GetAssayData(object = seurat_data, slot = slot_key,assay=assay_key)
    # write matrix 
    writeMM(counts_mtx, file.path(dir_path,mtx_name))
    # write barcodes 
    write.table(seurat_data@meta.data[,bc],file.path(dir_path,"barcodes.tsv"),row.names = FALSE,col.names=FALSE,quote=FALSE)
    # write feature names
    write.table(rownames(counts_mtx),file.path(dir_path,paste0(feature,".tsv")),row.names = FALSE,col.names=FALSE,quote=FALSE)
}

read_mtx_folder <- function(dir_path,assay_key,var_list,obs_list,atac=FALSE,frag_path=NULL){
    files = list.files(dir_path)
    mtx_name = files[grepl(".mtx",files)]
    mtx <- readMM(file.path(dir_path,mtx_name))

    md_list <- list()
    for (i in obs_list){
        md_list[[i]] <- data.frame(read.csv(file.path(dir_path, files[grepl(i,files)]),header=FALSE))
    }
    md <- do.call(cbind, md_list)
    names(md) <- obs_list
    md

    features <- read.csv(file.path(dir_path, files[grepl(var_list[1],files)]),header=FALSE)
    rownames(mtx) <- features[,1]
    # parse peak naming to identify the two separators 
    peak_example =  features[,1][grep("chr",features[,1])[1]]
    seps = gsub("[0-9]","",gsub("chr[0-9]+","",peak_example))
    sep1 = substr(seps,1,1)
    sep2 = substr(seps,2,2)
  
    colnames(mtx) <- md$barcodes

    rownames(md) <- md$barcodes
    if(atac){
        #atac_ref = SeuratDisk::LoadH5Seurat(atac_seurat_ref_path)
        if(is.null(frag_path)){
            fragments = file.path(dir_path,"fragments.tsv.gz")
        }else{
            fragments = frag_path
        }
        chrom_assay <- CreateChromatinAssay(
            counts = mtx,
            sep = c(sep1,sep2),
            fragments = fragments,
            min.cells = 10
        )
        seurat_obj <- CreateSeuratObject(chrom_assay,assay=assay_key,meta.data=md)
        DefaultAssay(seurat_obj) <- "ATAC"
        annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
        ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
        seqlevels(annotation) <- ucsc.levels
        Annotation(seurat_obj) <- annotation
    }
    else{
        seurat_obj <- CreateSeuratObject(mtx,assay=assay_key,meta.data=md)
    }

    return(seurat_obj)
}


read_paired_data <- function(outer_path, frag_path=NULL){
    require(tidyverse)
    require(EnsDb.Hsapiens.v86)
    paired_rna <- read_mtx_folder(file.path(outer_path,"paired_RNA"),
                                 "RNA",c("gene"),c("barcodes"))
    paired_atac <- read_mtx_folder(file.path(outer_path,"paired_ATAC"),
                                  "ATAC",c("peak"),c("barcodes"))
    if(is.null(frag_path)){
        fragments = file.path(outer_path,"paired_RNA","fragments.tsv.gz")
    }else{
        fragments = frag_path
    }
    chrom_assay <- CreateChromatinAssay(
        counts = paired_atac@assays$ATAC@counts,
        sep = c(":", "-"),
        fragments = fragments,
        min.cells = 10
    )
    paired_rna[["ATAC"]] <-  chrom_assay
    
    DefaultAssay(paired_rna) <- "ATAC"
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
    seqlevels(annotation) <- ucsc.levels
    Annotation(paired_rna) <- annotation
    
    DefaultAssay(paired_rna) <- "RNA"
    
    return(paired_rna)
}


#@ bc: list(string) -  list of barcodes specifying cells to be used (should be available in the fragments.tsv
#@ ratio: [0,1] - percentage of reads in the fragment file to be kept [NOTE: if ratio = 1, no shuffling or downsampling is done. But a fragment file with the selected barcodes will be saved to the location.]
#@ frag_out: string - path to store the downsampled fragment file
#@ frag_path: string - input path to original fragment file. Default is the PBMC one stored on server
#@ pgk_source: string -  path to where tabix and bgzip are installed
#@ feature_path: string - input path to the feature
#@ ncore: int (default=4) - number of cores to be used [improving FeatureMatrix speed]
#@ cal_mat: boolean(default=TRUE) - if true, re-tabulate the peak matrix with the specified peaks and new downsampled fragments.tsv 

# ASSUME: 'tmp.csv' and 'tmp_bc.csv' are present to indicate what peak and cell should be counted 
downsampleFragments <- function(bc,ratio,frag_out,
                                frag_path = NULL,
                                pkg_source="/home/myylee/anaconda3/envs/scib2/bin/",
                                feature_path=NULL,
                                ncore=1, cal_mat=TRUE){
    require(dplyr)
    require(Signac)
    require(Matrix)
    require(future)
    #plan("multicore", workers = ncore)
    #plan()
    
    bc = unlist(bc)
    fragments = data.table::fread(frag_path)
    fragments_bc_sel <- fragments %>% filter(V4 %in% bc)
    downsample_atac <- function(fragments_df,ratio){
        target = round(dim(fragments_df)[1]*ratio)
        print(paste0("target number of reads in fragment is ", target))
        keep_idx = sample(1:dim(fragments_df)[1],target)
        return(fragments_df[keep_idx,])
    }
    

    fragments_bc_sel_down <- downsample_atac(fragments_bc_sel,ratio)
    orig_dir = getwd()
    setwd(frag_out)
    out_path = "fragments_unsorted.tsv"
    data.table::fwrite(fragments_bc_sel_down,out_path,sep='\t',col.names=FALSE)
    
    print(paste0("writing files in ",frag_out))
    print("saving result to fragment.tsv")
    system("sort -k 1,1 -k2,2n fragments_unsorted.tsv > fragments.tsv",intern=TRUE)
    print("ziping to fragment.tsv.gz")
    system2(paste0(pkg_source,"bgzip"), c("-f","fragments.tsv"))
    print("tabix fragment.tsv.gz")
    system2(paste0(pkg_source,"tabix"), c("-f -p","bed","fragments.tsv.gz"))
    print("remove fragments_unsorted.tsv")
    system("rm fragments_unsorted.tsv",intern=TRUE)

    setwd(orig_dir)

    if(cal_mat){
        fragments <- CreateFragmentObject(file.path(frag_out,"fragments.tsv.gz"))
        features = read.csv(feature_path,header=FALSE)
        print("StringToGRanges")
        features = StringToGRanges(features[,1])
        print("FeatureMatrix")
        peak_matrix_down = FeatureMatrix(
            fragments = fragments,
            features = features
        )
        print("select based on bc")
        print(length(bc))
        print(head(bc))
        print(head(colnames(peak_matrix_down)))
        bc_shared = length(intersect(colnames(peak_matrix_down),bc))
        print(paste0("Shared barcodes: ",bc_shared))
        peak_matrix_down = peak_matrix_down[,bc]
        print("done preparing peak_matrix_down")
        return(peak_matrix_down)
    }else{
        print("not recalculating peak matrix")
        return()
    }
    
}

# uses leiden by default
find_resolution <- function(seurat_data, n_clusters, seed = 0, algorithm=4, max_iter=50,graph_name=NULL,resolution_start=50){
  obtained_clusters = -1
  iteration = 0
  resolutions = c(0, resolution_start)
  while(obtained_clusters != n_clusters & iteration < max_iter){
    current_res = sum(resolutions)/2
    if(!is.null(graph_name)){
      seurat_data <- FindClusters(object = seurat_data, 
                                  verbose = FALSE, 
                                  algorithm = algorithm, 
                                  resolution = current_res,
                                  random.seed = seed,
                                  graph.name = graph_name)
    }else{
      seurat_data <- FindClusters(object = seurat_data, 
                                  verbose = FALSE, 
                                  algorithm = algorithm, 
                                  resolution = current_res,
                                  random.seed = seed)
    }

    labels = Idents(seurat_data)
    obtained_clusters = length(unique(labels))
    
    if (obtained_clusters > n_clusters){
      resolutions[2] = current_res
      print(resolutions)
    }else{
      resolutions[1] = current_res
    }
    iteration = iteration + 1
    print(paste0("iteration = ", iteration, "; current res = ",current_res, "; number of clusters = ",obtained_clusters))
  }
  return(current_res)
}

find_resolution_liger <- function(liger_data, n_clusters, max_iter=50,resolution_start=50){
    obtained_clusters = -1
    iteration = 0
    resolutions = c(0, resolution_start)
    while(obtained_clusters != n_clusters & iteration < max_iter){
        current_res = sum(resolutions)/2
 
        liger_data_clust <- louvainCluster(liger_data, resolution = current_res,random.seed=1234)

        labels = liger_data_clust@clusters
        obtained_clusters = length(unique(labels))

        if (obtained_clusters > n_clusters){
            resolutions[2] = current_res
            print(resolutions)
        }else{
            resolutions[1] = current_res
        }
        iteration = iteration + 1
        print(paste0("iteration = ", iteration, "; current res = ",current_res, "; number of clusters = ",obtained_clusters))
    }
    return(current_res)
}

# load data: nested folder system that contains paired_RNA, paired_ATAC, unpaired_RNA, unpaired_ATAC
# assumption: fragments.tsv.gz is present in in_dir/unpaired_ATAC and out_in_dirdir/paired_ATAC folder 
load_datasets <- function(in_dir,obs=c("barcodes")){

    print("------ Loading scRNA ------")
    unpaired_rna <- read_mtx_folder(file.path(in_dir,"unpaired_RNA"),
                                     "RNA",c("gene"),obs)
    print("------ Loading snATAC ------")
    frag_path_unpaired = file.path(in_dir,"unpaired_ATAC","fragments.tsv.gz")
    if(!file.exists(frag_path_unpaired)){
        stop("No path to fragment file found for unpaired_ATAC")
    }
    print(paste0("using fragment file for unpaired ATAC: ",frag_path_unpaired))
    unpaired_atac <- read_mtx_folder(file.path(in_dir,"unpaired_ATAC"),
                               "ATAC",c("peak"),obs,atac=TRUE,
                                frag_path = frag_path_unpaired)
    print("------ Loading multiome-RNA ------")
    paired_rna <- read_mtx_folder(file.path(in_dir,"paired_RNA"),
                                     "RNA",c("gene"),obs)
    print("------ Loading multiome-ATAC ------")
    frag_path_paired = file.path(in_dir,"paired_ATAC","fragments.tsv.gz")
    if(!file.exists(frag_path_paired)){
        stop("No path to fragment file found for paired_ATAC")
    }
    print(paste0("using fragment file for paired ATAC: ",frag_path_paired))
    paired_atac <- read_mtx_folder(file.path(in_dir,"paired_ATAC"),
                               "ATAC",c("peak"),obs,atac=TRUE,
                                frag_path = frag_path_paired)
    return(list(paired_rna=paired_rna,paired_atac=paired_atac,
                unpaired_rna=unpaired_rna, unpaired_atac=unpaired_atac))
}






