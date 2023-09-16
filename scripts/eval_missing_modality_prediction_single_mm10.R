#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require(Seurat)
require(BSgenome.Mmusculus.UCSC.mm10)
require(stringr)
require(Signac)
require(Matrix)
require(future)
require(dplyr)
source("r_utils.R")
require(future)
plan(multisession)
options(future.rng.onMisue = "ignore")
# 3GB
options(future.globals.maxSize = 3000*1024^3)
# loading saved result
# predict_link.R

# assumes that input RNA data has been normalized 
predict_links <- function(in_dir,out_dir,method_key,links_path,expression_slot="counts",predict_folder="predicted"){   
    require(BSgenome.Mmusculus.UCSC.mm10)
    require(EnsDb.Mmusculus.v79)
    require(future)
    plan(multisession)
    options(future.rng.onMisue = "ignore")
    nbrOfWorkers()
    
    print("reading ATAC")
    print(paste0("ATAC path: ",file.path(out_dir,method_key,predict_folder,"ATAC")))
    if(! file.exists(file.path(out_dir,method_key,predict_folder,"ATAC")) || ! file.exists(file.path(out_dir,method_key,predict_folder,"/RNA/"))){
        print("no imputation found, skip gene-peak association evaluation")
        return()
    }
    atac <- read_mtx_folder(file.path(out_dir,method_key,predict_folder,"/ATAC/"),
                        "ATAC",c("peak"),c("barcodes"),atac=TRUE,frag_path="")
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
    seqlevels(annotation) <- ucsc.levels
    Annotation(atac) <- annotation

    print("reading RNA")
    paired <- read_mtx_folder(file.path(out_dir,method_key,predict_folder,"/RNA/"),
                            "RNA",c("gene"),c("barcodes"),atac=FALSE)
    print("reading RNA done")
    paired[["ATAC"]] <- atac@assays$ATAC
    
    ## subset for cells from single-modality only ("snATAC")
    unpaired_atac_barcode <- read.csv(file.path(in_dir,"unpaired_ATAC","barcodes.tsv"),header = F)[,1]
    barcodes_predicted <- intersect(unpaired_atac_barcode,colnames(paired))
    paired <-paired[,barcodes_predicted]
    print(paste0("Number of cells with RNA predicted: ",ncol(paired)))
          
    links_truth <- read.csv(links_path)
    print("loaded links_truth")
    DefaultAssay(paired)<-"RNA"
    if(expression_slot=="data"){
        paired <- NormalizeData(paired)
    }else if (expression_slot != "counts"){
        stop("wrong expression.slot inputted")
    }
    
    DefaultAssay(paired)<-"ATAC"
    main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
    keep.peaks <- which(as.character(seqnames(granges(paired))) %in% main.chroms)
    paired[["ATAC"]] <- subset(paired[["ATAC"]], features = rownames(paired[["ATAC"]])[keep.peaks])

    paired <- RegionStats(paired, genome = BSgenome.Mmusculus.UCSC.mm10)

    # 50kb list, with all genes 
    system.time({
        paired_links <- LinkPeaks(
            object = paired,
            peak.assay = "ATAC",
            expression.assay = "RNA",
            peak.slot = "counts",
            expression.slot = expression_slot,
            genes.use = rownames(paired@assays$RNA),
            distance = 50000,
        )
    })
    gene_link_50kb<-as.data.frame(paired_links@assays$ATAC@links)

    gene_link_50kb_unique<- gene_link_50kb %>%
        arrange(peak, -pvalue) %>%
        dplyr::filter(duplicated(peak) == FALSE) 

    truth = paste0(links_truth$gene,"_",links_truth$peak)
    pred = paste0(gene_link_50kb_unique$gene,"_",gene_link_50kb_unique$peak)
    
    tp = length(intersect(truth,pred))
    fp = length(pred) - tp
    fn = length(truth) - tp
    precision = tp/(tp+fp) 
    recall = tp/(tp+fn)
    f1 = 2 * (precision * recall) / (precision + recall)
    percent_recovered_50kb = tp/dim(links_truth)[1]
    print(paste0("percent_recovered_50kb: ",percent_recovered_50kb))
    print(paste0("f1: ",f1))
    write.table(c("percent_recovered_50kb" = percent_recovered_50kb,"f1"=f1),
            file = file.path(out_dir,method_key,paste0(method_key,"_prediction_eval.txt")), 
            sep = "\t",
            col.names = FALSE)
    
}

if (length(args) < 4) {
  stop("Insufficient number of arguments are supplied (input file).n", call.=FALSE)
}

print(paste0("argument 1: ",args[1]))
print(paste0("argument 2: ",args[2]))
print(paste0("argument 3: ",args[3]))
print(paste0("argument 4: ",args[4]))

#links_path="dataset/multiome_pbmc_10k/pbmc_10x_pmat_sig_links_50kb_unique.csv"

predict_links(in_dir = args[1],
              out_dir = args[2],
              method_key = args[3],
              links_path = args[4])

print("----eval_missing_modality_prediction.R Done----")