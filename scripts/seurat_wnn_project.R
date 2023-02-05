require(Seurat)
require(Signac)
require(dplyr)
require(ggplot2)
require(patchwork)


process_paired <- function(pdata, plot=T, file_path="plot.pdf",nclust=NULL,assay1="RNA",assay2="ATAC",cluster_algorithm = 2){
    source("r_utils.R")
    # RNA analysis
    DefaultAssay(pdata) <- assay1
    print("SCTransform")
    pdata <- SCTransform(pdata, verbose = T) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(pdata) <- assay2
    pdata <- RunTFIDF(pdata)
    pdata <- FindTopFeatures(pdata, min.cutoff = 5)
    pdata <- RunSVD(pdata)
    pdata <- RunUMAP(pdata, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    # WNN
    pdata <- FindMultiModalNeighbors(pdata, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    pdata <- RunUMAP(pdata, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",return.model=TRUE)
    # Supervised projection
    pdata <- RunSPCA(pdata, assay = 'SCT', graph = 'wsnn')
    n_clusters = nclust
    # louvain clutering by default
    res <- find_resolution(pdata, 
                           n_clusters = n_clusters,
                           algorithm = cluster_algorithm,
                           graph_name="wsnn")
    pdata <- FindClusters(pdata, graph.name = "wsnn", algorithm = cluster_algorithm, verbose = FALSE, resolution = res)

    if(plot){
        p1 <- DimPlot(pdata, reduction = "umap.rna",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
        p2 <- DimPlot(pdata, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
        p3 <- DimPlot(pdata, reduction = "wnn.umap",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
        pdf(file_path,width=8,height=11)
        plot(p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
        dev.off()
    }

    pdata@meta.data$ct <- as.character(Idents(pdata))

    return(pdata)
  
}



process_paired_multidonor <- function(pdata, plot=T, file_path="plot.pdf",nclust=NULL,assay1="RNA",assay2="ATAC",cluster_algorithm = 2, pergroup=NA){
    source("r_utils.R")
    # RNA analysis
    DefaultAssay(pdata) <- assay1

    # Split and RNA
    pdata.list <- SplitObject(pdata, split.by = "batch")
    pdata.list <- lapply(X = pdata.list, FUN = SCTransform, method = "glmGamPoi")
    features <- SelectIntegrationFeatures(object.list = pdata.list, nfeatures = 3000)
    pdata.list <- PrepSCTIntegration(object.list = pdata.list, anchor.features = features)
    pdata.list <- lapply(X = pdata.list, FUN = RunPCA, features = features)

    rna.anchors <- FindIntegrationAnchors(object.list = pdata.list, normalization.method = "SCT",
    anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)
    rna.combined.sct <- IntegrateData(anchorset = rna.anchors, normalization.method = "SCT", dims = 1:50,new.assay.name ='rnaintegrated')

    rna.combined.sct <- rna.combined.sct %>% 
        RunPCA() %>% 
        RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

    # UMAP based on integrated RNA
    p1 <- DimPlot(rna.combined.sct, group.by = "batch")+
        ggtitle("UMAP with integrated SCT datasets") + 
        theme(plot.title = element_text(hjust = 0.5))
    
    # ATAC
    atac.list <- SplitObject(rna.combined.sct, split.by = "batch")
    atac.list <- lapply(X = atac.list, FUN = function(obj){
        DefaultAssay(obj) <- "ATAC"
        obj <- FindTopFeatures(obj, min.cutoff = 5)
        obj <- RunTFIDF(obj)
        obj <- RunSVD(obj)
        obj@assays$RNA<-NULL
        obj@assays$SCT<-NULL
        obj@assays$rnaintegrated<-NULL
        obj
    })
    # merge
    atac.combined<- atac.list[[1]]
    for(i in 2:length(atac.list)){
        atac.combined <- merge(atac.combined,atac.list[[i]])
    }
    # process the combined dataset
    atac.combined <- FindTopFeatures(atac.combined, min.cutoff = 5)
    atac.combined <- RunTFIDF(atac.combined)
    atac.combined <- RunSVD(atac.combined)
    #atac.combined <- RunUMAP(atac.combined, reduction = "lsi", dims = 2:50)
    #p2 <- DimPlot(atac.combined, group.by = "batch",reduction = 'umap') + ggtitle("ATAC merged")

    # compute LSI
    # find integration anchors
    integration.anchors <- FindIntegrationAnchors(
      object.list = atac.list,
      anchor.features = rownames(atac.list[[1]]),
      reduction = "rlsi",
      dims = 2:50
    )

    # integrate LSI embeddings
    atac.integrated <- IntegrateEmbeddings(
      anchorset = integration.anchors,
      reductions = atac.combined[["lsi"]],
      new.reduction.name = "integrated_lsi",
      dims.to.integrate = 1:50
    )
    atac.integrated <- RunUMAP(atac.integrated, reduction = 'integrated_lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    p2 <- DimPlot(atac.integrated, group.by = "batch",reduction = 'umap.atac') + 
        ggtitle("ATAC integrated") +
        theme(plot.title = element_text(hjust = 0.5))
    
    integreated_multiome <- rna.combined.sct
    integreated_multiome[["ATAC"]]<- NULL
    # ATAC related counts should be from 'merge' workflow above
    integreated_multiome[["ATAC"]]<- atac.combined[['ATAC']]
    # Two reductions: integrated_lsi and umap.atac are added
    integreated_multiome@reductions<- append(integreated_multiome@reductions,atac.integrated@reductions)

    # WNN
    integreated_multiome <- FindMultiModalNeighbors(integreated_multiome, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:50))
    integreated_multiome <- RunUMAP(integreated_multiome, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",return.model=TRUE)
    # Supervised projection
    integreated_multiome <- RunSPCA(integreated_multiome, assay = 'rnaintegrated', graph = 'wsnn')
    n_clusters = nclust
    # louvain clutering by default
    res <- find_resolution(integreated_multiome, 
                           n_clusters = n_clusters,
                           algorithm = cluster_algorithm,
                           graph_name="wsnn")
    integreated_multiome <- FindClusters(integreated_multiome, graph.name = "wsnn", algorithm = cluster_algorithm, verbose = FALSE, resolution = res)

    if(plot){
        p3 <- DimPlot(integreated_multiome, reduction = "umap.rna",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p4 <- DimPlot(integreated_multiome, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p5 <- DimPlot(integreated_multiome, reduction = "wnn.umap",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p6 <- DimPlot(integreated_multiome, group.by = "batch", reduction = "umap.rna") + ggtitle("RNA") + theme(plot.title = element_text(hjust = 0.5))
        p7 <- DimPlot(integreated_multiome, group.by = "batch", reduction = "umap.atac") + ggtitle("ATAC") + theme(plot.title = element_text(hjust = 0.5))
        p8 <- DimPlot(integreated_multiome, group.by = "batch", reduction = "wnn.umap") + ggtitle("WNN") + theme(plot.title = element_text(hjust = 0.5))
        pdf(file_path,width=8,height=11)
        plot(wrap_plots(p1, p2, nrow = 2))
        plot(wrap_plots(p3, p4, p5, p6, p7, p8, nrow = 2))
        dev.off()
    }

    integreated_multiome@meta.data$ct <- as.character(Idents(integreated_multiome))

    return(integreated_multiome)
  
}


project_rna<-function(pdata,rna_data,plot=T,file_path="plot.pdf",return_anchor=F,reference_assay="SCT"){
    rna_data <- SCTransform(rna_data, verbose = FALSE,assay = 'RNA')
    DefaultAssay(rna_data) <- 'SCT'

    anchors_u1 <- FindTransferAnchors(
        reference = pdata,
        query = rna_data,
        normalization.method = "SCT",
        reference.assay = reference_assay,
        query.assay  = "SCT",
        reference.reduction = "spca",
        dims = 1:50
    )
  
    rna_data <- MapQuery(
        anchorset = anchors_u1, 
        query = rna_data,
        reference = pdata, 
        refdata = list(ct = "ct"),
        reference.reduction = "spca",
        reduction.model = "wnn.umap"
    )

    if(plot){
        p1 = DimPlot(rna_data, reduction = "ref.umap", group.by = "predicted.ct", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
        pdf(file_path,width=8,height=11)
        plot(p1)
        dev.off()
    }
    
    if(return_anchor){
        return(list(rna_data,anchors_u1))
    }else{
        return(rna_data)
    }

}

project_atac <- function(pdata,atac_data,plot=T,file_path="plot.pdf",return_anchor=F){
    
    DefaultAssay(pdata) <- 'SCT'
    gene.activities <- GeneActivity(atac_data, features = VariableFeatures(pdata))
    atac_data[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

    # normalize gene activities
    DefaultAssay(atac_data) <- "ACTIVITY"
    atac_data <- SCTransform(atac_data, verbose = FALSE,assay = 'ACTIVITY')

    anchors_u2 <- FindTransferAnchors(
        reference = pdata, 
        query = atac_data,
        normalization.method = 'SCT',
        features = VariableFeatures(object = pdata),
        reference.assay = "SCT", 
        query.assay ="ACTIVITY", 
        reduction = "cca"
    )

    atac_data <- MapQuery(
        anchorset = anchors_u2, 
        query = atac_data,
        reference = pdata, 
        refdata = list(ct = "ct"),
        reference.reduction = atac_data[["lsi"]],
        reduction.model = "wnn.umap",
        integrateembeddings.args = list(k.weight = 20)
    )
    if(plot){
        p1 <- DimPlot(atac_data, group.by = "predicted.ct", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
        pdf(file_path,width=8,height=11)
        plot(p1)
        dev.off()
    }

    if(return_anchor){
        return(list(atac_data,anchors_u2))
    }else{
        return(atac_data)
    }

  
}

project_atac_slsi <- function(pdata,atac_data,plot=T,file_path="plot.pdf",return_anchor=F,reference_assay="ATAC"){
    DefaultAssay(atac_data)<-"ATAC"
    atac_data <- RunTFIDF(atac_data)
    atac_data <- FindTopFeatures(atac_data, min.cutoff = 5)
    atac_data <- RunSVD(atac_data)
    
    DefaultAssay(pdata) <- 'ATAC'
    pdata <- RunSLSI(pdata, assay = reference_assay, graph = 'wsnn')

    anchors_u2 <- FindTransferAnchors(
        reference = pdata,
        query = atac_data,
        reference.assay = reference_assay,
        query.assay  = "ATAC",
        reference.reduction = "slsi",
        reduction = "lsiproject",
        dims = 2:50
    )


    atac_data <- MapQuery(
        anchorset = anchors_u2, 
        query = atac_data,
        reference = pdata, 
        refdata = list(ct = "ct"),
        reference.reduction = "slsi",
        reduction.model = "wnn.umap"
    )

  
    if(plot){
        p1 <- DimPlot(atac_data, group.by = "predicted.ct", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
        pdf(file_path,width=8,height=11)
        plot(p1)
        dev.off()
    }

    if(return_anchor){
        return(list(atac_data,anchors_u2))
    }else{
        return(atac_data)
    }

  
}

# for situations where there are multiple samples of ATAC in the reference, the reference-ATAC samples are first integrated among each other using rlsi, creating a reduction field called "integrated_lsi". Therefore, the reference_reduction is different, and the new data could be projected to lsi (instead of using slsi). Currently, the if reference_reduction = 'integrated_lsi' seem to be casuing error. Not sure why. 
project_atac_ref_lsi <- function(pdata,atac_data,plot=T,file_path="plot.pdf",return_anchor=F,reference_assay="ATAC",reference_reduction="lsi"){
    DefaultAssay(atac_data)<-"ATAC"
    atac_data <- RunTFIDF(atac_data)
    atac_data <- FindTopFeatures(atac_data, min.cutoff = 5)
    atac_data <- RunSVD(atac_data)
    
    DefaultAssay(pdata) <- reference_assay

    anchors_u2 <- FindTransferAnchors(
        reference = pdata,
        query = atac_data,
        reference.assay = reference_assay,
        query.assay  = "ATAC",
        reference.reduction = reference_reduction,
        reduction = "lsiproject",
        dims = 2:50
    )


    atac_data <- MapQuery(
        anchorset = anchors_u2, 
        query = atac_data,
        reference = pdata, 
        refdata = list(ct = "ct"),
        reference.reduction = reference_reduction,
        new.reduction.name = "ref.lsi",
        reduction.model = "wnn.umap"
    )

  
    if(plot){
        p1 <- DimPlot(atac_data, group.by = "predicted.ct", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
        pdf(file_path,width=8,height=11)
        plot(p1)
        dev.off()
    }

    if(return_anchor){
        return(list(atac_data,anchors_u2))
    }else{
        return(atac_data)
    }

  
}

