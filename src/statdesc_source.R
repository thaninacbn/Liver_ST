library(ggplot2)
library(ggpubr)
library(Seurat)
library(ggvenn)
library(ggrepel)


# Descriptive statistics plots  ------------------------------------------------


DescriptiveStatsPlots <- function(sample, verbose = FALSE, 
                              plot.detection.rate = TRUE,
                              plot.UMI = TRUE,
                              plot.QC = TRUE,
                              plot.scatter = TRUE,
                              mt_prefix= "^mt-"){
  
  if (!exists("Spatial", where = sample@assays)){
    stop("Data is not Spatial Transcriptomics Data.")
  }
  
  
  if (!exists("SCT", where = sample@assays)){
    message("SCT assay does not exist, normalizing the data using seurat::SCTransform")
    sample <- SCTransform(sample, assay = "Spatial")
  }
  
  message("Entry data is SCT normalized. Using SCT assay as default.")
  
  if (!exists("counts", where = sample@assays$SCT@SCTModel.list)){
    UpdateSCTAssays(sample)
    sample <- UpdateSeuratObject(sample)
    message("Seurat Object updated to V5")
  }

  
  if (!exists("percent.mt", where = sample@meta.data)){ 
    # Checks that the percentage of mt was included in metadata
    message("Percentage of Mitochondrial Genes is not included in metadata. Calculating.")
    #TODO : prévoir possiblement d'autres organismes (même si je suis pas sûre de quoi)
    sample <- PercentageFeatureSet(sample, mt_prefix, col.name = "percent.mt")
  }
  

  if(plot.QC) {
    p <- VlnPlot(sample, 
                 features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
                 ncol = 3) 
    p1 <- VlnPlot(sample, 
                  features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), 
                  ncol = 3) 
    #Stacks and prints the plots
    #combinevln <- ggarrange(p, p1, ncol = 1, nrow = 2)
    #print(combinevln)
    print(p)
    print(p1)

  }
  
  
  if (plot.scatter) {
    p1<-FeatureScatter(sample, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
    p2<-FeatureScatter(sample, feature1 = "nCount_SCT", feature2 = "nFeature_SCT")
    combine_scatter <- ggarrange(p1, p2, ncol = 2, nrow = 1)
    
    print(combine_scatter)
  }
  
  
  if (plot.detection.rate) {
    
    df_genes = sample@assays$SCT@SCTModel.list$counts@feature.attributes
    print(ggplot(df_genes,aes(x=detection_rate))+
            geom_histogram(aes(y = ..density..),bins=100)+
            ggtitle("Distribution of gene detection rate"))
  
  }
  
  if (plot.UMI){
    
    
    df_UMI = sample@assays$SCT@SCTModel.list$counts@cell.attributes
    print(ggplot(df_UMI,aes(x=umi))+
            geom_histogram(aes(y = ..density..),bins=100)+ 
            ggtitle("Distribution of UMI per spot"))
  }
  
  
  return(sample)
  
}

# Workaround for seurat::SpatiallyVariableFeatures -----------------------------

SpatiallyVariableFeatures_workaround <- function(object, assay="SCT", selection.method = "moransi") {
  #' This is work around function to replace SeuratObject::SpatiallyVariableFeatures function.
  #' return ranked list of Spatially Variable Features
  
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  
  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  
  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[
    data[[paste0(selection.method, ".spatially.variable")]] & 
      (!is.na(data[[paste0(selection.method, ".spatially.variable")]])), 
    moransi_cols
  ]
  
  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method,
                                                           ".spatially.variable.rank")]]), ]
  
  # Return row names of the sorted data frame
  rownames(sorted_data)
}


# Variable Genes computation ---------------------------------------------------

VariableGenes <- function(sample, # Seurat Object
                          hvg.selection = "vst", # Can  also be "mvp" or "disp"
                          nfeatures = 2000, #nb of features to select as top variable
                          top_nfeatures = 10, #nb of feature among those to show on plot
                          search.svg = TRUE, #whether to run FindSpatiallyVariableFeatures
                          svg.selection = "moransi") { #can also be "markvariogram"
  
  if (!exists("SCT", where = sample@assays)){
    stop("Data must be SCT normalized")
  }
  
  if (!exists("counts", where = sample@assays$SCT@SCTModel.list)){
    sample <- UpdateSCTAssays(sample)
    message("SCT Assay updated to V5")
  }
  
  sample <- FindVariableFeatures(sample,
                                 selection.method = hvg.selection, 
                                 nfeatures = nfeatures)
  sample[["HVG"]] <- VariableFeatures(sample)[1:ncol(sample)]
  
  top_n <- head(VariableFeatures(sample, top_nfeatures))
  hvg.plot <- VariableFeaturePlot(sample)
  hvg.plot<- LabelPoints(plot = hvg.plot, points = top_n, repel = TRUE)
  print(hvg.plot) 
  
  if (search.svg){
    sample <- FindSpatiallyVariableFeatures(sample,
                                            assay = "SCT",
                                            features = VariableFeatures(sample),
                                            selection.method = svg.selection)
   sample[["SVG"]] <- SpatiallyVariableFeatures_workaround(sample,
                                                           selection.method = svg.selection)[1:ncol(sample)]
  }
  
  return(sample)
}


# Run Dimension Reductions (PCA and UMAP) -------------------------------------

DimReduc <- function(sample, 
                     verbose = TRUE,
                     umap_dims = 1:5) {
  
  if (!exists("Spatial", where = sample@assays)){
    stop("Data is not Spatial Transcriptomics Data.")
  }
  
  if (!exists("SCT", where = sample@assays)){
    stop("Data must be SCT normalized.")
  }
  
  if (!exists("counts", where = sample@assays$SCT@SCTModel.list)){
    sample <- UpdateSCTAssays(sample)
    message("SCT Assay updated to V5")
  }
  
  #PCA
  sample <- RunPCA(sample, assay = "SCT", verbose = verbose)
  print(VizDimLoadings(sample, dims = 1:2, reduction = "pca",nfeatures=10))
  print(DimHeatmap(sample, dims = 1, balanced = TRUE))
  
  # UMAP
  sample <- RunUMAP(sample, reduction = "pca", dims = umap_dims, verbose = verbose) 
  sample<- ProjectDim(sample, reduction = "umap")#TODO ajouter les dims en option
  print(DimHeatmap(sample, reduction = "umap", dims = 1, projected = T, balanced = T))
  return(sample)
}



# Compare variable genes for multiple samples (still TODO) ---------------------

VariableGenesComparison <- function(sample,
                                    top_hvg = 20,
                                    top_svg = 20){
    
    if (!exists("HVG", where = sample@meta.data)){
      stop("HVG were not computed for sample")
    }
  
    
      if (!exists("SVG", where = sample@meta.data)){
        stop("SVG were not computed for sample")
      }
    
  list_genes <- list(HVG = sample[["HVG"]][1:top_hvg,], 
                     SVG = sample[["SVG"]][1:top_svg,])
  
  
  vennplot <- ggvenn(list_genes, c("HVG", "SVG"),
                    fill_color = c("Blue","Green"),
                    show_elements = TRUE,
                    label_sep = "\n",
                    set_name_size = 5,
                    set_name_color = c("Blue","Green"),
                    text_size = 3,
                    text_color = "black",
                    fill_alpha=.5) + 
    
    labs(title = "Intersection of Highly variable and Spatially variable genes") 
  
  print(vennplot)

}


GeneVisualization <- function(sample, 
                              gene,
                              feature.plot = TRUE,
                              reduction = "umap",
                              spatial.plot = TRUE,
                              alpha = 1) {
  
    if (feature.plot) {
      
      if (!exists(reduction, where = sample@reductions)){
        stop("Dimention Reduction was not computed for sample")
      }
      
      perc_spot = round((sum(GetAssayData(object = sample, layer = "data")[gene,]>0)/nrow(sample@meta.data))*100)
      subtitle = paste("Expressed in ",perc_spot, "% of spots", sep = "")
      title = paste(gene, " Expression")
      
      p1 <- FeaturePlot(sample, features=gene, reduction = reduction)+ 
        labs(Title = title, subtitle = subtitle)
      print(p1)
    }
  
  if (spatial.plot) {
    
    p2 <- SpatialFeaturePlot(sample, 
                             features= gene, 
                             slot = "scale.data", 
                             alpha = alpha)
    
    print(p2)
  }
  
  
}



# Testing ----------------------------------------------------------------------

testing = F

if (testing)
  
{
  test <- readRDS("data/clean/samp_1340_Normalized.rds")
  
  test <- DescriptiveStatsPlots(test)
  
  test <- VariableGenes(test)
  
  test <- DimReduc(test)
  
  VariableGenesComparison(test)
  
  GeneVisualization(test, "Glul")
  
  
  DimPlot(test, reduction = "umap")
  
}

#saveRDS(test, "data/clean/samp_1340_DimReduc.rds")

