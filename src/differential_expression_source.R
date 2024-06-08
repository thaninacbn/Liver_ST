################################################################################
#             Differential Expression Analysis - Source functions              #
################################################################################


#Function that merges two samples from two conditions
MergeData <- function(reference, treatment, refID, treatmentID, UMAP_dims = 1:5){
  
  #Assigns origin ID to each spot
  reference$spotID <- paste(reference$Bio.Order,refID, sep = "_")
  treatment$spotID <- paste(treatment$Bio.Order,treatmentID, sep = "_")
  Idents(reference) <- "spotID"
  Idents(treament) <- "spotID"
  
  IDs <- c(refID, treatmentID)
  
  #Merge data
  samp_merged <- merge(reference, treatment, merge.data = TRUE, add.cell.id = IDs)
  ind1 = ncol(reference)
  ind2 = ind1 + ncol(treatment)
  
  #Assign condition of origin to each spot
  Idents(samp_merged, cells = 1:ind1) <- refID
  Idents(samp_merged, cells = (ind1+1):ncol(samp_merged)) <- treatmentID
  
  #Creates a metadata slot for the condition
  samp_merged[["ID"]] <- rev(Idents(samp_merged))
  # Ren-normalizing (crucial: merging does not save the SCT layer and variable features??)
  samp_merged <- SCTransform(samp_merged, assay = "Spatial", vst.flavor = "v1")
  #PCA and UMAP on merged data
  samp_merged <- RunPCA(samp_merged)
  samp_merged <- RunUMAP(samp_merged, dims = UMAP_dims)
  
  #Reorder the leves according to the correct ID
  samp_merged[["ID"]] <- factor(x= samp_merged@meta.data$ID, levels = IDs)
  
  #Show dimplots of the merging. 
  #If samples dont superpose well (for example due to batch effects), user should
  #Integrate the data (CCA integration for example). Else, PrepSCTFindMarkers is sufficient
  print(DimPlot(samp_merged,
                reduction = "umap",
                group.by = "Bio.Order",
                split.by = "ID"))
  print(DimPlot(samp_merged,
                reduction = "umap",
                group.by = "ID"))
        
  return(samp_merged)  
}


#Function that corrects the data for library size -----------------------------
CorrectData <- function(merged_sample, refID, treatmentID, gene = "Alb"){
  
  #TODO : sanity check that this is indeed a merged sample ?
  
  #Correct for library size difference across samples
  samp_corrected <- PrepSCTFindMarkers(merged_sample)
  
  #reorder here as well
  condition_order <- c(refID, treatmentID)
  samp_corrected[["ID"]] <- factor(x= samp_corrected@meta.data$ID, levels = condition_order)
  
  #Compare before and after correction
  C1C2pal = c("#94d2bd","#0a9396")
  p1<- VlnPlot(merged_sample,
               split.by = "ID",
               split.plot = T,
               group.by = "Bio.Order",
               features = gene,
               cols = C1C2pal,
               pt.size = 0) + ggtitle(paste0(gene, " uncorrected"))
  p2 <- VlnPlot(samp_corrected,
                split.by = "ID",
                split.plot = T,
                group.by = "Bio.Order",
                features = gene,
                cols = C1C2pal,
                pt.size = 0)+ ggtitle(paste0(gene,  " corrected"))
  
  print(p1+p2)
  
  return(samp_corrected)
}






