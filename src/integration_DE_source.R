################################################################################
#                        Integration - DE source code                          #
################################################################################
library(Seurat)
library(scCustomize)
library(DESeq2)
library(citcdf)
library(dplyr)


#Function that takes 2 samples, merges them and performs UMI recorrection -------
MergeSamples <- function(sample1, sample2, 
                         sample_idents,#must be concatenation of idents e.g. c("CTRL", "STIM")
                         gene = "Alb"){
  
  sample1$spotID <- paste(sample1$Bio.Order, sample_idents[1], sep = "_")
  sample2$spotID <- paste(sample2$Bio.Order, sample_idents[2], sep = "_")
  
  samp12_merged <- merge(sample1, sample2, merge.data = TRUE, add.cell.id = sample_idents)
  
  ind1 = ncol(sample1)
  ind2 = ind1 + ncol(sample2)
  
  Idents(samp12_merged, cells = 1:ind1) <- sample_idents[1]
  Idents(samp12_merged, cells = (ind1+1):ncol(samp12_merged)) <- sample_idents[2]
  
  table(Idents(samp12_merged))
  samp12_merged[["ID"]] <- rev(Idents(samp12_merged))
  samp12_merged <- SCTransform(samp12_merged, assay = "Spatial", vst.flavor = "v1")
  samp12_merged <- RunPCA(samp12_merged)
  samp12_merged <- RunUMAP(samp12_merged, dims = 1:5)
  
  slot(object = samp12_merged@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Spatial"
  SCTResults(object=samp12_merged, slot="umi.assay")
  
  samp12_corrected <- PrepSCTFindMarkers(samp12_merged)
  
  condition_order <- sample_idents
  samp12_corrected[["ID"]] <- factor(x= samp12_corrected@meta.data$ID, levels = condition_order)
  samp12_merged[["ID"]] <- factor(x= samp12_merged@meta.data$ID, levels = condition_order)
  
  
  print(DimPlot(samp12_merged, reduction = "umap", group.by = "Bio.Order", split.by = "ID"))
  
  C1C2pal = c("#94d2bd","#0a9396")
  p1<- VlnPlot(samp12_merged,
               split.by = "ID",
               split.plot = T,
               group.by = "Bio.Order",
               features = gene,
               cols = C1C2pal,
               pt.size = 0) + ggtitle(paste0(gene, " sans correction"))
  p2 <- VlnPlot(samp12_corrected,
                split.by = "ID",
                split.plot = T,
                group.by = "Bio.Order",
                features = gene,
                cols = C1C2pal,
                pt.size = 0)+ ggtitle(paste0(gene,  " avec correction"))
  
  print(p1+p2)
  
  return(samp12_corrected)
}



#Function that subsamples a dataset and pseudobulks it --------------------------------------------
PreparePseudobulkSamples <- function(sample1, sample2, samp_id){
  
  subsamp1 <- c()
  subsamp2 <- c()
  subsamp3 <- c()
  
  for (i in levels(sample1@meta.data$Bio.Order)){
    tmp <- sample1@meta.data[which(sample1@meta.data$Bio.Order == i),]
    shuffled <- tmp[sample(1:nrow(tmp)), ] 
    #letters because r is annoying
    n <- nrow(tmp)/3
    m <- n+1
    o <- 2*n
    p <- o+1
    q <- 3*n
    
    subsamp1 <- c(subsamp1, rownames(shuffled[1:n,]))   
    subsamp2 <- c(subsamp2, rownames(shuffled[m:o,]))
    subsamp3 <- c(subsamp3, rownames(shuffled[p:q,]))   
  }
  
  
  subsample1 <- subset(sample1, cells = subsamp1)
  subsample2 <- subset(sample1, cells = subsamp2)
  subsample3 <- subset(sample1, cells = subsamp3)
  
  
  subsamp1 <- c()
  subsamp2 <- c()
  subsamp3 <- c()
  
  for (i in levels(sample2@meta.data$Bio.Order)){
    tmp <- sample2@meta.data[which(sample2@meta.data$Bio.Order == i),]
    shuffled <- tmp[sample(1:nrow(tmp)), ] 
    #inside variables because r is annoying
    n <- nrow(tmp)/3
    m <- n+1
    o <- 2*n
    p <- o+1
    q <- 3*n
    
    subsamp1 <- c(subsamp1, rownames(shuffled[1:n,]))   
    subsamp2 <- c(subsamp2, rownames(shuffled[m:o,]))
    subsamp3 <- c(subsamp3, rownames(shuffled[p:q,]))   
  }
  
  subsample4 <- subset(sample2, cells = subsamp1)
  subsample5 <- subset(sample2, cells = subsamp2)
  subsample6 <- subset(sample2, cells = subsamp3)
  
  object_list <- list(subsample1, subsample2, subsample3, subsample4, subsample5, subsample6)
  
  cell_ids <- c()
  for (i in 1:3){
    cell_ids <- c(cell_ids, paste0(samp_id[1], i))
  }
  for (i in 1:3){
    cell_ids <- c(cell_ids, paste0(samp_id[2], i))
  }
  
  merged_object <- Merge_Seurat_List(list_seurat = object_list, merge.data = T,
                                     add.cell.ids = cell_ids)
  
  ind1 = ncol(subsample1)
  ind2 = ind1 + ncol(subsample2)
  ind3 = ind2 + ncol(subsample3)
  ind4 = ind3 + ncol(subsample4)
  ind5 = ind4 + ncol(subsample5)
  ind6 = ind5 + ncol(subsample6)
  
  Idents(merged_object, cells = 1:ind1) <- cell_ids[1]
  Idents(merged_object, cells = (ind1+1):ind2) <- cell_ids[2]
  Idents(merged_object, cells = (ind2+1):ind3) <- cell_ids[3]
  Idents(merged_object, cells = (ind3+1):ind4) <- cell_ids[4]
  Idents(merged_object, cells = (ind4+1):ind5) <- cell_ids[5]
  Idents(merged_object, cells = (ind5+1):ind6) <- cell_ids[6]
  
  table(Idents(merged_object))
  merged_object[["ID"]] <- Idents(merged_object)
  
  merged_object <- SCTransform(merged_object, assay = "Spatial", vst.flavor = "v1")
  merged_object <- RunPCA(merged_object)
  merged_object <- RunUMAP(merged_object, dims = 1:5)
  
  print(DimPlot(merged_object, group.by = "ID"))
  print(DimPlot(merged_object, split.by = "ID", group.by = "Bio.Order"))
  
  
  pseudo_object<- AggregateExpression(merged_object,
                                      assays = "Spatial",
                                      return.seurat = T, 
                                      group.by = "ID")
  
  return(pseudo_object)
  
}


#Function that takes a pseudobulk object,
#the order in which the cells are and then performs DE -------------------------

PseudobulkDE <- function(pseudobulk_object, condition_order, res_df = T){
  
  
  pseudobulk_object[["condition"]]<- c(rep(condition_order[1], 3), rep(condition_order[2], 3))
  Idents(pseudobulk_object)<-pseudobulk_object@meta.data$condition
  
  counts <- as.matrix(GetAssayData(object = pseudobulk_object, layer = "counts"))
  metadata <- pseudobulk_object@meta.data
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = ~ condition)
  
  dds$condition <- relevel(dds$condition, ref = condition_order[2])
  
  dds <- DESeq(dds)
  cdt_name <- paste("condition_", condition_order[1], "_vs_",condition_order[2], sep = "")
  res <- results(dds, name = cdt_name, alpha = 0.05)
  
  if (res_df){
    return(data.frame(res))
  }
  else( return(res))
}


#Function that runs conditional independence testing ---------------------------

CitDE <- function(sample){
  
  Y <-data.frame(t(as.data.frame(LayerData(sample, assay = "SCT", layer = "data"))))
  
  X<-sample@meta.data$ID
  
  Z <- as.factor(sample@meta.data$Bio.Order)
  
  res <- do.call("rbind", pbapply::pblapply(1:ncol(Y), function(i) {             
    cit_asymp(Y[, i],data.frame(X=X), data.frame(Z=Z))} ))
  df<- data.frame(gene = colnames(Y), raw_pval = res$raw_pval, 
                  adj_pval = p.adjust(res$raw_pval,method = "BH"),
                  test_statistic = res$Stat)  
  
  return(df)
  
}


#Function that finds clusters by class -----------------------------------------

FindDEByClass <- function(merged_sample, ref, comp){
  
  Idents(merged_sample) <- "spotID"
  Idents(merged_sample) <- as.factor(Idents(merged_sample))
  
  final_df <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "v4_log2FC", "gene"))))
  
  merged_sample@meta.data$Bio.Order <- as.factor(merged_sample@meta.data$Bio.Order)
  
  
  
  list_cells = list()
  
  for (i in levels(merged_sample@meta.data$Bio.Order)){
    
    id2 <- as.character(paste(i, comp, sep="_"))
    id1 <- as.character(paste(i, ref, sep = "_"))
    
    print(class(id1))
   
      list_cells[[i]]$gp1 = which(merged_sample@meta.data$spotID == id1)
      list_cells[[i]]$gp2 = which(merged_sample@meta.data$spotID == id2)
  }
  
  print(levels(merged_sample@meta.data$Bio.Order))
   for (i in levels(merged_sample@meta.data$Bio.Order)){
    id1 <- paste(i, comp, sep="_")
    id2 <- paste(i, ref, sep = "_")
    tmp <- FindMarkers(merged_sample, assay="SCT",
                       ident.1 = id1,
                       ident.2 = id2,
                       test.use = "wilcox",
                       recorrect_umi = FALSE,
                       logfc.threshold = 0)
    tmp$cluster <- rep(i, nrow(tmp))
    
    current_cluster <- i
    genes <- rownames(tmp)
    
    # Extract gene expression data for the current gene and the cluster
    gp1 <-as.data.frame(merged_sample@assays$SCT$data[genes, list_cells[[current_cluster]]$gp1])
    gp1 <- expm1(gp1)+1
    gp1$mean <- rowMeans(gp1)
    gp1$logmean <- log(gp1$mean, base =2)
    
    gp2 <-as.data.frame(merged_sample@assays$SCT$data[genes, list_cells[[current_cluster]]$gp2])
    gp2 <- expm1(gp2)+1
    gp2$mean <- rowMeans(gp2)
    gp2$logmean <- log(gp2$mean, base =2)
    
    tmp$v4_log2FC <- gp2$logmean - gp1$logmean
    
    tmp$gene <- rownames(tmp)
    
    
    final_df <- rbind(final_df, tmp)
  
   }
  
  
  return(final_df)
}


#Function that takes a result table and compares it with the bulk table -------

CompareToBulk <- function(restable, 
                          bulktable,
                          logFC_col,
                          padj_col, 
                          FC_threshold = 1,
                          padj_threshold = 0.05){
  bulk_table <- bulktable[,c("GeneName", logFC_col, padj_col)]
  
  smolbulk <- bulk_table[which(abs(bulk_table[, logFC_col])>FC_threshold & bulk_table[, padj_col]<padj_threshold),]
  bulkgenes <- smolbulk$GeneName
  pseudobulk <- rownames(restable[which(abs(restable$log2FoldChange)>FC_threshold & restable$padj<padj_threshold),])
  
  list_DE <- list(bulk = bulkgenes,
                  pseudobulk = pseudobulk)
  
  print(ggVennDiagram(list_DE, color= "blue") + 
    scale_x_continuous(expand = expansion(mult = .2)))
  
  unfiltered_common <- base::intersect(bulk_table$GeneName, restable$X)
  pseudobulk_common <- restable[unfiltered_common, ]
  pseudobulk_common$GeneName <-rownames(pseudobulk_common)
  
  print(length(unfiltered_common))
  
  #Merging the dataframes
  smoler_bulk <- bulk_table %>%
    filter(GeneName %in% unfiltered_common)
  merged_df <- merge(smoler_bulk, pseudobulk_common, by = "GeneName")
  
  merged_df <- na.omit(merged_df)
  
  merged_df <- merged_df %>%
    mutate(
      condition = case_when(
        abs(merged_df[, logFC_col]) > FC_threshold & merged_df[, padj_col] < padj_threshold & abs(log2FoldChange) > FC_threshold & padj < padj_threshold ~ "both",
        abs(merged_df[, logFC_col]) > FC_threshold & merged_df[, padj_col] < padj_threshold ~ "bulk",
        abs(log2FoldChange) > FC_threshold & padj < padj_threshold ~ "pseudobulk",
        TRUE ~ "no" # if none of the conditions are met
      )
    )
  
  #Plotting the p-values
  p1 <- ggplot(merged_df, aes(x = merged_df[, padj_col], y = padj, color = condition, alpha = condition)) + geom_point()+
    xlab("pvalue ajustée bulk") +
    ylab("pvalue ajustée pseudobulk") +
    scale_alpha_manual(values = c("bulk" = 0.7, "pseudobulk" = 0.7, "both" = 0.8, "no"= 0.1))+
    geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
    geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
    xlim(0,0.5)+
    ylim(0,0.5)
  
  #Plotting the p-values
   p2 <- ggplot(merged_df, aes(x = merged_df[, logFC_col], y = log2FoldChange, color = condition, alpha = condition)) + geom_point()+
    xlab("Log2FC bulk") +
    ylab("log2FC pseudobulk") +
    scale_alpha_manual(values = c("bulk" = 0.6, "pseudobulk" = 0.6, "both" = 0.6, "no"= 0.4))
  
  print(p1)
  print(p2)
  
}


#Function that gets Pseudobulk DE genes ----------------------------------------
GetDEGenes <- function(restable, FC_threshold = 1,
                       padj_threshold = 0.05, up = T, down = T, filename = "DE_genes.txt"){
  
  if (!up){
    DE_genes <- rownames(restable[which(restable$log2FoldChange<FC_threshold & restable$padj<padj_threshold),])
  }
  
  if (!down){
    DE_genes <- rownames(restable[which(restable$log2FoldChange>FC_threshold & restable$padj<padj_threshold),])
  }
  
  if (up & down){
    DE_genes <- rownames(restable[which(abs(restable$log2FoldChange)>FC_threshold & restable$padj<padj_threshold),])
  }
  
  capture.output(cat(sapply(DE_genes, toString), sep="\n"), file = filename)
}


#Fonction qui de toute évidence ne MARCHE PAS AAARGGHHHH -----------------------
GetCommonDEGenes <- function(bulktable, restable,
                             logFC_col, padj_col, FC_threshold = 1,
                             padj_threshold = 0.05, up = T, down = T, pb_unique = F,
                             intersection = F,
                             filename = "DE_genes.txt"){
  bulk_table <- bulktable[,c("GeneName", logFC_col, padj_col)]
  rownames(restable)<- restable[,1]
  
  unfiltered_common <- base::intersect(bulk_table$GeneName, restable$X)
  pseudobulk_common <- restable[unfiltered_common, ]
  pseudobulk_common$GeneName <-rownames(pseudobulk_common)
  
  #Merging the dataframes
  smoler_bulk <- bulk_table %>%
    filter(GeneName %in% unfiltered_common)
  merged_df <- merge(smoler_bulk, pseudobulk_common, by = "GeneName")
  
  merged_df <- na.omit(merged_df)
  
  merged_df <- merged_df %>%
    mutate(
      condition = case_when(
        merged_df[, logFC_col] > FC_threshold & merged_df[, padj_col] < padj_threshold & log2FoldChange > FC_threshold & padj < padj_threshold ~ "both_up",
        merged_df[, logFC_col] > FC_threshold & merged_df[, padj_col] < padj_threshold ~ "bulk_up",
        log2FoldChange > FC_threshold & abs(merged_df[, logFC_col]) < FC_threshold & padj < padj_threshold  ~ "pseudobulk_up",
        merged_df[, logFC_col] < -FC_threshold & merged_df[, padj_col] < padj_threshold & log2FoldChange < -FC_threshold & padj < padj_threshold ~ "both_down",
        merged_df[, logFC_col] < -FC_threshold & merged_df[, padj_col] < padj_threshold ~ "bulk_down",
        log2FoldChange < -FC_threshold & abs(merged_df[, logFC_col]) < FC_threshold & padj < padj_threshold ~ "pseudobulk_down",
        TRUE ~ "no" # if none of the conditions are met
      )
    )
  
  merged_df$condition <- as.factor(merged_df$condition)
  
  
  if(up & down & !pb_unique){
    
    DE_genes <- merged_df[which(merged_df$condition == "pseudobulk_up" | merged_df$condition == "pseudobulk_down" | merged_df$condition == "both_up" | merged_df$condition == "both_down"),]$GeneName
    
   # DE_genes1 <- merged_df[which(merged_df$log2FoldChange > FC_threshold & merged_df[, logFC_col] > FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes2 <- merged_df[which(merged_df$log2FoldChange > FC_threshold & merged_df[, logFC_col] > FC_threshold & merged_df$padj < padj_threshold),]$GeneName
   # DE_genes3 <- merged_df[which(merged_df$log2FoldChange < FC_threshold & merged_df[, logFC_col] < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
   # DE_genes4 <- merged_df[which(merged_df$log2FoldChange < FC_threshold & merged_df[, logFC_col] < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes <- c(DE_genes1, DE_genes2, DE_genes3, DE_genes4)
  }
  
  if( up & down & pb_unique){
    
    DE_genes <- merged_df[which(merged_df$condition == "pseudobulk_up" | merged_df$condition == "pseudobulk_down"),]$GeneName
    
    #DE_genes1 <- merged_df[which(merged_df$log2FoldChange > FC_threshold & abs(merged_df[, logFC_col]) < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes2 <- merged_df[which(merged_df$log2FoldChange < FC_threshold & abs(merged_df[, logFC_col]) < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes <- c(DE_genes1, DE_genes2)
    }
  
  if (up & pb_unique){
    DE_genes <- merged_df[which(merged_df$condition == "pseudobulk_up" ),]$GeneName
    #DE_genes <- merged_df[which(merged_df$log2FoldChange > FC_threshold & abs(merged_df[, logFC_col]) < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
  }
  
  if (down & pb_unique ){
    DE_genes <- merged_df[which(merged_df$condition == "pseudobulk_down"),]$GeneName
    #DE_genes <- merged_df[which(merged_df$log2FoldChange < FC_threshold & abs(merged_df[, logFC_col]) < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
  }
  
  if (!down & !pb_unique ){
    
    DE_genes <- merged_df[which(merged_df$condition == "both_up" | merged_df$condition == "pseudobulk_up"),]$GeneName
    
    #DE_genes1 <- merged_df[which(merged_df$log2FoldChange > FC_threshold & merged_df[, logFC_col] > FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes2 <- merged_df[which(merged_df$log2FoldChange > FC_threshold & merged_df[, logFC_col] > FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes <- c(DE_genes1, DE_genes2)
  }
  
  if (!up & !pb_unique ){
    
    DE_genes <- merged_df[which(merged_df$condition == "both_down" | merged_df$condition == "pseudobulk_down"),]$GeneName
    
   # DE_genes1 <- merged_df[which(merged_df$log2FoldChange < FC_threshold & merged_df[, logFC_col] < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes2 <- merged_df[which(merged_df$log2FoldChange < FC_threshold & merged_df[, logFC_col] < FC_threshold & merged_df$padj < padj_threshold),]$GeneName
    #DE_genes <- c(DE_genes1, DE_genes2)
  }
  

  capture.output(cat(sapply(DE_genes, toString), sep="\n"), file = filename)
  


}

#testing zone waaahhhhhhhh -----------------------------


