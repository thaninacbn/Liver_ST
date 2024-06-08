################################################################################
#                       Marker Genes - Source functions                        #
################################################################################


library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggvenn)
library(patchwork)
library(gridExtra)


# Reorders UMAP classes according to bio-order --------------------------------

BioOrder <- function(sample, bio.order){
  # TODO: sécurité pour vérifier que StrongForms a bien été compute
  
  sample <- subset(sample, subset = StrongForms == "0", invert = T)
  
  ordered_level_list <- as.factor(bio.order)
  sample@meta.data[["Bio.Order"]] = sample@meta.data[["StrongForms"]]
  for (i in 1:length(ordered_level_list)) {
    sample@meta.data[["Bio.Order"]][which(sample@meta.data[["StrongForms"]] == i)] = ordered_level_list[i]
  }
  
  sample@meta.data$Bio.Order <- droplevels(sample@meta.data$Bio.Order)
  
  Idents(sample) <- sample@meta.data[["Bio.Order"]]
  return(sample)
  
  
}



#Find Markers with or without Roc AUC test-------------------------------------

FindStrongMarkers <- function(sample, 
                              strong.markers = T,
                              min_AUC = 0.7,
                              ntop = 20) {
  
  
  PreMarkers <- FindAllMarkers(sample, only.pos = T)
  classes <- levels(Idents(sample))
  temp.Markers <- subset(PreMarkers, p_val_adj < 0.05)
  
  v4_log2FC <- numeric(nrow(temp.Markers))
  
  list_cells = list()
  
  for (i in classes){
    list_cells[[i]]$gp1 = which(Idents(sample) == i)
    list_cells[[i]]$gp2 = which(Idents(sample) != i)
  }
  
  for (i in 1:nrow(temp.Markers)) {
    gene <- temp.Markers[i, 7]
    current_cluster <- temp.Markers[i, 6]
    
    # Extract gene expression data for the current gene and the cluster
    gp1 <-sample@assays$SCT$data[gene, list_cells[[current_cluster]]$gp1]
    exp_gp1 <- log(mean(expm1(gp1) + 1), base = 2)
    
    gp2 <- sample@assays$SCT$data[gene, list_cells[[current_cluster]]$gp2]
    exp_gp2 <- log(mean(expm1(gp2) + 1), base = 2)
    
    # Compute log2fc for this gene
    v4_log2FC[i] <- exp_gp1 - exp_gp2
  }
  
  temp.Markers$v4_log2FC <- v4_log2FC
  ref.Markers = subset(temp.Markers, 0.25 < avg_log2FC)
  
  # Add ROC test 
  if(strong.markers){
    
    genes_to_test <- unique(ref.Markers$gene)
    
    samp_levels <- levels(Idents(sample))
    
    pre.roc.markers <- FindAllMarkers(sample,
                                  verbose = F,
                                  test.use = "roc",
                                  features = genes_to_test)
    tmp.roc.markers <- subset(pre.roc.markers,
                              min_AUC < myAUC)
    
    
    res = data.frame()
    
    for (i in samp_levels){
      
      Wilcox.Markers <- ref.Markers  %>% group_by(cluster) %>% select(cluster,gene)  %>% filter(cluster==i)
      Roc.Markers <- tmp.roc.markers  %>% group_by(cluster) %>% select(cluster,gene) %>% filter(cluster==i)
      inter = as.data.frame(intersect(Wilcox.Markers, Roc.Markers)) 
      
      tmp.roc.markers = tmp.roc.markers[order(tmp.roc.markers$gene),]
      ref.Markers = ref.Markers[order(ref.Markers$gene),]
      
      temp = cbind("myAUC" = tmp.roc.markers[which(tmp.roc.markers$cluster == i & tmp.roc.markers$gene %in% inter$gene),]$myAUC, ref.Markers[which(ref.Markers$cluster == i & ref.Markers$gene %in% inter$gene),])
      res = rbind(res, temp)
      
    }
    
    res = res %>% arrange(cluster) %>% group_by(cluster) %>% top_n(n= ntop, wt = myAUC)
    
    ref.Markers <- res
    
  }
  
  return(ref.Markers)
  
}



# match cluster from markers, gaetan version ----------------------------------
MatchClustersFromMarkers <- function(Ref = Markers.Louvain,
                                     Comp = Markers.Leiden,
                                     Comp2 = NULL,
                                     Comp3 = NULL,
                                     ntop = 5,
                                     sorted = FALSE) {
  
  top5 <- Ref %>% group_by(cluster) %>% top_n(n = ntop, wt = v4_log2FC)
  NewOrder = Comp
  missed.value = NULL
  index.list = c() 
  if (!sorted) {
    for (i in levels(Ref$cluster)){
      topgene = top5[which(top5$cluster==i),]$gene
      tab = Comp[Comp$gene%in%topgene,]
      index = which.max(table(tab$cluster)) #vote majoritaire
      
      test.ambiguity = sort(table(tab$cluster), decreasing = TRUE)
      if(test.ambiguity[1] == test.ambiguity[2]) {
        print(table(tab$cluster))
        print(paste("indice initial :", index, "pas modifié pour", i))
        print("Vote ambigu")
        missed.value = i 
      }
      else {
        NewOrder[which(Comp$cluster==index),]$cluster = i
        index.list = append(index.list, value = index)
        print(table(tab$cluster))
        print(paste("indice initial :", index, "modifié pour", i))
      }
      
    }
    if (!is.null(missed.value)) {
      missing.index = setdiff(1:max(levels(Ref$cluster)), index.list)[1]
      print(paste("indice manquant :", missing.index, "remplacé par", missed.value))
      NewOrder[which(Comp$cluster==missing.index),]$cluster = missed.value
      
    }
  }
  argnames <- sys.call()
  arg.name.list = unlist(lapply(argnames[-1], as.character))
  
  if(!is.null(Comp2)){
    NewOrder2 = Comp2
    if (!sorted) {
      for (i in levels(Ref$cluster)){
        topgene = top5[which(top5$cluster==i),]$gene
        tab = Comp2[Comp2$gene%in%topgene,]
        index = which.max(table(tab$cluster)) #vote majoritaire
        NewOrder2[which(Comp2$cluster==index),]$cluster = i
        #print(table(tab$cluster))
        #print(paste("indice initial :", index, "modifié pour", i )
      }
    }
    #si il y a 4 éléments à comparer
    if(!is.null(Comp3)){
      NewOrder3 = Comp3
      if (!sorted) {
        for (i in levels(Ref$cluster)){
          topgene = top5[which(top5$cluster==i),]$gene
          tab = Comp3[Comp3$gene%in%topgene,]
          index = which.max(table(tab$cluster)) #vote majoritaire
          NewOrder3[which(Comp3$cluster==index),]$cluster = i
        }
      }
      for (i in levels(Ref$cluster)){
        VennAlgo = list("Ref"= Ref[which(Ref$cluster==i),]$gene,
                        "Comp1" = NewOrder[which(NewOrder$cluster==i),]$gene,
                        "Comp2" = NewOrder2[which(NewOrder2$cluster==i),]$gene,
                        "Comp3" = NewOrder3[which(NewOrder3$cluster==i),]$gene)
        
        mytable<-cbind(c("Ref","Comp1","Comp2"),c(arg.name.list[1],arg.name.list[2],arg.name.list[3]))
        
        vennplot <- ggvenn(VennAlgo, c("Ref", "Comp1", "Comp2", "Comp3"), fill_color = c("Blue","Green", "Yellow","Red"), show_elements = TRUE, label_sep = "\n",  set_name_size = 5, set_name_color = c("Blue","Green","Yellow", "Red"), text_size = 3, fill_alpha=.5) + labs(title=paste0("Intersection des markers genes dans le cluster ",i), caption = paste0("ref : ", arg.name.list[1], "   comp1 : ", arg.name.list[2], "    comp2 : ", arg.name.list[3], "     comp3 : ", arg.name.list[4] )) + theme(plot.caption.position = "plot",plot.caption = element_text(hjust = 0.5))
        
        print(vennplot)
      }
    }  
    
    else if (is.null(Comp3)) {
      for (i in levels(Ref$cluster)){
        VennAlgo = list("Ref"= Ref[which(Ref$cluster==i),]$gene,
                        "Comp1" = NewOrder[which(NewOrder$cluster==i),]$gene,
                        "Comp2" = NewOrder2[which(NewOrder2$cluster==i),]$gene)
        
        vennplot <- ggvenn(VennAlgo, c("Ref", "Comp1", "Comp2"), fill_color = c("Blue","Green", "Orange"), show_elements = TRUE, label_sep = "\n",  set_name_size = 5, set_name_color = c("Blue","Green","Orange"), text_size = 3, text_color = "black", fill_alpha=.5) + theme(plot.caption.position = "plot",plot.caption = element_text(hjust = 0.5)) + labs(title = paste0("Intersection des markers genes dans le cluster ",i), caption = paste0("ref : ", arg.name.list[1], "   comp1 : ", arg.name.list[2], "    comp2 : ", arg.name.list[3] ))  
        print(vennplot) }
    }
  }
  
  
  else if (is.null(Comp2) & is.null(Comp3)){
    for (i in levels(Ref$cluster)){
      VennAlgo = list("Ref"= Ref[which(Ref$cluster==i),]$gene,
                      "Comp" = NewOrder[which(NewOrder$cluster==i),]$gene)
      vennplot <- ggvenn(VennAlgo, c("Ref", "Comp"), fill_color = c("Blue","Green"), show_elements = TRUE, label_sep = "\n",  set_name_size = 5, set_name_color = c("Blue","Green"), text_size = 3, text_color = "black", fill_alpha=.5)  + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + labs(title=paste0("Intersection des markers genes dans le cluster ",i), caption = paste0("ref : ", arg.name.list[1], "   comp : ", arg.name.list[2]))  
      print(vennplot)
      
    }
    
  }
  
}
  

#gaetan function: reorders clusters from marker genes
ReorderClustersFromMarkers <- function(Ref.cluster,
                                       Comp.cluster,
                                       Ref.Markers = NULL,
                                       Comp.Markers = NULL,
                                       Ref.ech = NULL,
                                       Comp.ech = NULL,
                                       ntop = 10) {
  
  
  top5 <- Ref.Markers %>% group_by(cluster) %>% top_n(n = ntop, wt = avg_log2FC)
  NewOrder = Comp.cluster
  missed.value = NULL
  index.list = c() 
  
  if(max(levels(Ref.Markers$cluster)) != max(levels(Comp.Markers$cluster))) {
    print("Clusterings do not have the same length, alorithm may crash.")
  }
  
  for (i in levels(Ref.Markers$cluster)){
    topgene = top5[which(top5$cluster==i),]$gene
    tab = Comp.Markers[Comp.Markers$gene%in%topgene,]
    index = which.max(table(tab$cluster)) #vote majoritaire
    
    test.ambiguity = sort(table(tab$cluster), decreasing = T)
    if(test.ambiguity[1] == test.ambiguity[2]){
      missed.value = i
    }
    
    else {
      NewOrder[which(Comp.cluster == index)] = i
      print(paste("initial index:", index, "Changed for", i ))
      index.list = append(index.list, value = index)
    }
  }
  
  if (!is.null(missed.value)){
    missing.index = setdiff(1:max(levels(Ref.Markers$cluster)), index.list)[1]
    NewOrder[which(Comp.cluster == missing.index)] = missed.value
  }
  return(NewOrder)
}
  

# Function that compares the classifications for n samples -------------------

CompareConditions <- function(sample_list){
  for (sample in sample_list){
    #Removes level 0 in BioOrder bc we subsampled it
    sample@meta.data$Bio.Order <- droplevels(sample@meta.data$Bio.Order)
    #Sets BioOrder as Ident in case it's not already
    sample <- SetIdent(sample, value = sample@meta.data$Bio.Order)
  }
  
  compare_df = data.frame()
  
  for (i in 1:length(sample_list)){
    sample_idents = Idents(sample_list[[i]])
    tmp_df = as.data.frame(table(Idents(sample_list[[i]])))
    print(tmp_df)
    tmp_df = cbind(tmp_df,
                   "Sample"= rep(i, length(levels(sample_idents))),
                   "Prop" = round(100*tmp_df$Freq/ncol(sample_list[[i]])))
    compare_df = rbind(compare_df, tmp_df)
  }
  colnames(compare_df)[1] = "ClustID"
  compare_df$Sample = as.factor(compare_df$Sample)
  
  p1 <- ggplot(compare_df)+
    aes(x=ClustID, fill=Sample) +
    geom_bar(stat="identity", aes(y=Prop), position=position_dodge())  +
    geom_text(aes(label=Freq, y = (Prop-1)), vjust=1.6, color="black",
              position = position_dodge(0.9), size=3)+
    scale_fill_brewer(palette = "Reds")+
    labs(y = "Distribution of spots in classes", x = "Clusters", colour = "Samples")+
    theme_minimal()
  
  print(p1)
}

# Function that makes a pretti heatmap -----------------------------------------

ClusterHeatmap <- function(sample, strong_markers, ordered_list = NULL, return_list = F){
  
  sample@meta.data$Bio.Order <- droplevels(sample@meta.data$Bio.Order)
  
  #Scale the data (not sure why we have to do it but anyways)
  scaled = ScaleData(sample,
                     features = rownames(sample),
                     do.scale = T,
                     do.center = T,
                     model.use = "negbinom")
  
  
  for (i in levels(scaled@meta.data$Bio.Order)) {
    meanval = rowMeans(scaled@assays$SCT@scale.data[, which(scaled@meta.data$Bio.Order==i)])
    scaled@assays$SCT@scale.data[, which(scaled@meta.data$Bio.Order==i)] = meanval
  }
  
  if (is.null(ordered_list)){
    
    ordered_list = list()
    
    mid = length(levels(scaled@meta.data$Bio.Order))-2
    
    #Orders the genes in highly -> lowly expressed order
    for (k in 1:mid) {
      list <- strong_markers[strong_markers$cluster == k,]$gene
      data.gene <- scaled@assays$SCT@scale.data[list, which(scaled$Bio.Order == k)]
      ordered_list[[k]] <- names(sort(rowMeans(data.gene), decreasing = T))
    }
    
    #class mid+1 ordered in highly -> lowly but without the genes in mid+2
    list_mid <- strong_markers[strong_markers$cluster == (mid+1),]$gene
    list_last <- strong_markers[strong_markers$cluster == (mid+2),]$gene
    
    kept_mid <- setdiff(list_mid, list_last)
    #updates list for mid+1
    data.gene <- scaled@assays$SCT@scale.data[kept_mid, which(scaled$Bio.Order == (mid+1))]
    ordered_list[[mid+1]] <- names(sort(rowMeans(data.gene), decreasing = T))
  
    #Except for the last class where the gens are ordered in lowly -> highly
    #This is so we have a gradient as nice as possible !!
    data.gene <- scaled@assays$SCT@scale.data[list_last, which(scaled$Bio.Order == (mid+2))]
    ordered_list[[mid+2]] <- names(sort(rowMeans(data.gene)))
    fully_ordered_list <- unique(unlist(ordered_list))
  }
  
  else {
    fully_ordered_list = ordered_list
  }
  
  #Create smth that generates a color grad vector depending on the nb of classes OR 
  # take color grad as input
  color_grad = rev(c("#ce1212","#f57a00","#e6d11d","#6FD577","#048FEB","#240CE1"))
  heatmap <- DoHeatmap(scaled,
            features = fully_ordered_list,
            group.by = "Bio.Order",
            group.colors = color_grad,
            disp.min = -2.5,
            disp.max = 2.5 ) + scale_fill_viridis_c(option = "magma")
  
  print(heatmap)
  #print(mid)
  
  if (return_list) {
    return(fully_ordered_list)
  }
  
}

# Function that plots the zonation of a gene -----------------------------------

GeneZonation <- function(sample_list, gene){
  
  VlnPlotTheme <- theme(legend.position = "none",
                        plot.subtitle = element_text(hjust = 0.5, face = "italic"))
  list_plot = list()
  
  for (i in 1: length(sample_list)) {
    p1 <- VlnPlot(sample_list[[i]], 
                  features = gene,
                  slot = "data",
                  pt.size = 0,
                  assay = "SCT") +
      ylim(0, 5) +
      VlnPlotTheme +
      ggtitle(paste(gene, " in condition ", i, sep = ""))+
      labs(x = "Zonation")
    
    list_plot[[i]]= p1
  }
  
  
p2 = wrap_plots(list_plot, ncol = 2) + 
  ggtitle(paste("Expression of gene ", gene, sep = ""))
  
print(p2)
  #return(list_plot)
  
}




# testing ----------------------------------------------------------------------

testing = F

if (testing) {
  
  cond1 <- readRDS("data/clean/samp_1340_Clustering_k6.Rds")
  cond2 <- readRDS("data/clean/samp_1352_Clustering_k6.Rds")
  cond3 <- readRDS("data/clean/samp_1353_Clustering_k6.Rds")
  cond4 <- readRDS("data/clean/samp_1354_Clustering_gamma15_k6.Rds")
  
  DimPlot(cond4, group.by = "StrongForms")
  
  cond1 <- BioOrder(cond1, c(3,6,5,2,1,4))
  
  DimPlot(cond1)
  
  #cond2 <- subset(cond2, subset = StrongForms == "0", invert = T)
  #cond3 <- subset(cond3, subset = StrongForms == "0", invert = T)
  #cond4 <- subset(cond4, subset = StrongForms == "0", invert = T)
  
 # Idents(cond2) <- cond2@meta.data$StrongForms
  #Idents(cond3) <- cond3@meta.data$StrongForms
  
  DimPlot(cond2, group.by = "StrongForms")
  cond2 <- BioOrder(cond2, c(3,2,1,4,6,5))
  #DimPlot(cond2)
  
  cond3 <- BioOrder(cond3, c(2,5,4,1,3,6))
  DimPlot(cond3)
  
  #cond4 <- BioOrder(cond4, c(3,5,6,2,4,1))
  #DimPlot(cond4)
  
  system.time(marks1 <- FindStrongMarkers(cond1, strong.markers = F))
  
  marks2 <- FindStrongMarkers(cond2, strong.markers = F)
  marks3 <- FindStrongMarkers(cond3, strong.markers = F)
  marks4 <- FindStrongMarkers(cond4, strong.markers = F)
  
  nb = 10
  
  top1 <- marks1 %>% group_by(cluster) %>% top_n(n = nb, wt = v4_log2FC)
  top2 <- marks2 %>% group_by(cluster) %>% top_n(n = nb, wt = v4_log2FC)
  top3 <- marks3 %>% group_by(cluster) %>% top_n(n = nb, wt = v4_log2FC)
  
  MatchClustersFromMarkers(Ref = top1, Comp = top3, ntop=nb)
  
  cond2@meta.data$Bio.Order <- ReorderClustersFromMarkers(Ref.cluster = cond1@meta.data$Bio.Order,
                                                          Comp.cluster = cond2@meta.data$StrongForms,
                                                          Ref.Markers = stmarks1,
                                                          Comp.Markers = stmarks2,
                                                          Ref.ech = cond1,
                                                          Comp.ech = cond2,
                                                          ntop = 10)
  
  
  DimPlot(cond2, group.by = "Bio.Order")
  
  cond3@meta.data$Bio.Order <- ReorderClustersFromMarkers(Ref.cluster = cond1@meta.data$Bio.Order,
                                                          Comp.cluster = cond3@meta.data$StrongForms,
                                                          Ref.Markers = marks1,
                                                          Comp.Markers = marks3,
                                                          Ref.ech = cond1,
                                                          Comp.ech = cond3,
                                                          ntop = 10)
  
  DimPlot(cond3, group.by = "Bio.Order")
  
  cond4@meta.data$Bio.Order <- ReorderClustersFromMarkers(Ref.cluster = cond1@meta.data$Bio.Order,
                                                          Comp.cluster = cond4@meta.data$StrongForms,
                                                          Ref.Markers = marks1,
                                                          Comp.Markers = marks4,
                                                          Ref.ech = cond1,
                                                          Comp.ech = cond4,
                                                          ntop = 5)
  
  DimPlot(cond4, group.by = "Bio.Order")
  
  stmarks1 <- FindStrongMarkers(cond1)
  stmarks2 <- FindStrongMarkers(cond2)
  stmarks3 <- FindStrongMarkers(cond3)
  stmarks4 <- FindStrongMarkers(cond4)
  
  test <- c(cond1, cond2, cond3, cond4)
  df <- CompareConditions(test)
  
  GeneZonation(test, "Mup11")
  
  cond2@meta.data$Bio.Order <- ReorderClustersFromMarkers(Ref.cluster = cond1@meta.data$Bio.Order,
                                                          Comp.cluster = cond2@meta.data$StrongForms,
                                                          Ref.Markers = stmarks1,
                                                          Comp.Markers = stmarks2,
                                                          Ref.ech = cond1,
                                                          Comp.ech = cond2,
                                                          ntop = 10)
  
  DimPlot(cond2, group.by = "Bio.Order")
  
  
  cond3@meta.data$Bio.Order <- ReorderClustersFromMarkers(Ref.cluster = cond1@meta.data$Bio.Order,
                                                          Comp.cluster = cond3@meta.data$StrongForms,
                                                          Ref.Markers = stmark1,
                                                          Comp.Markers = stmarks3,
                                                          Ref.ech = cond1,
                                                          Comp.ech = cond3,
                                                          ntop = 5)
  
  DimPlot(cond3, group.by = "Bio.Order")
  
  cond4@meta.data$Bio.Order <- ReorderClustersFromMarkers(Ref.cluster = cond1@meta.data$Bio.Order,
                                                          Comp.cluster = cond4@meta.data$StrongForms,
                                                          Ref.Markers = marks1,
                                                          Comp.Markers = stmark4,
                                                          Ref.ech = cond1,
                                                          Comp.ech = cond4,
                                                          ntop = 5)
  
  DimPlot(cond4, group.by = "Bio.Order")
  
  DimPlot(cond2, group.by = "StrongForms")
  
  
  list_mark_1 <- ClusterHeatmap(sample = cond1, strong_markers = stmarks1)
  list_mark_2 <- ClusterHeatmap(sample = cond2, strong_markers = stmarks2)
  list_mark_3 <- ClusterHeatmap(sample = cond3, strong_markers = stmarks3)
  list_mark_4 <- ClusterHeatmap(sample = cond4, strong_markers = stmarks4)
  
  idk <- ClusterHeatmap(sample = cond1, strong_markers = stmarks1, ordered_list = list_mark_2)
  idk <- ClusterHeatmap(sample = cond1, strong_markers = stmarks1, ordered_list = list_mark_3)
  idk <- ClusterHeatmap(sample = cond1, strong_markers = stmarks1, ordered_list = list_mark_4)
  
  idk <- ClusterHeatmap(sample = cond2, strong_markers = stmarks1, ordered_list = list_mark_1)
  idk <- ClusterHeatmap(sample = cond2, strong_markers = stmarks1, ordered_list = list_mark_3)
  idk <- ClusterHeatmap(sample = cond2, strong_markers = stmarks1, ordered_list = list_mark_4)
  
  idk <- ClusterHeatmap(sample = cond3, strong_markers = stmarks1, ordered_list = list_mark_1)
  idk <- ClusterHeatmap(sample = cond3, strong_markers = stmarks1, ordered_list = list_mark_2)
  idk <- ClusterHeatmap(sample = cond3, strong_markers = stmarks1, ordered_list = list_mark_4)
  
  idk <- ClusterHeatmap(sample = cond4, strong_markers = stmarks1, ordered_list = list_mark_1)
  idk <- ClusterHeatmap(sample = cond4, strong_markers = stmarks1, ordered_list = list_mark_2)
  idk <- ClusterHeatmap(sample = cond4, strong_markers = stmarks1, ordered_list = list_mark_3)
  
  #--------------- comparer un peu les genes
  
  top_marks1 <- marks1 %>% group_by(cluster) %>% top_n(n = 5, wt = v4_log2FC)
  top_marks2 <- marks2 %>% group_by(cluster) %>% top_n(n = 5, wt = v4_log2FC)
  top_marks3 <- marks3 %>% group_by(cluster) %>% top_n(n = 5, wt = v4_log2FC)
  top_marks4 <- marks4 %>% group_by(cluster) %>% top_n(n = 5, wt = v4_log2FC)
  
  
  for (i in 1:6){
    tmp <- list("cond1"= top_marks1[which(top_marks1$cluster==i),]$gene,
                "cond2"= top_marks2[which(top_marks2$cluster==i),]$gene,
                "cond3"= top_marks3[which(top_marks3$cluster==i),]$gene,
                "cond4"= top_marks4[which(top_marks4$cluster==i),]$gene)
    
    p1 <- ggvenn(tmp,
                 fill_color = c("Blue","Green","Orange","Red"),
                 show_elements = TRUE,
                 label_sep = "\n",
                 set_name_size = 5,
                 set_name_color = c("Blue","Green","Orange","Red"),
                 text_size = 3,
                 text_color = "black",
                 fill_alpha=.5) +
      theme_void() +
      ggtitle(paste0("Intersection of top 5 wilcoxon marker genes in cluster ", i))
    
    print(p1)
  }
  
  
  top_stmarks1 <- stmarks1 %>% group_by(cluster) %>% top_n(n = 15, wt = v4_log2FC)
  top_stmarks2 <- stmarks2 %>% group_by(cluster) %>% top_n(n = 15, wt = v4_log2FC)
  top_stmarks3 <- stmarks3 %>% group_by(cluster) %>% top_n(n = 15, wt = v4_log2FC)
  top_stmarks4 <- stmarks4 %>% group_by(cluster) %>% top_n(n = 15, wt = v4_log2FC)
  
  for (i in 1:6){
    tmp <- list("cond1"= top_stmarks1[which(top_stmarks1$cluster==i),]$gene,
                "cond2"= top_stmarks2[which(top_stmarks2$cluster==i),]$gene,
                "cond3"= top_stmarks3[which(top_stmarks3$cluster==i),]$gene,
                "cond4"= top_stmarks4[which(top_stmarks4$cluster==i),]$gene)
    
    p1 <- ggvenn(tmp,
                 fill_color = c("Blue","Green","Orange","Red"),
                 show_elements = TRUE,
                 label_sep = "\n",
                 set_name_size = 5,
                 set_name_color = c("Blue","Green","Orange","Red"),
                 text_size = 3,
                 text_color = "black",
                 fill_alpha=.5) +
      theme_void() +
      ggtitle(paste0("Intersection of top 5 strong marker genes in cluster ", i))
    
    print(p1)
  }
  
  
}



