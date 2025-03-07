################################################################################
#                      Clustering gene expression profiles                     #
################################################################################

library(Seurat)
library(coseq)
library(dplyr)


#Function that gets mean counts for each gene in each cluster -----------------
GetMeanCounts <- function(sample){
  
  cellsidents <- CellsByIdentities(sample)
  mean_counts <- NULL
  
  for (i in 1:length(cellsidents)) {
    tmp <- data.frame(LayerData(sample, assay = "SCT", layer = "counts", cells = cellsidents[[i]]))
    tmp$mean <- apply(tmp, 1, mean)
    mean_counts <- cbind(mean_counts, tmp$mean)
  }
  
  mean_counts <- data.frame(mean_counts)
  rownames(mean_counts) <- rownames(tmp)
  
  return(mean_counts)
}


#Function that runs coseq and filters the "trash" cluster
FilterGenes <- function(mean_counts){
  
  #Running coseq on all counts
  res1 <- coseq(mean_counts, K= 2:30, model = "kmeans", normFactors = 'none', seed = 7654)
  cl1 <- clusters(res1)
  print(plot(res1, graphs = c("profiles", "boxplots")))
  cl1 <- as.factor(cl1)
  
  #fetching the yprofiles
  profiles1 <- as.data.frame(profiles(res1))
  mean_profiles1 <- data.frame(matrix(ncol=ncol(profiles1),nrow=0, dimnames=list(NULL, colnames(profiles1))))
 
   #getting mean profile for each cluster
  for (i in levels(cl1)){
    tmp_profile <- apply(profiles1[which(cl1==i),], 2, mean)
    mean_profiles1 <- rbind(mean_profiles1, tmp_profile)
  }
  colnames(mean_profiles1) <- colnames(profiles1)
  
  #Computing variation between clusters
  n <- ncol(mean_profiles1)
  for (i in 1:(n-1)){
    mean_profiles1[, (i+n)] <- abs(mean_profiles1[,(i+1)]-mean_profiles1[,(i)])
  }
  mean_profiles1$maxvar <- apply(mean_profiles1[, 7:11], 1, max)
  
  #If maximum variation in the cluster is very small, this is a cluster to remove
  new_vec <- c()
  for (i in levels(cl1)){
    if(mean_profiles1[i,]$maxvar < 0.01) {
      new_vec <- c(new_vec, "poubelle")
      print(paste0("Removed cluster ", i))
    }
    else{
      new_vec <- c(new_vec, i)
    }
  }
  
  todelete = which(new_vec == "poubelle")
  genelist = c()
  for (i in 1:length(todelete)){
    genelist <- c(genelist, which(cl1 == todelete[i]))
  }
  
  #Removing genes belonging to clusters to remove
  clean <- mean_counts[-genelist, ]
  
  return(clean)
}

#Function that takes a clean gene set and labels them -------------------------


LabelGenes <- function (clean, proba_threshold = 0.8, seed = 7294 ){
  
  #Launching coseq a second time on the cleaned up gene set
  res <- coseq(clean, K= 2:30, model = "kmeans", normFactors = 'none', seed = seed)
  cl <- clusters(res)
  print(plot(res, graphs = c("profiles", "boxplots")))
  
  #Filtering out genes with a low posterior probability
  proba_post <- as.data.frame(coseq::kmeansProbaPost(as.matrix(clusters(res)), as.matrix(tcounts(res))))
  low_prob_genes <- c()
  for (gene in rownames(clean)) {
    cluster_number <- cl[gene]
    probability <- proba_post[gene, cluster_number]

    if (probability < proba_threshold) {
      low_prob_genes <- c(low_prob_genes, gene)
    }
  }
  
  filtered <- clean[!rownames(clean) %in% low_prob_genes, ]
  profiles <- as.data.frame(profiles(res))
  profiles <- profiles[!rownames(profiles) %in% low_prob_genes, ]
  cl <- cl[!names(cl) %in% low_prob_genes ]
  
  #Computing mean profile for each gene
  mean_profiles <- data.frame(matrix(ncol=ncol(profiles),nrow=0, dimnames=list(NULL, colnames(profiles))))
  
  cl <- as.factor(cl)
  
  for (i in levels(cl)){
    tmp_profile <- apply(profiles[which(cl==i),], 2, mean)
    mean_profiles <- rbind(mean_profiles, tmp_profile)
  }
  colnames(mean_profiles) <- colnames(profiles)
  
  n <- ncol(mean_profiles)
  for (i in 1:(n-1)){
    
    mean_profiles[, (i+n)] <- abs(mean_profiles[,(i+1)]-mean_profiles[,(i)])
  }
  
  #Finding the maximum difference between clusters
  mean_profiles$maxvar <- apply(mean_profiles[, 7:11], 1, max)
  mean_profiles$sumvar <- apply(mean_profiles[, 7:11], 1, sum)
  #The column in which the max difference is
  mean_profiles$columnmax <- apply(mean_profiles[, 7:11], 1, which.max)
  #Finding the column of max expression
  mean_profiles$max_expr <- apply(mean_profiles[, 1:6], 1, which.max)
  
  #Attributing the labels
  mean_profiles <- mean_profiles %>%
    mutate(
      label = case_when(
        sumvar< 0.04   ~ "flat",
        max_expr == 1 & maxvar > 0.4  ~ "strong_portal",
        max_expr == 1 & maxvar < 0.4  ~ "portal",
        max_expr == 2 & (max_expr+columnmax)/2 <3 ~ "periportal",
        max_expr == 2 & (max_expr+columnmax)/2 >=3 ~ "mid",
        max_expr == 3 | max_expr == 4 ~ "mid",
        max_expr == 5 ~ "pericentral",
        max_expr == 6 & maxvar > 0.4 ~ "strong_central",
        max_expr == 6 ~ "central",
        TRUE ~ "flat" # if none of the conditions are met
      )
    )
  
  #Creating results df
  mean_profiles$cluster <- rownames(mean_profiles)
  clustering_df <- data.frame(gene = names(cl), cluster = as.vector(cl))
  filtered$gene <- rownames(filtered)
  
  result_df <- filtered %>% 
    left_join(clustering_df, by = "gene") %>%
    left_join(mean_profiles[,16:17], by = "cluster")
  
  result_df$label <- as.factor(result_df$label)
  
  result_list <- list(result_df = result_df,
                      coseq_result = res)
  
  return(result_list)
  
}


#testing zone ------------------------------------------------------------------

testing = F

if testing { 
  samp <- readRDS("data/clean/samp_1340_bio_ordered.Rds")
  meancounts <- GetMeanCounts(samp)
  filtered <- FilterGenes(meancounts)
  labels <- LabelGenes(filtered)
  df_ctrl <- labels[["result_df"]]
  table(df_ctrl$label)


  samp <- readRDS("data/clean/samp_1352_bio_ordered.Rds")
  meancounts <- GetMeanCounts(samp)
  filtered <- FilterGenes(meancounts)
  labels <- LabelGenes(filtered, seed = 1412)
  df_stim <- labels[["result_df"]]
  table(df_stim$label)


  samp <- readRDS("data/clean/samp_1353_bio_ordered.Rds")

  meancounts <- GetMeanCounts(samp)
  filtered <- FilterGenes(meancounts)
  labels <- LabelGenes(filtered)
  df_dako <- labels[["result_df"]]
  table(df_dako$label)

  samp <- readRDS("data/clean/samp_1354_bio_ordered.Rds")
  meancounts <- GetMeanCounts(samp)
  filtered <- FilterGenes(meancounts)
  labels <- LabelGenes(filtered)
  df_dakostim <- labels[["result_df"]]
  table(df_dakostim$label)


  write.csv(df_ctrl, "labels_ctrl.csv")
  write.csv(df_stim, "labels_stim.csv")
  write.csv(df_dako, "labels_dako.csv")
  write.csv(df_dakostim, "labels_dakostim.csv")
  test <- read.csv("labels_ctrl.csv")

}

