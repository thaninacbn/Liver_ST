################################################################################
#                         Clustering - Source functions                        #
################################################################################

library(plyr)
library(dplyr)
library(Seurat)
library(BayesSpace)
library(ggplot2)
library(reshape2)
library(doParallel)
library(parallel)
library(foreach)
library(mclust)
library(microbenchmark)



#Prepare the sample for BayesSpace Clustering----------------------------------
# returns SCE !!
PrepSample <- function(sample, k = 20,
                        skip.PCA = T,
                        log.normalize = F,
                        platform = "Visium"){
  
  if (!exists("umap", where = sample@reductions) |
      !exists("pca", where = sample@reductions) ){
    stop("Dimention reductions were not computed.")
  }
  
  sample <- FindNeighbors(sample,
                          k.param= k,
                          reduction="pca",
                          dims = 1:30)
  
  diet.seurat = Seurat::DietSeurat(sample,
                                   graphs = "SCT_nn",
                                   dimreducs = c("pca")) #slim down Seurat obj prior to conversion
  sce = as.SingleCellExperiment(diet.seurat, assay = "SCT") #convert seurat to SCE
  colData(sce) = cbind(colData(sce),sample@images$slice1@coordinates) 
  
  sce = spatialPreprocess(sce, 
                          platform = platform,
                          skip.PCA = skip.PCA, 
                          log.normalize = log.normalize) 
  
  return(sce) # BEWARE: this returns the sce object !!!
}


#Run BayesSpace on prepped sample ----------------------------------------------
# returns samp
RunBayesSpace <- function(sample,
                          sce, 
                          seed = 42, 
                          niter = 20,
                          nrep = 10000,
                          ncluster = 7, 
                          gamma = 2, 
                          parallel = T,
                          nrun = 0) {
  
  if (!exists("umap", where = sample@reductions) |
      !exists("pca", where = sample@reductions) ){
    stop("Dimention reductions were not computed.")
  }
  
  # Create a random vector of numbers to be used as seed for BayesSpace random start
  set.seed(seed)
  random_vec <- round(runif(n=niter, min =1, max = 100), 0)
  
  if (parallel){
    #Initiallize parallel computation
    Ncpus <- parallel::detectCores()-1
    cl <- parallel::makeCluster(Ncpus)
    doParallel <- registerDoParallel(cl)
  

  # Create list of clusterings 
    list_clusters<- foreach(i = 1:length(random_vec),
                        .packages = c("BayesSpace")) %dopar% {
                          set.seed(random_vec[i])
                          sce = spatialCluster(sce,
                                               nrep = nrep,
                                               q= ncluster ,
                                               burn.in = 10,
                                               model = "t",
                                               gamma = gamma)
                          sample@meta.data[[paste0("bayes.space.", (i+nrun))]] <- as.factor(sce$spatial.cluster)
                        }
  # Stop parallel computation
    parallel::stopCluster(cl)
    
    # Add clusterings to metadata
    for (i in 1:niter){
      sample@meta.data[[paste0("bayes.space.", (i+nrun))]] <- list_clusters[[i]]
    }
  }
  
  else {
    
    for (i in 1:length(random_vec)){
      set.seed(random_vec[i])
      sce = spatialCluster(sce, 
                           nrep = nrep,
                           q=ncluster , 
                           burn.in = 10, 
                           model = "t", 
                           gamma = gamma)
      sample@meta.data[[paste0("bayes.space.", (i+nrun))]] <- as.factor(sce$spatial.cluster)
    }
  }

  
  return(sample)
}


#Create ARI matrix -------------------------------------------------------------

ARI_matrix<-function(cluster){
  #library(mclust)  
  if (class(cluster)=="list"){
    ari_matrix<-matrix(rep(0,length(cluster)^2),length(cluster),length(cluster))
    for (i in 1: length(cluster)){
      ari_matrix[i,i]=1
    }    
    for (i in 1:(length(cluster)-1)){
      for (j in (i+1):length(cluster)){
        ari_matrix[i,j]=adjustedRandIndex(cluster[[i]],cluster[[j]])
        ari_matrix[j,i]=ari_matrix[i,j]
      }
    }
  }else{
    ari_matrix<-matrix(rep(0,ncol(cluster)^2),ncol(cluster),ncol(cluster))
    for (i in 1: ncol(cluster)){
      ari_matrix[i,i]=1
    }    
    for (i in 1:(ncol(cluster)-1)){
      for (j in (i+1):ncol(cluster)){
        ari_matrix[i,j]=adjustedRandIndex(cluster[,1],cluster[,j])
        ari_matrix[j,i]=ari_matrix[i,j]
      }
    }    
  }
  
  return(ari_matrix)
}


#Make heatmap from ARI matrix --------------------------------------------------

heatmapari = function(ari){
  #library(reshape2)
  #library(scales)
  co=melt(ari)
  co$Var1=as.factor(co$Var1)
  co$Var2=as.factor(co$Var2)
  gari = ggplot(co, aes(Var1, Var2)) + # x and y axes => Var1 and Var2
    geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
    geom_text(aes(fill = co$value, label = round(co$value, 2)),size=2) + # write the values
    scale_fill_gradient2(low = "white", 
                         mid = "red", 
                         high = "blue", 
                         midpoint = 1) + # determine the colour
    theme(panel.grid.major.x=element_blank(), #no gridlines
          panel.grid.minor.x=element_blank(), 
          panel.grid.major.y=element_blank(), 
          panel.grid.minor.y=element_blank(),
          panel.background=element_rect(fill="white"), # background=white
          axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
          plot.title = element_text(size=20,face="bold"),
          axis.text.y = element_text(size = 12,face = "bold")) + 
    #ggtitle("ARI Plot") + 
    theme(legend.title=element_text(face="bold", size=14)) + 
    scale_x_discrete(name="") +
    scale_y_discrete(name="") +
    labs(fill="ARI value")
  
  #gari=heatmap.2(ari,trace="none",dendrogram="none",scale="none")
  return(gari)
}



#Create ordered ARI matrix for one sample --------------------------------------
# returns list with tree and ari matrix
CreateMatrix <- function(sample, method = "complete") {
  
  clusts <- sample@meta.data %>% select(starts_with("bayes.space."))
  nclust <- length(colnames(clusts))
  
  #Create a list of clusterings from the metadata 
  clusters <- list()
  for (i in 1:nclust){
    tmp <-paste0("bayes.space.", i)
    to_append <- sample@meta.data[[tmp]] 
    clusters <- append(clusters, list(to_append))
  }
  
  
  #TODO : mettre une sécurité si jamais clusters est vide
  
  #Create ari matrix from list of clustering
  options(digits=2)
  ari_matrix<-ARI_matrix(clusters)
  
  #Add row and column names to the matrix
  names <- colnames(clusts)
  rownames(ari_matrix)=names
  colnames(ari_matrix)=names
  
  #Security: turning the negative ARI values to 0
  ari_matrix[ari_matrix<0] <- 0
  
  #Hierarchical clustering of the clusterings
  distance_mat <- 1- ari_matrix
  clusters_clustering <- hclust(as.dist(distance_mat), method = method)
  print(plot(clusters_clustering))
  
  #Order the matrix according to hclust order
  cluster_matrix <- ari_matrix[clusters_clustering$order, clusters_clustering$order]
  print(heatmapari(cluster_matrix))
  
  clustering_info <- list(tree = clusters_clustering, 
                          ARI = cluster_matrix)
  
  return(clustering_info)
}



# Read through clustree and find an appropriate cluster of clusterings ---------

ClimbTree <- function(clust_list, min_ari = 0.85){
  #Set the constants from clustlist
  clustree <- clust_list$tree
  ari_mat <- clust_list$ARI
  n <- length(clustree$labels)
  tau = 1-min_ari # tau variable is 1-tau mathematically but whatever
  
  #Cutree according to given height (height == 1- ARI threshold)
  clusters <- cutree(clustree, h= tau)
  
  #Find most represented cluster 
  counts <- plyr::count(clusters)
  maj = which.max(counts$freq)
  
  #Check that this cluster contains enough members
  if (counts$freq[maj]>= (n/4)){
    #Check that the minimum ARI within this cluster is > to threshold
    consensus <- names(which(clusters==maj))
    min_arimat <- min(ari_mat[consensus,consensus])
    if(min_arimat> min_ari){
      print(consensus)
      #Return BayesSpace clusterings that are selected
      return(consensus)
    }
  }
  
  #Return 0 if conditions are not satisfied
  return(0)
  
}


# Find best group of clusterings------------------------------------------------

FindClusters <- function(sample,
                         sce, 
                         clustree,
                         seed = 1412, 
                         niter = 5,
                         nrep = 3000,
                         ncluster = 7, 
                         gamma = 2, 
                         method = "complete",
                         max.runs = 4) {
  
  
  
  nclusterings = length(clustree$tree$labels)
  majcluster <- ClimbTree(clustree)
  
  if (length(majcluster)==1){# If best clustering has not been found yet
     # bayesspace loop
    set.seed(seed)
    seed_vec <- round(runif(n=max.runs, min =1, max = 100), 0)
    
    for (i in 1:length(seed_vec)){
      sample <- RunBayesSpace(sample = sample,
                              sce = sce,
                              niter = niter,
                              nrep = nrep, 
                              ncluster = ncluster, 
                              nrun = nclusterings,
                              seed = seed_vec[i])
      nclusterings = nclusterings + niter
      new_ari <- CreateMatrix(sample = sample, nclust = nclusterings)
      majcluster <- ClimbTree(new_ari)
      
      if(length(majcluster) != 1)
        maj_ari <- new_ari[majcluster, majcluster]
        samp_list <- list(sample = sample, ari = maj_ari)
        return(samp_list)
    }
    
  }
  
  maj_ari <- clustree$ARI[majcluster, majcluster]
  samp_list <- list(sample = sample, ari = maj_ari)
  return(samp_list) }



# Rematch cluster labels -------------------------------------------------------

RematchClusters <- function(samp_list){
  
  sample <- samp_list$sample
  maj_cluster <- samp_list$ari
  
  #TODO sécurité au cas ou maj_cluster == 0 ou pas été calculé
  
  clusters <- rownames(maj_cluster)
  ref_cluster <- clusters[1]
  nclust <- length(levels(sample@meta.data[[ref_cluster]])) #number of clusters that were computed using bayesspace

  
  for (cluster in clusters) {
    comp_table <- table(sample@meta.data[[cluster]], sample@meta.data[[ref_cluster]])
    tmp_list <- c()
    
    for (i in 1:nclust){
      to_rename <- which.max(comp_table[i,])
      tmp_list <- c(tmp_list, to_rename)
    }
    
    
    tmp_cluster <- sample@meta.data[[cluster]]
    for (i in 1:length(tmp_list)) {
      tmp_cluster[which(sample@meta.data[[cluster]] == i)] = tmp_list[i]
    }
    
    sample@meta.data[[cluster]] <- tmp_cluster
    
  }
  
  return(sample)
  
}



#Create consensus clustering ---------------------------------------------------

ConsensusCluster <- function(sample_list){
  
  sample <- RematchClusters(sample_list)
  
  clusters_list <- list()
  for (clustering in rownames(sample_list$ari)){
    to_append <- sample@meta.data[[clustering]] 
    clusters_list <- append(clusters_list, list(to_append))
  } 
  
  final_clusters <- numeric(length(clusters_list[[1]]))
  
  # Iterate over each spot
  for (i in seq_along(final_clusters)) {
    
    # Extract the cluster assignments for the current spot from all vectors
    spot_assignments <- sapply(clusters_list, `[`, i)
    
    # Check if all cluster assignments are the same for this spot
    if (length(unique(spot_assignments)) == 1) {
      # If all vectors agree, assign the common cluster number
      final_clusters[i] <- unique(spot_assignments)
    } else {
      # If not, assign 0
      final_clusters[i] <- 0
    }
  }
  
  sample[["StrongForms"]] <- as.factor(final_clusters)
  
  
  return(sample)
  
}



#testing -----------------------------------------------------------------------

testing = F
benchmark = F

if (testing) {
  
  test <- readRDS("data/clean/samp_1340_DimReduc.rds")
  sce <- PrepSample(test)
  test <- RunBayesSpace(test, sce, niter = 10, nrep= 3000)
  clust <- CreateMatrix(test)
  
  
  majclust<- FindClusters(test, sce, clust)
  
  test2 <- ConsensusCluster(majclust)
  DimPlot(test2, group.by = "StrongForms")
  
}

if (benchmark) {
  library(microbenchmark)
  
  time <- microbenchmark(parallel = RunBayesSpace(test, sce, parallel = T, nrep = 1000, niter = 5),
                         normal = RunBayesSpace(test, sce, parallel = F, nrep = 1000, niter = 5),
                         times = 10)
  
}

