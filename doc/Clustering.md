# Clustering

### PrepSample
A function that prepares a Seurat Object for BayesSpace clustering: Computes FindNeighbors, converts the object to SingleCellExperiment class and runs SpatialPreprocess.

#### Arguments
- `sample`: a Seurat Object
- `skip.PCA`: boolean, wether to skip PCA in BayesSpace::SpatialPreprocess. Default to True
- `log.normalize`: Boolean, whether to log normalize the data. Default to false (data should be SCT normalized).
- `platform`: Spatial Transcriptomics platform, used to determine spot neighborhood structure.  Can be "Visium" (hexagonal) or "ST" (square). 

#### Returns
Sample converted to sce

### RunBayesSpace
A function that performs multiple runs of BayesSpace clustering

#### Arguments: 
- `sample`: A Seurat Object
- `sce`: sce object associated to the seurat object
- `seed`: random seed for reproductibility of the results. Default is 42
- `niter`: Number of runs of BayesSpace, default to 20
- `nrep`: Number of MCMC iterations for a single run of BayesSpace, default to 10000 (recommended number of iterations)
- `gamma`: Smoothing parameter. Default to 2, but values between 1-3 seem to work well.
- `parallel`: Boolean; wether to run BayesSpace in parallel. Set to True by default
- `nrun`: Number of prior calls to this function, default is 0 (initial call to RunBayesSpace). Used for saving the results in seurat object metadata.

#### Returns
Seurat Object with clusterings stored in metadata under "bayes.space.n", n being the nth run of BayesSpace

### CreateMatrix

A function that creates an ordered ARI matrix for one sample. Clusters the results of BayesSpace using hierarchical clustering

#### Arguments: 
- `sample`: A Seurat Object on which BayesSpace clustering has been performed through `RunBayesSpace`
- `method`: method to use for hierarchical clustering. Default is complete.

#### Returns
List containing:
- `tree`: The result of hierarchical clustering
- `ARI`: Ordered ARI matrix

### ClimbTree

Function that prunes a tree and returns consensus BayesSpace clustering if one is found, and 0 otherwise.

#### Arguments: 
`clust_list`: The list returned by `CreateMatrix`
`min_ari`: ARI threshold used to prune the tree. 

#### Returns
Names of clusters belonging to the consensus if they are found, 0 otherwise.

### FindClusters
Function that searches for a consensus clustering, and runs BayesSpace as long as none is found.

#### Arguments: 
`sample`: A Seurat Object
`sce`: SingleCellExperiment object associated with sample
`clustree`: Clustering list returned by `CreateMatrix`
`seed`: random seed for reproductibility of the results. Default to 1412
`niter`: number of re-iterations of BayesSpace 
`nrep`: number of MCMC iterations of one BayesSpace run
`ncluster`: number of clusters 
`gamma`:
`method`: 
`max_runs`: maximum of times BayesSpace should be re-run.
#### Returns
List containing :
`sample`: Seurat object
`ARI`: ARI matrix of the consensus clustering


### RematchClusters 
Function that rematches the label of clusters belonging to consensus

#### Arguments: 
`samp_list`: List containing sample and ARI matrix of consensus clustering (output of `FindClusters`)

#### Returns:
`sample`: Seurat Object with rematched labels for clusters belonging to consensus

### ConsensusCluster
Function that computes the StrongForms of consensus Clustering

#### Arguments:
`sample_list`: List containing sample and ARI matrix of consensus clustering (output of `FindClusters`)

#### Returns:
`sample`: Seurat Object containing `StrongForms`column in metadata

