# Descriptive statistics

### FindCorrelatedGenes
A function that looks for genes that are highly correlated with mouse Hemoglobin genes. It is useful in case of blood contamination in certain samples in order to remove genes that are highly correlated with Hemoglobins. The genes that are identified should be checked in a database such as GTEx Portal for _Whole Blood_ expression.

**Arguments :**
- `sample`: a Seurat Object
- `samp_num`: string, the name of the sample, usually a numeric identifier (e.g. "1353")
- `expr`: threshold for Hb expression (this value depends on the default assay for Seurat Object). used to filter the high-hb expressing spots
- `corr`: correlation threshold, default: 0.40. The threshold should be adapted in case known hb-correlated genes are not captured

**Returns:**
- `corr_genes`:a list.
Contains:--  `Correlations` : list of dataframes of gene correlations for each hemoglobin gene.
--  `Top_genes`: list of dataframes containing genes that are above the correlation threshold for each hemoglobin gene.


### DescriptiveStatsPlots 

 A function that plots various informative descriptive statistics. The plots are specified as options by the user. By default, all plots are printed by the function.

This function checks that the data is normalized using SCTransform by verifying the existence of the SCT assay within the object and, if not, performs SCT normalization on the sample. It also calculates the percentage of mitochondrial genes if a “`percent.mt`” column (default Seurat nomenclature) is not found in the metadata.

Seurat V4 SCT assay structure will be tolerated by the function.

**Arguments:**
-`sample`: A Seurat object
-`plot.detection.rate` : Boolean, if TRUE, will generate a histogram of the detection rate of genes  by bins of 100
-`plot.UMI`: Boolean, if TRUE, will generate a histogram of distribution of UMI among spots
-`plot.QC`: Boolean, if TRUE, will generate violinplots of the number of features, number of UMIs and percentage of mitochondrial genes in raw counts and normalized data
-`plot.scatter`: Boolean, if TRUE, will generate a scatterplot of the number of counts as a function of the number of features in raw data and SCT Normalized data.
-`mt_prefix`: str, prefix for mitochondrial genes

**Returns:**
Prints the plots specified as arguments.




### VariableGenes
 A function that computes Highly Variables Genes and/or Spatially Variable Genes in a sample. Adds the list of identified genes to sample metadata.

**Arguments:**
-`sample`: A Seurat object
-`hvg.selection` : Selection method for HVG selection. Can be “vst”, “mvp” or “disp”. “vst” by default.
-`nfeatures` : Number of features to select as top variable genes
-`top_nfeatures`: number of features to show on featureplot
- `search_svg`: boolean, if TRUE, will search for Spatially Variable Genes
- `svg.selection`: selection method for SVG selection. “moransi” by default, can also be “markvariogram”.

**Returns:**
Seurat object with HVG and SVG dataframes in metadata.


### DimReduc
A function that performs necessary dimention reductions. Runs PCA and then runs UMAP on the PCA embeddings.

**Arguments:**
- `sample`: a Seurat object
- `verbose`: boolean, wether or not to print in console the avancement of PCA and UMAP runs

**Returns:** 
Seurat object with pca and umap dimention reduction slots
Dimheatmaps of the first dimension for PCA and UMAP


### VariableGenesComparaison
A function that looks for an overlap between the computed HVG and SVG. Returns a venn diagram of the overlap.

**Arguments:**

-sample: a Seurat object for which HVG and SVG have been computed and stored in HVG and SVG slots of metadata using VariableGenes() function.
- `top_hvg:` number of HVG to take into account for comparaison, default to 20
- `top_svg:` number of SVG to take into account for comparaison, default to 20

**Returns:** 
Prints the Venn diagram of the HVG and SVG for a given sample


### GeneVisualization
Creates various plots for showing the expression of a gene in a sample.

**Arguments:**

- `sample`: a Seurat object
- `gene`: Gene to visualize
- `feature.plot:` boolean, if true, plots gene expression on  selected dimention reduction
- `reduction`: dimention reduction to plot gene expression on. can be “umap” or “pca”,
- `spatial.plot`: boolean, if True, plots gene expression on slide
- `alpha` :between 0-1, transparency of points on spatial plot

**Returns:**
Prints spatial and/or feature plot of selected gene


