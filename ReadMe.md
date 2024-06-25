
# LiverST

A pipeline for the pre-processing, clustering and differential expression analysis of *(liver)* Spatial Transcriptomics data. 

Student project carried out in collaboration between the I2MC and INSA Toulouse for the Master's Thesis titled "*Analysis of Liver Spatial Transcriptomics Data*" written for the MSc in Bioinformatics at the Université of Paris Cité. 

This pipeline aims to simplify and automate the steps for Spatial Transcriptomics analysis. It is written in R and built upon open-source packages such as Seurat, BayesSpace, DEseq2 and more!

## Usage


This pipeline is built using Seurat V5 (Hao et. al, 2023) and uses the Seurat V5 assay. Users using objects stored as older Seurat V3  (used in Seurat V4 and below) assays are encouraged to update their objects using this command :
```
# convert a v3 assay to a v5 assay
brain[["Spatial5"]] <- as(object = brain[["Spatial3"]], Class = "Assay5")
```

The entry data should be pre-processed (i.e. filtered for low quality spots and genes such as hemoglobins) and normalized using SCTransform. 

The pipeline should be used in this order :
 - Descriptive statistics (`Statdesc_source.R`)
 - Clustering (`clustering_source.R`)
 - Marker genes identification (`marker_genes_source.R`)
 - Profile reconstruction (`profile_reconstruction_source.R`)
 - Gene labelling (`gene_labelling_source.R`)
 - Differential Gene Expression across conditions (`Integration_DE_source.R`)

