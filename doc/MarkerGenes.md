# Marker Genes

### BioOrder
A function that takes a sample and and the order of said classes on the UMAP projection, reorders them according to said order and uses this order as the Ident in the seurat object.


#### Arguments
- `sample` : A seurat object
- `bio.order`: order in which the classes should be rearranged. This order should be given uch as the index is the label of the class in the Strong Forms, and the number at the index position is the new label to give to this class, e.g its position from left to right on the umap. 

#### Returns

Seurat object with additional "Bio.Order" column in metadata, which is also used as idents.

### FindStrongMarkers
Function that computes marker genes for 


