## Processed datasets

See `_data_loaders/prep_data.Rmd`

## Multiome data matrices 

* `Multiome_ATAC_clean_SCE[.RDS|.h5ad]` cleaned annotated matrix objects for Multiome scATAC assay, all cell types 
* `Multiome_ATAC_ext_trajectory_SCE[.RDS|.h5ad]` cleaned annotated matrix objects for Multiome scATAC assay, subset to excitatory neurons (GluN) trajectory
* `Multiome_RNA_clean_SCE[.RDS|.h5ad]` cleaned annotated matrix objects for Multiome scRNA assay, all cell types 
* `Multiome_RNA_ext_trajectory_SCE[.RDS|.h5ad]` cleaned annotated matrix objects for Multiome scRNA assay, subset to excitatory neurons (GluN) trajectory

## Single omic  data matrices
* `scRNA_clean_SCE[.RDS|.h5ad]` cleaned annotated matrix objects for scRNA assay (subset to one batch), all cell types 
* `scATAC_clean_SCE[.RDS|.h5ad]` cleaned annotated matrix objects for scATAC assay (subset to one batch), all cell types 

## Annotations 
* `TrevinoEtAl_cluster_annotations.csv` collapsed cluster annotations for all datasets (already incorporated in cleaned files colData/obs)