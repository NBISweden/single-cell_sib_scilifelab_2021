# <img border="0" src="../logos/spatial_transcriptomics.png" width="40" height="40" style="vertical-align:middle;"> Tips for Spatial Transcriptomics projects
***

<br/>

Here, we will try to list some issues, tips related to different programs etc.


## Loading data

### Loading h5ad in R

To convert a scanpy AnnData object to a Seurat object in R, you need to have SeuratDisk installed. It should then be easy to read it in in R, however, it is very sensitive to having the correct formats, naming conventions etc for it to work. So it may at times require some troubleshooting.

Loading the VISp scRNAseq dataset should work with:

```
library(Seurat)
library(SeuratDisk)
library(hdf5r)

Convert("VISp-sc_data.h5ad", dest = "h5seurat", overwrite = TRUE)
data <- LoadH5Seurat("VISp-sc_data.h5seurat")

```

