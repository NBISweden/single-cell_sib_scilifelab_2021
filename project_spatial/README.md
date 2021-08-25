# <img border="0" src="../logos/spatial_transcriptomics.png" width="40" height="40" style="vertical-align:middle;"> Spatial Transcriptomics
***

<br/>

- __*Alma Andersson*__, KTH Royal Institute of Technology, Stockholm, 🇸🇪 Sweden
- __*Åsa Björklund*__, National Bioinformatics Infrastructure (NBIS), 🇸🇪 Sweden

<br/>

### Background

Spatial transcriptomics (ST) allows for visualization and quantitative analysis of gene expression in tissue sections. There is a multitude of recent methods for ST, but here we will focus on data from the Visium platform.

Spatial transcriptomic data with the Visium platform is in many ways similar to single cell transcriptomics (SC). It contains UMI counts for about 5-20 cells per spot instead of single cells, but is still quite sparse in the same way as SC data is, but with the additional information about spatial location in the tissue.

ST is a fast evolving field with a fast development in analysis tools, in this course the aim is to get familiar with a few of these tools.

<br/>


### Useful links

* [Tutorials](https://nbisweden.github.io/workshop-scRNAseq/exercises) for analysis of ST data and integration with SC data using Seurat, Scran or Squidpy.
* [Museum of Spatial Transcriptomics](https://pachterlab.github.io/LP_2021/) Extensive overview and lists of tools for Spatial transcriptomics. 


<br/>

### Project topics

Each task will be done in a group of four students, but many of the steps involves testing different tools or different analysis parameters. So we suggest that you subdivide the tasks within your group as you see fit.

In both projects we will work with Visium mouse brain sections. 

#### Topic 1: Analysis of multiple visium datasets

In this project you will work with mouse brain sections from [Hasel et al.](https://www.nature.com/articles/s41593-021-00905-6) which can be downloaded [here](fill in link). 


#### Topic 2: Annotating celltypes in visium data using single cell data

In this project you will focus on the integration of single cell data with spatial data in the Visceral cortex  which can be downloaded [here](fill in link).  


### Practical information

We have provided conda recipies for some of the most common tools that you may want to use in the projects in this [folder](https://github.com/NBISweden/single-cell_sib_scilifelab_2021/tree/main/project_spatial/conda/)

If you are not familiar with conda, please have a look at the [Precourse material](../precourse). There are separate environments for R or python and we suggest you use whatever language you feel most comfortable with. 

Please refer to our [Tips and tricks](tips) where we try to list common issues with some of the most common tools.

<br/>

### Milestone 0: Spatial tutorials

Start by having a look at the [Tutorials](https://nbisweden.github.io/workshop-scRNAseq/exercises) for analysis of ST data and integration with SC data using Seurat, Scran or Squidpy.

Select wichever pipeline you feel is most relevant for you to use.


### Milestone 1: Load data

1.1. Load the relevant ST data for your project into the data object of your choice (AnnData, SeuratObject, SCE etc.).

1.2. Have a look at the object and describe in your report briefly what the different slots of the objects are.

1.3 (Task 1) Concatenate the different sections into one object if relevant for that data type.

<br/>

### Milestone 2: Quality control

2.1 Calculate quality metrics such as number of genes and mitochondrial reads per spot.

2.2 Visualize the QC stats per spot, do you see even distribution across spots/sections?

2.3 Possibly filter out low quality spots, discuss with your group and motivate how the data should be filtered. 

2.4 Filter genes, we suggest to always remove mitochondrial genes, but possibly also other problematic genes.

<br/>

### Milestone 3: Dimensionality reduction and clustering

3.1 Select a relevant gene set for dimensionality reduction, motivate why you think this is a relevant choice.

3.2 Start with PCA and make a reasonable selection of number of principal components.

3.3 Run a graph based dimensionality reduction, like UMAP or another method of choice.

3.4 Cluster the spots with your method of choice, motivate why you are using that method.

<br/>


### Milestone 4 (Task 1): Integration

4.1 Look at the dimensionality reduction and clustering in step 3. Do you see any batch effects? What would be the best way to remove that effect do you think?

4.2 Select integration features, how is this done best do you think? Discuss with your group and test a few different methods.

4.3 Select a method for integration. A good idea is that everyone in the group tests different tools and then you compare the results.

4.4 Rerun clustering and dimensionality reduction in the integrated space.

4.5 Evaluate the integration, is there still some batch effects visible or are you happy with the results? It may be important to go back and tweak gene selection, number of components etc in the integration util you are happy with the results.

<br/>

### Milestone 4 (Task 2): Subset the section

4.1 Identify the area of interest (VISp) region and create a new object with only the selected region.

4.2 Visualize the region

4.3 Rerun dimensionality reduction and clustering with the selected region.

<br/>

### Milestone 5 (Task 1): Cluster sections individually

5.1 Instead of working with all sections combined, instead run the analysis one section at a time.

5.2 How well does the clustering on section level agree with the clustering after integration? If they do not agree, how can you explain the differences?


### Milestone 5 (Task 2): Load SC data

5.1 Perform clustering and QC of the data and compare to the annotations provided by the paper.

5.2 Perform Differential Expression analysis, this may be a relevant gene set to use for aligning the SC data to the ST data.

<br/>

### Milestone 6 (Task 1): Image features

6.1 Use a tool that can calculate image features

6.2 Use the image features to cluster the sections. How does it compare to the expression based clusters?

6.3 Do you see more/less batch effects in clustering based on image features? What could the reasons be?

<br/>

### Milestone 6 (Task 2): Integrate or deconvolve

6.1 Look at existing methods for aligning SC and ST data and select a few that you think are relavant. Discuss with the group what to select and why. 

6.2 Select features for integration/deconvolution, this can be done using variable genes in SC data only, both datasets or DEGs from SC data. Discuss what you think is the best method to use and possibly try out a few different ones if time allows.

6.3 Run the methods you selected.

6.4 Compare the results, what are main similiarities/differeces between the tools?

<br/>

### Bonus tasks:

If time allows we have listed a few additional things you could try out if time allows.

* Image segmentation, can you use the information on number of cells per spot to enhance the deconvolution/clustering?
*
* 

<br/>

