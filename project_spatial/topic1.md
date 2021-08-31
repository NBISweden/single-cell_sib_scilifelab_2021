# <img border="0" src="../logos/spatial_transcriptomics.png" width="40" height="40" style="vertical-align:middle;"> Spatial Transcriptomics - Topic1
***

<br/>



###  Milestone 0: Spatial tutorials

Start by having a look at the
[Tutorials](https://nbisweden.github.io/workshop-scRNAseq/exercises) for
analysis of ST data and integration with SC data using Seurat, Scran or Squidpy.

Select wichever pipeline you feel is most relevant for you to use.


###  Milestone 1: Load data

1.1. Load the relevant ST data for your project into the data object of your choice (AnnData, SeuratObject, SCE etc.).

1.2. Have a look at the object and describe in your report briefly what the different slots of the objects are.

1.3  Concatenate the different sections into one object if relevant for that data type.

<br/>

###  Milestone 2: Quality control

2.1 Calculate quality metrics such as number of genes, UMIs, and mitochondrial reads per spot.

2.2 Visualize the QC stats on a per spot basis, do you see even distribution across spots/sections?

2.3 Possibly filter out low quality spots, discuss with your group and motivate how the data should be filtered.

2.4 Filter genes, we suggest to always remove mitochondrial genes, but possibly also other problematic genes (e.g., ribosomal and Hb-genes).

<br/>


###  Milestone 3: Dimensionality reduction and clustering

3.1 Select a relevant set of genes for dimensionality reduction, motivate why you think this is a relevant choice.

3.2 Start with PCA and make a reasonable selection of number of principal components.

3.3 Run a graph based dimensionality reduction, like UMAP or another method of choice.

3.4 Cluster the spots with your method of choice, motivate why you are using that method.

<br/>


###  Milestone 4: Integration

4.1 Look at the dimensionality reduction and clustering in step 3. Do you see any batch effects per section? The sections are from 2 different mice, 3 and 2 consecutive sections per mouse, do you see any batch effects per animal? What would be the best way to remove that effect do you think?

4.2 Select integration features, how is this done best do you think? Discuss with your group and test a few different methods.

4.3 Select a method for integration. A good idea is that everyone in the group
tests different tools and then you compare the results, here it's also fit to
discuss what _metrics_ that should be used to gauge "performance" of each
method, i.e., how do you quantitatively measure how well the data has been integrated?

4.4 Rerun clustering and dimensionality reduction in the integrated space.

4.5 Evaluate the integration, is there still some batch effects visible or are
you happy with the results? It may be important to go back and tweak gene
selection, number of components etc in the integration util you are happy with
the results.

<br/>

###  Milestone 5: Cluster sections individually

5.1 Instead of working with all sections combined, instead run the analysis one section at a time.

5.2 How well does the clustering on section level agree with the clustering after integration? If they do not agree, how can you explain the differences?

<br/>

###  Milestone 6: Differential expression and spatially variable genes

6.1 Run a few different methods for differential expression among you in the group. How well do the methods agree? 

6.2 Run a few different methods for spatially variable genes and compare the results.

6.3 What is the main difference in the gene sets you get when you select variable genes, spatially variable genes or run DEG detection? Compare across the 3 methods and discuss with your group why the results are different.

<br/>

###  Milestone 7: Image features

7.1 Use a tool that can calculate image features

7.2 Use the image features to cluster the sections. How does it compare to the expression based clusters?

7.3 Do you see more/less batch effects in clustering based on image features? What could the reasons be?

<br/>


### Bonus tasks:

If time allows we have listed a few additional things you could try out if time allows.

* Nuclei segmentation. Try out a method for segmenting nuclei per spot.  Can you use the information on number of cells per spot to enhance the clustering? How well does predicted number of nuclei correlate with nUMI/nFeatures?



