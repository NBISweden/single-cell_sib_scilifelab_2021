# <img border="0" src="../logos/omics_integration.png" width="40" height="40" style="vertical-align:middle;"> Omics Integration
***

<br/>

- __*Emma Dann*__, Wellcome Sanger Institute, Cambridge, 🇬🇧 United Kingdom
- __*Charlotte Soneson*__, Friedrich Miescher Institute, SIB/Basel, 🇨🇭 Switzerland

<br/>

### Background

Corticogenesis is the dynamic process that results in the formation of the cerebral cortex, and is characterized by the generation of excitatory glutamatergic neurons from cortical progenitors, and the differentiation of astrocytes and oligodendrocytes. Dynamic changes in the activity of cis-regulatory DNA elements underlie the complex phenotypic transformations that occur during development.

Here we will be analyzing Human fetal brain cortex data from [Trevino et al. 2021](https://www.sciencedirect.com/science/article/abs/pii/S0092867421009429) ([source](https://github.com/GreenleafLab/brainchromatin)) ([OA preprint](https://www.biorxiv.org/content/10.1101/2020.12.29.424636v2.full)) to study the interplay between chromatin accessibility and gene expression in early corticogenesis.

The main goal will be to identify non-coding genomic regions where chromatin accessibility is associated with expression of genes involved in excitatory neuron development.

* Group 1 will address the question with diagonal integration of the unmatched assay data (scRNA-seq data and scATAC-seq data from different cells)
* Group 2 will address the question with vertical integration of the multiome data (scRNA-seq and scATAC-seq for the same cells)

### Multiomics topic 1 - diagonal integration
<img border="0" src="https://www.svgrepo.com/show/165459/business-presentation.svg" width="30" height="30" style="vertical-align:middle;"> [SLIDES with analysis results](Multiomics_group1.pdf)


### Multiomics topic 2 - vertical integration
<img border="0" src="https://www.svgrepo.com/show/165459/business-presentation.svg" width="30" height="30" style="vertical-align:middle;"> [SLIDES with analysis results](Multiomics_group2.pdf)

### Pre-course

#### Software and data set-up 

**With AWS:** 

- TBA

- clone the course code repository locally
```
git clone https://github.com/NBISweden/single-cell_sib_scilifelab_2021.git
```

- Open the Jupyter notebooks from the JupyterHub GUI

**Without AWS:** if you prefer you can set up your own working environment locally

- Download the preprocessed data from [GDrive](https://drive.google.com/drive/folders/1YjHfhxk2Z62pTEOTu27G-AgKqQawKEBT?usp=sharing) 

- clone the course code repository locally
```
git clone https://github.com/NBISweden/single-cell_sib_scilifelab_2021.git
```

- create a clone of our conda environment named `sc2021-multiomics` (if needed [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) first):
```
cd single-cell_sib_scilifelab_2021/project_omics
conda env create --file multiomics-environment.yml
```

- activate the newly created conda environment
```
conda activate sc2021-multiomics
```

You may have to repeat the activation after a new login or after deactivating the environment with `conda deactivate`.

- open the Jupyter notebooks from the terminal (make sure the conda enviroment is active which is indicated by `(sc2021-multiomics) ` at the beginning of the command line prompt):
```
jupyter notebook ./multiomics_unmatched.ipynb
```


#### Knowledge requirements

We will have a warm-up session to start exploring our datasets before the project starts. 

You can start familiarizing yourself with some of the tools we will be using, trying out the examples in vignettes:

* Python framework for multi-modal data analysis - muon ([vignette](https://muon-tutorials.readthedocs.io/en/latest/single-cell-rna-atac/index.html))
* Dimensionality reduction for co-embedding of matched data:
  * Multi-Omics Factor Analysis - ([vignette in R](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html))([implementation in muon](https://muon.readthedocs.io/en/latest/omics/multi.html#multi-omics-factor-analysis))
  * Weighted Nearest Neighbors - ([vignette in R](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1))([implementation in muon](https://muon.readthedocs.io/en/latest/omics/multi.html#weighted-nearest-neighbours))
* Dimensionality reduction for co-embedding of unmatched data:
  * Seurat CCA - ([vignette in R](https://satijalab.org/seurat/articles/atacseq_integration_vignette.html))
* Tools for scATAC analysis
  * Motif enrichment analysis - ChromVAR ([vignette](https://greenleaflab.github.io/chromVAR/articles/Introduction.html))
  * Working with BED like genomic locations - GenomicRanges ([vignette](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html))
  * Extracting Ensembl gene annotations ([vignette](https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#1_Introduction))

### Dataset background

Corticogenesis is the dynamic process that results in the formation of the cerebral cortex, and is characterized by the generation of excitatory glutamatergic neurons from cortical progenitors, and the differentiation of astrocytes and oligodendrocytes. Dynamic changes in the activity of cis-regulatory DNA elements underlie the complex phenotypic transformations that occur during development.

Here we will be analyzing Human fetal brain cortex data from [Trevino et al. 2021](https://www.sciencedirect.com/science/article/abs/pii/S0092867421009429) ([source](https://github.com/GreenleafLab/brainchromatin)) ([OA preprint](https://www.biorxiv.org/content/10.1101/2020.12.29.424636v2.full)) to study the interplay between chromatin accessibility and gene expression in early corticogenesis.

You will find the preprocessed datasets in the `/data/multiomics/` directory on the AWS server (or alternatively, [on Google Drive](https://drive.google.com/drive/folders/1YjHfhxk2Z62pTEOTu27G-AgKqQawKEBT?usp=sharing))

* `gr1_unmatched_diagonal` contains data from unmatched scRNA-seq (19373 cells x 33197 genes) and scATAC-seq (6423 cells x 657930 peaks) assays.  
* `gr2_matched_vertical` contains data from matched scRNA-seq (8981 cells x 34104 genes) and scATAC-seq (8981 cells x 467315 peaks) assays.

You will find the same datasets saved both in anndata format for use in python (`*.h5ad`) and in SingleCellExperiment format for use in R (`*.RDS`). 

For each scATAC dataset we have also provide precomputed "gene activities", counting ATAC fragments over gene bodies and promoters, as implemented by the `Signac` function [`GeneActivity`](https://satijalab.org/signac/reference/geneactivity) (authors of this dataset did not share raw data because of patient privacy).

In the template notebooks we demonstrate how to preprocess and merge the single modality objects in MuData objects from the python package `muon`.

### Practical information

Your main goal will be to identify non-coding genomic regions where chromatin accessibility is associated with expression of genes involved in excitatory neuron development.

* Group 1 will address the question with diagonal integration of the unmatched assay data (scRNA-seq data and scATAC-seq data from different cells)
* Group 2 will address the question with vertical integration of the multiome data (scRNA-seq and scATAC-seq for the same cells)

In the project folder, you will find a template Jupyter Notebook guiding you through the steps for the integration project:

* Group 1: [`multiomics_unmatched.ipynb`](https://github.com/NBISweden/single-cell_sib_scilifelab_2021/blob/main/project_omics/multiomics_unmatched.ipynb)
* Group 2: [`multiomics_matched.ipynb`](https://github.com/NBISweden/single-cell_sib_scilifelab_2021/blob/main/project_omics/multiomics_matched.ipynb)

Because we will need to use both tools in R and in python, we provide an additional notebook illustrating how to use R code in jupyter environment using the [RPy2](https://rpy2.github.io/) framework - [`rpy2_interoperability_examples.ipynb`](https://github.com/NBISweden/single-cell_sib_scilifelab_2021/blob/main/project_omics/rpy2_interoperability_examples.ipynb)). 

<!-- 
### Milestone 1:

1.1.

1.2.

1.3.

1.4.

<br/>

### Milestone 2:

2.1.

2.2.

2.3.

2.4.

2.5.

<br/>

### Milestone 3:

3.1.

3.2.

3.3.

<br/>

### Milestone 4:

4.1.

4.2.

4.3.

4.4.

4.5.

<br/>

### Milestone 5:

5.1.

5.2.

5.3.

<br/>

### Milestone 6:

6.1.

6.2.

6.3.

<br/>

### Milestone 7:

7.1.

7.2.

7.3.

<br/>

### Milestone 8:

8.1.1.

8.1.2.

8.1.3.

8.1.4.

8.1.5.
 -->
