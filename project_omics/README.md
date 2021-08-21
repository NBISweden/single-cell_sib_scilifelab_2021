# <img border="0" src="../logos/omics_integration.png" width="40" height="40" style="vertical-align:middle;"> Omics Integration
***

<br/>

- __*Emma Dann*__, Wellcome Sanger Institute, Cambridge, 🇬🇧 United Kingdom
- __*Charlotte Soneson*__, Friedrich Miescher Institute, SIB/Basel, 🇨🇭 Switzerland

<br/>

### Pre-course

Clone the code repository locally
```
git clone https://github.com/NBISweden/single-cell_sib_scilifelab_2021.git
```

Create a new conda environment and install the required python packages (if needed [install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) first):
```
cd single-cell_sib_scilifelab_2021/project_omics
conda env create --file sc_sib_scilifelab_2021_multiomics.yml
```

Install some additional R packages. From R:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("ensembldbr", "AnnotationDbi", 'EnsDb.Hsapiens.v86',
    "SingleCellExperiment", "GenomicRanges", "scran",
    'BSgenome.Hsapiens.UCSC.hg38', "chromVAR", "MOFA2",
    "MultiAssayExperiment"))

install.packages(c('tidyverse', "Seurat", "Signac", "Matrix"))
```

Download the data
```
mkdir processed_data
## TBD ##
```

### Background

Corticogenesis is the dynamic process that results in the formation of the cerebral cortex, and is characterized by the generation of excitatory glutamatergic neurons from cortical progenitors, and the differentiation of astrocytes and oligodendrocytes. Dynamic changes in the activity of cis-regulatory DNA elements underlie the complex phenotypic transformations that occur during development.

Here we will be analyzing Human fetal brain cortex data from [Trevino et al. 2021](https://www.sciencedirect.com/science/article/abs/pii/S0092867421009429) ([source](https://github.com/GreenleafLab/brainchromatin)) ([OA preprint](https://www.biorxiv.org/content/10.1101/2020.12.29.424636v2.full)) to study the interplay between chromatin accessibility and gene expression in early corticogenesis.


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
