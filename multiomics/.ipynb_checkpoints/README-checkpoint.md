## Integrative analysis of single-cell multi-omic data

The development of single-cell assays for multiple molecular modalities other than gene expression provides a powerful tool for investigating multiple dimensions of cellular heterogeneity. Here we will be analyzing gene expression and chromatin accessibility profiles generated with scRNA-seq and scATAC-seq from the same tissue to identify genomic regions involved in cellular differentiation. We will use data integration techniques for experimental designs where both data modalities are simultaneously profiled from the same cells (vertical integration) and for designs where scRNA-seq and scATAC-seq are applied to separate groups of cells from the same tissue (diagonal integration). 

### Project description

**Background:** Corticogenesis is the dynamic process that results in the formation of the cerebral cortex, and is characterized by the generation of excitatory glutamatergic neurons from cortical progenitors, and the differentiation of astrocytes and oligodendrocytes. Dynamic changes in the activity of cis-regulatory DNA elements underlie the complex phenotypic transformations that occur during development. 

**Main research question:** which non-coding regions are dynamically associated with expression of corticogenesis factors?

- Group 1: address the research question with diagonal integration of the unmatched assay data
- Group 2: address the research question with vertical integration of the multiome data

**Data:** human fetal brain cortex data from [Trevino et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.29.424636v2.full) ([source](https://github.com/GreenleafLab/brainchromatin)). 
