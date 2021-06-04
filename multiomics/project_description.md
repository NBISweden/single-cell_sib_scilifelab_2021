# Project description

__to be expanded in template for report__

#### Background 
Corticogenesis is the dynamic process that results in the formation of the cerebral cortex, and is characterized by the generation of excitatory glutamatergic neurons from cortical progenitors, and the differentiation of astrocytes and oligodendrocytes. Dynamic changes in the activity of cis-regulatory DNA elements underlie the complex phenotypic transformations that occur during development. 

#### Main research question 
Which non-coding regions are dynamically associated with expression of corticogenesis factors?

#### Data 
Human fetal brain cortex data from [Trevino et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.29.424636v2.full) ([source](https://github.com/GreenleafLab/brainchromatin)). 

- Group 1: address the research question with diagonal integration of the unmatched assay data
- Group 2: address the research question with vertical integration of the multiome data

#### Objectives

1. Cell-level integration: Infer a common embedding/pseudotime axis for both modalities 
    * How can you evaluate that the embedding is sensible?
    * Do you gain resolution in cell-type clustering using both modalities?
    
2. Feature-level integration: Using the common embedding, find non-coding regions that are associated to gene expression of genes involved in cortical neuron development
    * Does it help to aggregate signal across multiple peaks?
    * What is gained with peak-to-gene mapping at the cluster level VS pseudotime or other continuous embedding strategy?
    * How can you filter out spurious correlations/associations?
    * What strategies can you use to corroborate your peak-to-gene associations? (e.g. TF motif enrichment, comparison with a negative control, co-accessibility) (check out-of-the-box solutions for peak to gene mapping)
    
3. Bonus: what genes are likely to be affected by mutations associated with in Autism Spectrum Disorder? annotations [here](https://github.com/GreenleafLab/Brain_ASD/tree/master/filtered_mutations)
