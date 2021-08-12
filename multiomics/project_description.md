# Project description

### Background 
Corticogenesis is the dynamic process that results in the formation of the cerebral cortex, and is characterized by the generation of excitatory glutamatergic neurons from cortical progenitors, and the differentiation of astrocytes and oligodendrocytes. Dynamic changes in the activity of cis-regulatory DNA elements underlie the complex phenotypic transformations that occur during development. 

### Data

Human fetal brain cortex data from [Trevino et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.29.424636v2.full) ([source](https://github.com/GreenleafLab/brainchromatin)). 

- Group 1: address the research question with diagonal integration of the unmatched assay data
- Group 2: address the research question with vertical integration of the multiome data

### Main goal

Identify non-coding genomic regions where chromatin accessibility is associated with expression of genes involved in excitatory neuron development.

We can sub-divide the work in the following steps:

- Step 1: Define a pseudotime ordering of the differentiating excitatory neurons from the nIPCs
- Step 2: Define the genes and accessibility features that you want to test for associations
- Step 3: Test for associations between chromatin accessibility and gene expression

### Step 1: Define a pseudotime ordering of the differentiating excitatory neurons from the nIPCs

From a K-nearest neighbor graph or other latent embedding, order the cells in pseudotime. For the purpose of this tutorial (i.e. not focused on trajectory inference), we can use simple trajectory inference algorithms, such as [diffusion pseudotime]()

**Unmatched assays:** this step will require you to find a meaningful common embedding for the cells profiled with different assays. The most common approach to this problem is to first summarise accessibility over gene promoters/bodies, then use a horizontal integration method (anchoring on the genes) to co-embed the cells from the two assays.

Different pipelines/papers use different strategies to reduce accessibility signal to a gene x cell matrix. 

- Counting fragments over gene bodies and promoters (implemented in [Signac](https://satijalab.org/signac/reference/GeneActivity.html))
- [ArchR gene scoring](https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html): uses counts over gene bodies and promoter + weighting for distance of regulatory elements around the gene
- [Cicero gene activity score](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/##cicero-gene-activity-scores): based on measuring co-accessibility between peaks
- [Gene scores from cisTopic](https://www.embopress.org/doi/full/10.15252/msb.20209438): taking the average de-noised accessibility signal from peaks around each gene ([example implementation](https://github.com/emdann/scATAC_prep/blob/master/N2_add_cistopic.ipynb))

As for integration methods, we recommend trying 

- Seurat CCA ([paper](), [vignette]())
- LIGER ([paper]()), [vignette]())

Alternatively, you could try out diagonal integration methods that don't apply any transformation to the features to generate a common embedding. These are more experimental, and will require more parameter tuning.

- UINMF ([vignette](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/UINMF_vignette.html)) - extension of LIGER, requires at least partial overlap (for example by generating a gene-level matrix of accessibility and analyze it together with the peak-accessibility matrix)
- SCOT ([code](https://github.com/rsinghlab/SCOT))
- bindSC ([paper](https://www.biorxiv.org/content/10.1101/2020.12.11.422014v1.full.pdf))
- MultiMAP ([code](https://github.com/Teichlab/MultiMAP))

N.B. This task is _hard_, and high-throughput multi-modal data hasn't been around enough for a single or a few methods to establish itself as the state-of-the-art! If all else fails, remember that you can always take the simpler approach of aggregating cells into cell types that are
confidently identified in both assays, for example by inspection of expression/accessibility of marker genes.

**Matched assays:** in principle, you could just generate a latent embedding and pseudotime ordering on the scRNA modality, and take that as a representative embedding of your cells. However, it might be beneficial to incorporate information from both modalities in the embedding, using what we defined as "vertical integration" strategies. We recommend trying:

- Multi-omics factors analysis: ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)) ([code](https://github.com/bioFAM/MOFA2)) ([vignette](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html))
- Seurat V4 Weighted Nearest Neighbor analysis ([paper](https://www.cell.com/cell/fulltext/S0092-8674%2821%2900583-3))([preprint](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1))([code](https://github.com/satijalab/seurat))([vignette](https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html))

Other interesting approaches: [j-UMAP](https://github.com/canzarlab/JVis-learn), Schema ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02313-2)) ([code](https://schema-multimodal.readthedocs.io/en/latest/overview.html)).

N.B. it will be interesting to assess what you gain from the join embedding, in comparison to just using the RNA. It's worth noting that while the differences might be very subtle in the case of RNA-ATAC integration, it has been shown that there is a tangible gain in clustering resolution using these vertical integration strategies on other modalities such as RNA-protein (CITE-seq data).

**Expected output for this step:** an assignment of a continuous value of pseudotime for all the cells (stored in `adata.obs` or `colData(sce)`)

### Step 2: Define the genes and accessibility features that you want to test for associations

Now that you have a common embedding and hopefully a common pseudotime axis, we now need to select the features that we will use to associate gene expression to chromatin accessibility. We could be greedy and try to test all the possible combinations of genes and peaks, but this would likely be prohibitively computational expensive (with > 30K genes and > 400k peaks!), not to mention the multiple testing burden we would have. In addition, we know that long range interactions on the genome are possible, but they are not very common, so it doesn't make too much sense to test for associations between genes and chromatin regions that are extremely far apart e.g. on different chromosomes. 

To start with, select a subset of interesting genes that seem to have a dynamic behaviour in the differentiation trajectory. You could select these with standard strategies for variable gene selection (in [R](https://bioconductor.org/books/release/OSCA/feature-selection.html), in [python](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html)), or you could correlate the log-normalized gene expression to the value of pseudotime, to select genes that increase or decrease in expression during differentiation.

Then, you'll need to decide _what type of features_ and _which features_ to use at the chromatin level:

- what type of feature: You could use all the peak identified by peak calling, or aggregate counts in peaks that contain the same transcription factor binding motifs (e.g. using [ChromVAR](https://github.com/GreenleafLab/chromVAR)), or aggregate counts over peaks in larger genomic windows (see the [working with genomic ranges]() example code)
- which features: especially if you pick features defined by genomic positions, it's best to subset the possible gene-region pairs to regions within a certain range of the gene. A window of 100k base pairs is often used.

**Expected output for this step:** a list, `data.frame` or other object storing possible gene-region pairs

### Step 3: Test for associations between chromatin accessibility and gene expression

At this step it's worth thinking again about **aggregating profiles from multiple cells**, using clustering or the KNN graph on the common embedding. This will be necessary for pairing features in the unmatched datasets. However, it will be advantageous also for the matched cells, as it helps with the high sparsity of the scATAC profiles. We suggest creating pseudobulk profiles of neighbouring cells in the KNN graph. As well as denoising, this step reduces the dimensionality of our dataset (on the cell axis), speeding up computations. You could sample _n_ cells and create pseudobulk profiles from the K-nearest neighbors of those cells (selecting neighbourhoods where there is a sufficient number of cells from both assays in the case of unmatched data). More robust algorithms to sample neighbourhoods on the KNN graph are implemented in [miloR](https://github.com/MarioniLab/miloR) or [MetaCell](https://github.com/tanaylab/metacell). Alternatively, overclustering with the Louvain or Leiden algorithm can be used.

To testing for correlation between gene expression and chromatin accessibility, a common approach is to compute a gene-peak correlation coefficient (Spearman ⍴ or Pearson R2) by correlating normalized scATAC-seq peak peak counts with the corresponding gene’s expression across all ATAC-RNA paired cells. Then the observed coefficient for each gene-peak pair is compared to the coefficients obtained matching the gene to _n_ "background peaks" with matching total accessibility and GC content to the tested one, accounting for technical biases between cells that are caused by PCR amplification or variable Tn5 tagmentation conditions. This results in a null peak-gene correlation distribution that is independent of peak-gene proximity, that can be used for significance testing (for reference see [Ma et al. 2020](https://www.cell.com/cell/fulltext/S0092-8674(20)31253-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420312538%3Fshowall%3Dtrue), [Kartha et al. 2021](https://www.biorxiv.org/content/10.1101/2021.07.28.453784v1.full)).

The chromVAR package implements simple functions to identify the background peaks:
```
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38) ## your reference genome

atac_sce <- readRDS("./processed_data/Multiome_ATAC_ext_trajectory_SCE.RDS")
atac_chromvar_sce <- addGCBias(atac_sce, genome=BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(atac_chromvar_sce, niterations = 50) ## This can take a while
```

The ArchR package also provides a [`Peak2GeneLinkage`](https://www.archrproject.com/bookdown/peak2genelinkage-with-archr.html) function to perform a similar mapping (however this only works if you follow their RNA-ATAC integration pipeline from start to finish).

**Test your workflow on a few genes first!**

### Evaluating your results

Take some time to assess whether your predicted associations make sense. Some ideas:
- Use positive control: e.g. we expect that peaks in the proximity of Transcription Start Sites should be more frequently linked to gene expression.
- Can you define a set of "negative control pairs" where you _don't_ expect to see association?
- Are peaks associated to genes enriched in Transcription Factor motifs? Have a look at [muon functions]() and [motifmatchr]()
- get creative!


### Bonus question
Can you find genes that are likely to be affected by mutations associated with in Autism Spectrum Disorder? Annotations [here].

---

<!-- ### Objectives

A. Cell-level integration: Infer a common embedding/pseudotime axis for both modalities 
    - How can you evaluate that the embedding is sensible?
    - Do you gain resolution in cell-type clustering using both modalities?
B. Feature-level integration: Using the common embedding, find non-coding regions that are associated to gene expression of genes involved in cortical neuron development
    - Does it help to aggregate signal across multiple peaks? And across multiple cells?
    - What is gained with peak-to-gene mapping at the cluster level VS pseudotime or other continuous embedding strategy?
    - How can you filter out spurious correlations/associations?
    - What strategies can you use to corroborate your peak-to-gene associations? (e.g. TF motif enrichment, comparison with a negative control, co-accessibility) (check out-of-the-box solutions for peak to gene mapping)
C. Bonus: what genes are likely to be affected by mutations associated with in Autism Spectrum Disorder? annotations here -->
<!-- 
### Multi-omics integration methods

A. **Horizontal integration:** multiple matrices with a common cell axis and different features 
    1. Multi-Omics Factor Analysis (MOFA)
    2. Weighted Nearest Neighbors analysis
    3. [j-UMAP](https://github.com/canzarlab/JVis-learn)
    4. Schema ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02313-2)) ([code](https://schema-multimodal.readthedocs.io/en/latest/overview.html))
   
   Others: [totalVI](https://docs.scvi-tools.org/en/stable/user_guide/notebooks/totalVI.html), [Multigrate](https://icml-compbio.github.io/2021/papers/WCBICML2021_paper_44.pdf),

B. **Diagonal integration:** multiple matrices with no common axis (different cells, different features)
    1. MultiMAP
    2. bindSC ([paper](https://www.biorxiv.org/content/10.1101/2020.12.11.422014v1.full.pdf))
    3. UINMF ([vignette](http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/UINMF_vignette.html)) - extension of LIGER, requires at least partial overlap

C. **_Verticalized_ integration:** transform multiple matrices with no common axis to match features (e.g. generate a score of gene accessibility from scATAC fragments/peaks, match proteins in CITE-seq to their coding gene in RNA)
    1. Seurat CCA
    2. LIGER
    3. Conos
    
    Others: the wealth of methods for batch correction...
  -->
