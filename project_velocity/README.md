# <img border="0" src="../logos/rna_velocity.png" width="40" height="40" style="vertical-align:middle;"> RNA velocity
***

<br/>

- __*Volker Bergen*__, Institute of Computational Biology, Munich, 🇩🇪 Germany
- __*Paulo Czarnewski*__, National Bioinformatics Infrastructure (NBIS), 🇸🇪 Sweden


### Background

The combination of velocities across genes is used to show the direction of movement of an individual cell in a dimensionality-reduced embedding. Incorrect visual directions can potentially arise from biases in the way velocities are projected. Simultaneously, the interpretation of the projected velocities is hampered by the difficulty in identifying individual gene dynamics that give rise to the projections. For instance, projections can be distorted due to multiple dynamic processes that occur simultaneously in a specific regime, such as cell cycle and differentiation. Here, we will analyse the impact of embedding choices on the vector field representation in the lower-dimensional space (Topic 1),
and systematically identify genes that giving rise to these projection, in order to address one of the major challenges to interpreting RNA velocity results (Topic 2).

<br/>

**RNA Velocity Topic 1:**

The choice of your embedding matters - its relevance and impact on interpretation of RNA velocity estimates. Final reports are available here:

<img border="0" src="https://www.svgrepo.com/show/165459/business-presentation.svg" width="30" height="30" style="vertical-align:middle;"> [SLIDES with analysis results](RNA_Velocity_group1.pdf)  

<img border="0" src="https://static.thenounproject.com/png/67360-200.png" width="30" height="30" style="vertical-align:middle;"> [HTML Project Report with step-by-step analysis](RNA_velocity_group1.html)  

<img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/38/Jupyter_logo.svg/1200px-Jupyter_logo.svg.png" width="30" height="30" style="vertical-align:middle;"> [Jupyter Notebook with step-by-step analysis](RNA_velocity_group1.ipynb)  

<br/>

**RNA Velocity Topic 2:**

Identifying genes that give rise to vector field representations, for better interpretation of RNA velocity results. Final reports are available here:

<img border="0" src="https://www.svgrepo.com/show/165459/business-presentation.svg" width="30" height="30" style="vertical-align:middle;"> [SLIDES with analysis results](RNA_Velocity_group2.pdf)  

<img border="0" src="https://static.thenounproject.com/png/67360-200.png" width="30" height="30" style="vertical-align:middle;"> [HTML Project Report with step-by-step analysis](RNA_velocity_group2.html)  

<img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/38/Jupyter_logo.svg/1200px-Jupyter_logo.svg.png" width="30" height="30" style="vertical-align:middle;"> [Jupyter Notebook with step-by-step analysis](RNA_velocity_group2.ipynb)  

<br/>

### Preparation

- Follow the instructions here [https://scvelo.org/installation](https://scvelo.org/installation) to install `scanpy` and `scvelo`.
- Run the tutorials at [https://scvelo.org](https://scvelo.org) and make yourself comfortable with the main steps.
- Single cell anlaysis pipeline: [https://nbisweden.github.io/workshop-scRNAseq/exercises](https://nbisweden.github.io/workshop-scRNAseq/exercises)

<br/>

### Recommended literature

Luecken and Theis (2019), Current best practices in single-cell RNA-seq analysis: a tutorial, [MSB](https://doi.org/10.15252/msb.20188746).

La Manno *et al.* (2018), RNA velocity of single cells, [Nature](https://doi.org/10.1038/s41586-018-0414-6).

Bergen *et al.* (2020), Generalizing RNA velocity to transient cell states through dynamical modeling, [Nature Biotech](https://doi.org/10.1038/s41587-020-0591-3).

Bergen *et al.* (2021), RNA velocity - current challenges and future perspectives, [MSB](https://doi.org/10.15252/msb.202110282).

Coenen and Pearce (2021), Understanding UMAP, [GitHub](https://pair-code.github.io/understanding-umap/).

<br/>

### Milestone 1: Read into AnnData object

1.1. Read the pancreas dataset

1.2. Why do we use an AnnData object?

1.3. Where are the un/spliced counts and annotations for observations and variables stored?  

1.4. Describe in form of text the rational for these steps in your markdown report.

<br/>

### Milestone 2: Preprocessing

2.1. Why do we need to filter out a number of genes?

2.2. What are meaningful gene selection criteria and why?

2.3. What is the rationale of normalizing the counts?

2.4. After gene selection and normalization the data is log-transformed. What is it necessary for?

2.5. Show how the PCA embeddings compare when you use log vs. non-log data.

<br/>

### Milestone 3: Neighbor graph, and moments

3.1. The PCA representation (computed on the log-normalized data) is used to compute a neighbor graph from cell-to-cell distances. What does this graph represent and what can it be used for?

3.2. Now, we use the graph to compute first-order moments, i.e., for each cell, the mean (first-order moments) and variance (second-order moments) of its neighboring cells is computed. This can also be regarded as a knn-imputation of the counts.

3.3. Show how the k parameter (number of neighbors) impacts the imputed counts, e.g., by comparing the imputed gene expression of two genes  with different values for k.

<br/>

### Milestone 4: Velocity estimation

4.1. For velocity estimation, per default the stochastic model is applied. It estimates the steady-state ratio for each gene, and computes velocities as the residuals of the observations from this steady-state line.

4.2. Demonstrate for a few genes, how the steady-state line is fit. Can you show mathematically why velocities are given as the residuals from this line?

4.3. Now, apply the dynamical model instead and show how the steady-state lines differ.

4.4. Now use a few gene phase portraits to demonstrate how inferred up- and down-regulation can be read from these portraits.

4.5. What are transient states as opposed to steady states? How many steady states exist?

<br/>

### Milestone 5: Embeddings

5.1. Compute various dimensionality reduction methods / embeddings (e.g., tsne, umap, diffmap) and show how they differ and what type of topology they capture (local vs. global similarity).

5.2. As you vary the main parameters (e.g., min_dist in umap), how does that impact the embedding?

5.3. Which of these embeddings would you use for visualising trajectories and why?

<br/>

### Milestone 6: Velocity graph

6.1. The velocity graph summarises all cell-to-cell transitions inferred from velocity. What exactly does the ij-th entry in the graph represent?  

6.2. Can you run a number of cell transitions using the graph to simulate a trajectory starting from a particular cell in one of the embeddings?

6.3. How do the trajectories compare, when run from different starting cells (e.g., cell cycle vs. differentiation)?

<br/>

### Milestone 7: Velocity projection

7.1. Now the high-dimensional velocity vector field is visualized by projection into the two-dimensional embedding. Show how the velocity graph (transition matrix) is being used to project the data. For instance, you can pick one particular cell, and show its potential transitions, and how the actual arrow is derived from these.

7.2. Vary with the n_neighbors parameter and show how that impacts the projected vector field.

7.3. Why is the vector field getting smoother/averaged with a large n_neighbors parameter?

<br/>

### Milestone 8: Topic 1 - Vector field representations depend on the embedding

8.1.1. Now show how the projected vector field is represented in the embedding as you vary some of the parameters from before (such as min_dist for UMAP).

8.1.2. Is there a specific parameter set, where the arrows and the embedding coordinates do not work together?

8.1.3. Which embedding provides the highest resolution into the vector field.

8.1.4. Using the simulated cellular transitions, can you show which embedding specification represent the high-dimensional vector field best?

8.1.5. Can you even systematically quantify, how well the vector field is represented in the embedding?
Possible metrics could, e.g., be correlation of embedded velocities with pseudo time / latent time, 'compactness' of transitions, correlation of transitions in the embedding weighted by their probabilities.

<br/>

### Milestone 8: Topic 2 - Identification of putative driver genes

8.1.1. Identification of dynamically relevant genes is motivated in three ways:

8.1.2. We would like to understand, which genes are the ones that drive the vector field representation in the embeeding, i.e., which are most impactful on the projection.

8.1.3. We would like to understand, which gene are transiently expressed, i.e., dynamically activating in the differentiation process.

8.1.4. Finally, we want to know which genes are actual driver of the underlying biological processes.

8.1.5. Using only a small number of selected genes, show how the velocity projection differs from the projection using all genes. Can you quantify that difference?

8.1.6. There are different ways of computing ‘relevant’ genes in scvelo. Show, which of them best represents the entire vector field.

8.1.7. Which genes are supporting / contradicting the embedded trajectories? How does the projection differ, as you subset to these genes or exclude these genes?

8.1.8. Find more ways of identifying important genes, e.g., using PCA on the velocity vector field.

<br/>

### Bonus Tasks - Topic 1:

9.1.1. Can you find other ways to projecting the data (e.g., parametric UMAP)?

9.1.2. Metrics: correlation with pseudo time / latent time, compactness of transitions, weighted correlation of trasition probabilities in embedding.

9.1.3. How ‘far’ are the transitions reaching, how gene-specific / cell-specific, how does that correlate with velocity?

9.1.4. Where does the displacement ‘land’ (by its closest sampled neighbors or parametric projection), and how does that vary with delta t.

<br/>

### Bonus Tasks - Topic 2:

9.1.1. How ‘far’ are the transitions reaching, how gene-specific / cell-specific, how does that correlate with velocity?

9.1.2. Can you predict the embedded arrows from high-dim velocities (e.g., using random forest or logistic regression) in order obtain weights as a means of interpreting gene's weights as their impact on the projection?

9.1.3. Can you think of ways to disentangle or decompose different processes that simultaneously occur (e.g., cell cycle vs. differentiation)?

<br/>
