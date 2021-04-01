---
title: <img border="0" src="https://static.thenounproject.com/png/67360-200.png" width="40" height="40"> Project report (single cell)
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 2
    toc_float:
        collapsed: false
    code_folding: show
  pdf_document: default
---

***

<br/>

# Research Analysis Tasks

The following milestones and tasks are what you will need to perform in order
to complete the project-based learning outcomes for this course. Go through
them in order and refer to the glossary for more detailed information on each
individual step.

## Milestone 1

### Download the data
The data for this project is available at the course GitHub repository, which
can be downloaded like so:

1. Run `git clone https://github.com/sib-swiss/SchoolRNA2020.git`
2. A new directory will be created named `SchoolRNA2020`; the data is residing
   in the `single_cell/data/` directory

You can either read the data from there (check the glossary for details) or
copy the data to some other location, if you prefer. You can run the above
command either on the command line interface or from within RStudio (more
details in the glossary).


### Load and merge datasets
- Consult the Glossary or additional sources for help
- Which file format do we have the data in?
- Why do we need to create a Seurat object?
- Where in the Seurat object is your counts stored?
- Describe in form of text the rational for this step in your markdown report.


### Compute QC
- Consult the Glossary or additional sources for help
- Which QC metrics should you calculate?
- What does each QC metric mean?
- Describe in form of text the rationale for this step in your markdown report


### Define appropriate filtering thresholds
- Consult the Glossary or additional sources for help
- Which QC metric needs filtering in your data?
- Which filtering parameters and respective thresholds suit your data?
- Describe in form of text the rationale for this step in your markdown report


<br/>


## Milestone 2

### Normalization and scaling
- Consult the Glossary or additional sources for help
- Why do we need to normalize the data? What exactly are we normalizing?
- Where in the Seurat object is your normalized data stored?
- Which data covariates could potentially influence the interpretation of the results?
- Following the question above, are there any covariates that need to be regressed out?
- Are all genes equally important for your analysis? Justify.
- Where in the Seurat object is your scaled data stored?
- Describe in form of text the rational for this step in your markdown report


### Data Visualization
- Consult the Glossary or additional sources for help
- Which method would you use to capture most significant information out of the data?
- Following the question above, which parameters would you choose? Why?
- Which method would you choose for visualization of the differences between your cells?
- Following the question above, how the parameters in this method influence your visual representation?
- How some of your QC parameters and datasets influence the separation of your cells?
- Where in the Seurat object is your reductions stored?
- Describe in form of text the rational for this step in your markdown report

<br/>


## Milestone 3

### Dataset integration
- Consult the Glossary or additional sources for help
- Are there any batch effects in the data? Is batch correction / dataset integration necessary?
- Following the question above, which parameters allow you to say that batch effects are present?
- Would a simple linear regression be sufficient for removing batch effects on your dataset?
- Which method for batch effect would you choose?
- How could you tell the batch correction procedure worked?
- Where in the Seurat object is your integrated data stored?
- After batch correction, do you have corrected data matrix (genes x samples) or a matrix in a reduction embedding (samples x dimensions).
- Where where should you put the embedding results?
- Visualize your results using the new matricies.
- Describe in form of text the rational for this step in your markdown report

<br/>


## Milestone 4

### Cell Clustering
- Consult the Glossary or additional sources for help
- What is a graph?
- Which kind of graph is the most robust to represent your data?
- Where in the Seurat object is your graph stored?
- Why do we need to cluster our cells?
- Which parameters have you chosen when clustering?
- How can you tell which clustering resolution is best?
- Do the clustering reflect the cell separation seen by the visualization method you are using?
- How are your clusters distributed across the samples, groups, experimental conditions, etc.?
- Where in the Seurat object is your clustering data stored?
- Describe in form of text the rationale for this step in your markdown report

<br/>


## Milestone 5

### Differential expression
- Consult the Glossary or additional sources for help
- Which biological question(s) do you want to answer with differential expression?
- Are you interested in comparing all cells or using a specific cluster?
- If you are interested in a particular cluster, which cluster, why?
- Which clustering resolution would you run your differential expression?
- Which test did you choose for differential expression?
- What parameters did you set for computing differential expression? Justify each one
- Which marker genes can separate each of the cell clusters in your data?
- Which cell types do they represent?
- How would you visualize the list of differentially expressed genes?
- Describe in form of text the rationale for this step in your markdown report

<br/>

### Trajectory inference analysis
- Consult the Glossary or additional sources for help
- Which biological question(s) do you want to answer with trajectory?
- Are you sure you have a developmental path in your data? 
- Are you interested in using all cells or using a specific cluster?
- Which embeddings will you use for computing trajectories? why?
- Which differential expression test are you interested in?
- How would you visualize your results?
- Describe in form of text the rationale for this step in your markdown report


