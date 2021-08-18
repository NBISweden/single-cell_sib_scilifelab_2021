### TO DO / DISCUSS

- Where should the data be stored?
- Should peak to gene search downsized to a set of known marker genes?
- How much example code / notebooks to provide? 
  - peak to gene mapping functions (see `peak2gene_example.Rmd`)
  - Peak accessibility to gene activity
  - KNN smoothing / aggregation for denoising in peak 2 gene mapping
- How much precomputed data to provide? 
  - dimensionality reduction on scATAC data (cisTopic, LSI, ChromVAR)
  - ChromVAR background peaks (just in case)
- Should the unmatched group have access to the matched dataset too? e.g. for benchmarking
- Prepare vignettes for:
    - GenomicRanges and a python equivalent
    - R to python conversions: using R in python, using python in R, 
    - MOFA usage (to test software requirements)
- Software requirements / environment 

---
- precompute latent dimensions for ATACX
- make clean folder for data sharing
- remove all zero genes/peaks after subsetting
