### TO DO / DISCUSS

- more prefiltering of the ATAC data (remove empty peaks)
- run dimensionality reduction methods on scATAC data (cisTopic, peakVI, LSI, ChromVAR)?
- Vignette and package requirements for MOFA
- Prepare vignettes and references for:
    - GenomicRanges and a python equivalent
    - R to python conversions: using R in python, using python in R, 
    - List of useful tools: chromVAR, Seurat vignettes, muon vignettes, MOFA
    
- Should the unmatched group have access to the matched dataset too? e.g. for benchmarking
- Should the data be in the github repo or not?
- Should peak to gene mapping function be provided?
- Compute and save ChromVAR background peaks? (just in case)
- Should peak to gene search downsized to a set of known marker genes?