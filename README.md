# Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses.

## Generate output documents

To generate the output pdfs from the .Rnw documents in this repository (eg):

```S
library("knitr")
knit("PartIIphyloseq.rnw")
```

The three .rnw files can be processed independently.

The bioinformatics processing is fairly time-consuming, and takes 3-6 hours on modern laptops. In addition, an internet connection is required for the bioinformatics section, as the fastq files being processed are not included in this repository but must be downloaded.

## Run interactively

If may be of more use to run this workflow interactively inside an R session. To do that, the .rnw files may be ignored, as all necessary code is contained in the .R script files.

Before running the commands in the .R files, your working directory must be set to the base directory of this repository (i.e. the directory containing the .rnw files):

```
setwd("/path/to/F1000_workflow/") # CHANGE ME
```

The code for the analysis portion of the workflow is broken up into a number of different component .R files: analysis-setup.R, preprocessing.R, ordinations.R, supervised.R, graph-testing.R, linear-modeling.R, hierarchical-test.R and multitable.R. After the code in the first two files (analysis-setup.R and preprocessing.R), the code in the remaining files can all be run independently of the others.
