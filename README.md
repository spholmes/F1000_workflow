

# Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses.

## Generate output documents

To generate the output pdfs from the .Rnw documents in this repository (eg):

```S
library("knitr")
knit("PartIIphyloseq.rnw")
```

The three Rnw files can be executed independently. 

`PartIdada.rnw` is the read/bioinformatics processing and is fairly time-consuming, (it can take 3-6 hours on modern laptops). In addition, an internet connection is required, as the fastq files being processed are not included in this repository but must be downloaded.


`PartIIphyloseq.rnw` performs a few analyses with phyloseq and is quite fast.

`PartIIIanalysis.rnw` performs all the statistical analyses and the machine learning components can take about 10
minutes.


## Run interactively

If may be of more use to run this workflow interactively inside an R session. To do that, the .rnw files may be ignored, as all necessary code is contained in the .R script files.

Before running the commands in the .R files, your working directory must be set to the base directory of this repository (i.e. the directory containing the .rnw files):

```
setwd("/path/to/F1000_workflow/") # CHANGE ME
```

The code for the analysis portion of the workflow is broken up into a number of different component .R files: analysis-setup.R, preprocessing.R, ordinations.R, supervised.R, graph-testing.R, linear-modeling.R, hierarchical-test.R and multitable.R. After the code in the first two files (analysis-setup.R and preprocessing.R), the code in the remaining files can all be run independently of the others.
