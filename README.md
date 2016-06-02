

# Bioconductor Workflow for Microbiome Data Analysis: from raw reads to community analyses.

To run the code in the .Rnw documents in this repository (eg):

```S
library("knitr")
knit("PartIIphyloseq.Rnw")
```

The three Rnw files can be executed independently. 

PartIdada.rnw is the read/bioinformatics preprocessing and is fairly time-consuming, 
(it can take 3-6 hours on modern laptops).

PartIIphyloseq.rnw performs a few analyses with phyloseq and is quite fast.

PartIIIanalysis.rnw performs all the statistical analyses and the machine learning components can take about 10
minutes.


