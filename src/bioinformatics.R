#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# This is the bioinformatics section of the F1000
# bioconductor workflow paper (Part I of the three parts that need execution).

## ---- init ----
.cran_packages  <-  c("ggplot2", "gridExtra")
.bioc_packages  <- c("dada2", "phyloseq", "DECIPHER", "phangorn")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
set.seed(100)

# you have to change this (unless we put the data on dropbox and
# automatically download

## ---- files ----
miseq_path <- file.path("data", "MiSeq_SOP")
filt_path <- file.path("data", "filtered")

if(!file_test("-d", miseq_path)) {
  dir.create(miseq_path)
  download.file("http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar",
               destfile = file.path(miseq_path, "StabilityNoMetaG.tar"))
  system(paste0("tar -xvf ", file.path(miseq_path, "StabilityNoMetaG.tar"),
               " -C ", miseq_path, "/"))
}

fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

## ---- profile ----
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

## ---- filter ----
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=10, truncLen=c(245, 160),
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)
}

## ---- derep ----
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

## ---- rates ----
ddF <- dada(derepFs[1:40], err=NULL, selfConsist=TRUE)
ddR <- dada(derepRs[1:40], err=NULL, selfConsist=TRUE)

## ---- plot-rates ----
plotErrors(ddF)
plotErrors(ddR)

## ---- dada ----
dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=TRUE)
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=TRUE)

## ---- merge ----
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

## ---- seqtab ----
seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

## ---- chimeras ----
seqtab <- removeBimeraDenovo(seqtab.all)

## ---- tax ----
ref_fasta <- "data/rdp_train_set_14.fa.gz"
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

## ---- msa ----
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

## ---- tree ----
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)


## ---- samdat ----
mimarks_path <- "data/MIMARKS_Data_combined.csv"
samdf <- read.csv(mimarks_path, header=TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove duplicate entries for reverse reads
rownames(seqtab) <- gsub("124", "125", rownames(seqtab)) # Fixing an odd discrepancy
all(rownames(seqtab) %in% samdf$SampleID) # TRUE

rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment", "host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass", "diet", "family_relationship", "genotype", "SampleID")
samdf <- samdf[rownames(seqtab), keep.cols]

## ---- phyloseq ----
ps <- phyloseq(tax_table(taxtab), sample_data(samdf),
               otu_table(seqtab, taxa_are_rows = FALSE), phy_tree(fitGTR$tree))
