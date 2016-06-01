#! /usr/bin/env Rscript

# File description -------------------------------------------------------------
# Functions for hierarchically testing hypotheses.

## ---- deseq-transform ----
setup_example(c("phyloseq", "structSSI", "plyr", "dplyr", "reshape2",
                "ggplot2", "DESeq2"))
ps_dds <- phyloseq_to_deseq2(ps,  ~ age_binned + family_relationship)
varianceStabilizingTransformation(ps_dds, blind = TRUE, fitType = "parametric")
ps_dds <- estimateSizeFactors(ps_dds)
ps_dds <- estimateDispersions(ps_dds)
abund <- getVarianceStabilizedData(ps_dds)

## ---- structssi-shorten-names ---- 
short_names <- substr(rownames(abund), 1, 5)%>%
  make.names(unique = TRUE)
rownames(abund) <- short_names

## ---- deseq-vis ----
abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(pslog)),
                               sample = rownames(otu_table(pslog)),
                               type = "log(1 + x)"))

ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")

## ---- structssi-unadjp ----
el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$age_binned)

## ---- structssi-test ----
hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)
plot(hfdr_res, height = 5000) # opens in a browser

## ---- structssi-tax ----
tax <- tax_table(pslog)[, c("Family", "Genus")] %>%
  data.frame()
tax$seq <- short_names

hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>%
  head(10)
