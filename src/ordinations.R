#! /usr/bin/env Rscript

# File description ---------------------------------------------------
# Code for performing ordinations on the processed data. You must
# specify base_dir; this must the path to the BiocWorkflow directory.

## ----ordinations-bray---------------------------------------------------------
setup_example(c("phyloseq", "ggplot2", "plyr", "dplyr", "reshape2",
                "ade4", "ggrepel"))
out.bc.log <- ordinate(pslog, method = "MDS", distance = "bray")

## ----ordinations-bray-plot----------------------------------------------------
evals <- out.bc.log$values$Eigenvalues
plot_ordination(pslog, out.bc.log, color = "age_binned") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age")

## ----ordinations-dpcoa---------------------------------------------------------
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")

## ----ordinations-dpcoa-plot----------------------------------------------------
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned",
                shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")

## ----ordinations-dpcoa-species ----
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))

## ----ordinations-wuf---------------------------------------------------------
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")

## ----ordinations-wuf-plot---------------------------------------------------
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")

## ---- ordinations-wuf-species ----
plot_ordination(pslog, out.wuf.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))

## ---- pca-rank-get-ranks ----
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))

## ---- pca-rank-truncate ----
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1

## ---- pca-rank-visualize-procedure ----
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = .7) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")

## ---- pca-rank-pca-setup ----
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)

main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))

row_scores <- row_scores %>%
  left_join(sample_data(pslog))
col_scores <- col_scores %>%
  left_join(tax)

## ---- pca-rank-pca-plot ----
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

## ---- ccpna-correspondence-analysis ----
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)

## ---- ccpna-join-data ----
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))

species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)

## ---- ccpna-plot-age ----
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                  aes(x = CCA1, y = CCA2, label = otu_id),
                  size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.33) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

## ---- ccpna-plot-litter ----
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                  aes(x = CCA1, y = CCA2, label = otu_id),
                  size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45  ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
