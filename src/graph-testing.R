#! /usr/bin/env Rscript

# File description ----------------------------------------------------
# Code for performing network testing using phyloseqGraphTest. The
# processed data in data/rds/ps_processed.rds is generated in the
# script preprocessing.R

## ---- network-setup ----------------------------------------------
setup_example(c("igraph", "phyloseq", "phyloseqGraphTest", "ggnetwork", "intergraph"))

## ---- ggnetwork --------------------------------------------------
net <- make_network(ps, max.dist=0.35)
sampledata <- data.frame(sample_data(ps))
V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$litter <- sampledata[names(V(net)), "family_relationship"]

## ---- ggnetwork-plot ---------------------------------------------
ggplot(net, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter)) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(0.5,"line")) +
  guides(col = guide_legend(override.aes = list(size = .25)))

## ---- mst -----------------------------------------------------------
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "mst")
gt$pval

## ---- mst-plot ------------------------------------------------------
plotNet1=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet1, plotPerm1)

## ---- knn-1 ---------------------------------------------------------
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)

## ---- knn-1-plot ----------------------------------------------------
plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)


## ---- knn-2 ---------------------------------------------------------
gt <- graph_perm_test(ps, "family_relationship",
                      grouping = "host_subject_id",
                      distance = "bray", type = "knn", knn = 2)

## ---- knn-2-plot ----------------------------------------------------
#plot_test_network(gt) +
#  theme(legend.text = element_text(size = 8),
#        legend.title = element_text(size = 9))
#plot_permutations(gt)

## ---- threshold-720 -------------------------------------------------
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "bray", type = "threshold.nedges", nedges = 720,
                      keep.isolates = FALSE)


## ---- threshold-720-plot --------------------------------------------
plotNet3= plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm3=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet3, plotPerm3)

## ---- threshold-2000-------------------------------------------------
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "bray", type = "threshold.nedges", nedges = 2000,
                      keep.isolates = FALSE)

## ---- threshold-2000-plot -------------------------------------------
#plot_test_network(gt) +
#  theme(legend.text = element_text(size = 8),
#        legend.title = element_text(size = 9))
#plot_permutations(gt)
