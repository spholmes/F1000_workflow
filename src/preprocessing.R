
## ---- preprocessing-setup ----
pkgs <- c("phyloseq", "ggplot2", "shiny", "miniUI")
invisible(sapply(pkgs,require, character = TRUE))
ps <- raw_data() 

## ---- preprocessing-plot ----
qplot(sample_data(ps)$age, geom = "histogram") + xlab("age")
qplot(log10(rowSums(otu_table(ps)))) +
  xlab("Logged counts-per-sample") 

## ---- outlier-detect-------------------------------------------------
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
sample_data(pslog)$age_binned <- cut(sample_data(pslog)$age,
                                     breaks = c(0, 100, 200, 400))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")

## ---- outlier-detect-plot ----
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))

## ----shiny-----------------------------------------------------------
select_ggplot = function(p) {
    points = c()
    p$data$id = 1:nrow(p$data)
    
    ui = miniPage(
        gadgetTitleBar("Brush to select points",
                       miniTitleBarButton("add", "Add points")),
        miniContentPanel(
            plotOutput("plot", height = "100%", brush = "brush")
        )
    )
    server = function(input, output, session) {        
        output$plot = renderPlot({
            p
        })
        observeEvent(input$add, {
            points <<- c(points, brushedPoints(p$data, input$brush)$id)
        })
        observeEvent(input$done, {
            points <<- c(points, brushedPoints(p$data, input$brush)$id)
            stopApp(subset(p$data, id %in% points))
        })
    }
    runGadget(ui, server)
}

## ----outlier-analyze------------------------------------------------------------
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram") +
  xlab("Relative abundance")

