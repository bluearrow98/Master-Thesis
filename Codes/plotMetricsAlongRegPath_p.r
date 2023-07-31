rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)

setwd(getwd())
p <- c(5, 10, 15, 20, 25)
n <- c(100, 1000, 10000, 1e+05, Inf)

# Colour blind friendly colours
cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Collect results
method <- list("ResultsLassoBIC", "ResultsWithRegEdgeBIC", "ResultsWith1EdgeBIC", "ResultsWithGlassoBIC", "ResultsWithLassoIniBIC")
data <- vector('list', length(method))

for (j in n) {
  for (i in 1:length(method)) {
    data[[i]] <- data.frame(p = p,
      t(sapply(p, function(x){
        load(paste0(getwd(), "/../Results/", method[[i]], "/completeResults_edgeProb25/"
            , j, "_", x, "_res.RData"))

        return (c(results$max_acc, results$max_f1, results$aucroc, results$aucpr))})))
    colnames(data[[i]]) <- c("p", "Max. Acc", "Max. F1 Score", 
        "AUC ROC", "AUC PR")
    data[[i]] <- melt(data[[i]], id.vars = "p", variable.name = "Metrics",
      value.name = "metrics")
  }

  # Combine data
  combineData <- bind_rows(data, .id = "Method")


  dir = getwd()
  # pdf(file = paste0(dir,"/../Results/Plots for Variable signal size/10 point grid/", 
  #       "MetricsAlongRegPath/d25/n", j, "_res.pdf"))
  p1 <- ggplot(combineData)+
    geom_point(aes(x = as.character(p),y = metrics,group = Method,color = Method),size = 2) + 
    geom_line(aes(x = as.character(p),y = metrics,group = Method, color = Method),linewidth = 0.5) +
    scale_x_discrete(limits = as.character(p)) +
    scale_color_manual(values = cbbPalette[1:5],
                      labels = c("Lasso", "MIP_reg", "MIP_1edge", "MIP_glasso", "MIP_lasso")) +
    labs(title = paste0("Performance measures for n = ", j),
        x = "Signal Size") + theme(axis.text = element_text(size = 8)) +
    facet_wrap(Metrics ~ ., ncol = 2)
  print(p1)
  dev.off()
}