rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
p = c(5,6,7,8,9,10)
n = c(100,200,500,1000,5000,10000,100000,Inf)

for (i in p) {
  
  # Best subset results
  databs <- data.frame(n=n,
                     bs = t(sapply(n,function(x){load(paste0(getwd(),"/../Results/SimResults/"
                                                        ,i,"_",x,"_res.RData"))
                        return (c(accMax,f1scoreAvg,tprAvg,fprAvg,aucprAvg,aucrocAvg))})))
  colnames(databs) <- c("n","Max. Acc","F1 Score Avg.","TPR Avg.","FPR Avg.", "AUC PR Avg.",
                      "AUC ROC Avg.")
  databs <- melt(databs,id.vars = "n", variable.name = "Metrics",value.name = "metrics")
  
  # Lasso results
  datalasso <- data.frame(n=n,
                     bs = t(sapply(n,function(x){load(paste0(getwd(),"/../Results/LassoResults/"
                                                             ,x,"_",i,".RData"))
                       return (c(accMax,f1scoreAvg,tprAvg,fprAvg,aucprAvg,aucrocAvg))})))
  colnames(datalasso) <- c("n","Max. Acc","F1 Score Avg.","TPR Avg.","FPR Avg.", "AUC PR Avg.",
                      "AUC ROC Avg.")
  datalasso <- melt(datalasso,id.vars = "n", variable.name = "Metrics",value.name = "metrics")
  
  # Combine data
  combineData <- bind_rows(databs,datalasso, .id = "Method")
  
  p1 <- ggplot(combineData)+
    geom_point(aes(x = as.character(n),y = metrics,group = Method,color=Method),size = 2) + 
    geom_line(aes(x = as.character(n),y = metrics,group = Method, color =Method),size = 1) +
    scale_x_discrete(limits = as.character(n)) +
    scale_color_manual(values = c("black", "orange"),
                       labels = c("MIP","Lasso")) + 
    labs(title = paste0("Performance measures for p = ",i),
         x = "Sample size") +
    facet_wrap(Metrics ~ .,ncol = 2) + theme_bw()
  
  
  dir = getwd()
  pdf(file = paste0(dir,"/../Results/Plots for Constant signal size/p",i,"_res.pdf"))
  print(p1)
  dev.off()
}