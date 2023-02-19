rm(list=ls())

library(ggplot2)
library(reshape2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
p = c(5,6,7,8,9,10)
n = 1000

# Best subset results
databs <- data.frame(p=p,
                   t(sapply(p,function(x){load(paste0(getwd(),"/../Results/SimResults-10 point grid/"
                                                      ,x,"_",n,"_res.RData"))
                     return (c(accMax,f1scoreAvg,tprAvg,fprAvg,aucprAvg,aucrocAvg))})))
colnames(databs) <- c("p","Max. Acc","F1 Score Avg.","TPR Avg.","FPR Avg.", "AUC PR Avg.",
                    "AUC ROC Avg.")
databs <- melt(databs,id.vars = "p", variable.name = "Metrics",value.name = "metrics")

# Lasso results
datalasso <- data.frame(p=p,
                        bs = t(sapply(p,function(x){load(paste0(getwd(),"/../Results/LassoResults/"
                                                                ,n,"_",x,".RData"))
                          return (c(accMax,f1scoreAvg,tprAvg,fprAvg,aucprAvg,aucrocAvg))})))
colnames(datalasso) <- c("p","Max. Acc","F1 Score Avg.","TPR Avg.","FPR Avg.", "AUC PR Avg.",
                         "AUC ROC Avg.")
datalasso <- melt(datalasso,id.vars = "p", variable.name = "Metrics",value.name = "metrics")

# Combine data
combineData <- bind_rows(databs,datalasso, .id = "Method")


dir = getwd()
pdf(file = paste0(dir,"/../Results/Plots for Variable signal size/10 point grid/n",n,"_res.pdf"))
p1 <- ggplot(combineData)+
  geom_point(aes(x = as.character(p),y = metrics,group = Method,color = Method),size = 2) + 
  geom_line(aes(x = as.character(p),y = metrics,group = Method, color = Method),size = 1) +
  scale_x_discrete(limits = as.character(p)) +
  scale_color_manual(values = c("black", "orange"),
                     labels = c("MIP","Lasso")) + 
  labs(title = paste0("Performance measures for n = ",n),
       x = "Signal Size") +
  facet_wrap(Metrics ~ .,ncol = 2) + theme_bw()
print(p1)
dev.off()