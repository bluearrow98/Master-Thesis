rm(list=ls())

library(ggplot2)
library(reshape2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
p = 5
n = c(100,200,500,1000,5000,10000,100000,Inf)

data <- data.frame(n=n,
                   t(sapply(n,function(x){load(paste0(getwd(),"/",p,"_",x,"_res.RData"))
                      return (c(accMax,f1scoreAvg,tprAvg,fprAvg))})))
colnames(data) <- c("n","Max. Acc","F1 Score","TPR Avg.","FPR Avg.")
data <- melt(data,id.vars = "n", variable.name = "Metrics",value.name = "metrics")


dir = getwd()
method = "bs"
pdf(file = paste0(dir,"/../Results/Plots for Constant signal size/p",p,"_bs_res.pdf"))
p1 <- ggplot(data,aes(x = as.character(n),y = metrics,group = 1))+
  geom_point(size = 2) + geom_line(size = 1) +
  scale_x_discrete(limits = as.character(n)) +
  labs(title = paste0("Performance measures for p = ",p),
       x = "Sample size") +
  facet_wrap(Metrics ~ .,ncol = 2) + theme_bw()
print(p1)
dev.off()