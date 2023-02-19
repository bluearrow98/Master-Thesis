rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("helpFunctions.R")
library(ROCR)
library(PRROC)

set.seed(200)
n = 6
comp = runif(n)
gt = c(1,1,0,1,0,1)


x = list()
# comp_sort <- c(0,sort(comp))
# Manual threshold
c <- seq(0,1,by=0.001)
for (i in 1:length(c)) {
  comp_i <- (comp > c[i])
  x$tpr[i] <- tpr(comp_i,gt)
  x$fpr[i] <- fpr(comp_i,gt)
}

pred = prediction(comp,gt)
perf = performance(pred,"auc")

aucrocOwn <- AUCROC(x)
aucROC <- roc.curve(comp,gt)
