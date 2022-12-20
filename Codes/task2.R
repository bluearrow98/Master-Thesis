rm(list=ls())

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Import required packages
library(parallel)
library(doParallel)
library(bestsubset)
source("helpFunctions.R")

# Parallel settings
numCores = makeCluster(detectCores())
registerDoParallel(numCores)
clusterExport(numCores,c("bs"))
clusterEvalQ(numCores,source("helpFunctions.R"))
clusterEvalQ(numCores,source("discreteFirstOrderAlgos.R"))
clusterEvalQ(numCores,library(glmnet))
clusterEvalQ(numCores,library(gurobi))


# Signal size
p = 5

# Number of signals 
N = 100


# Method for best subset
method = "bs"

set.seed(99)


for (s in p) {
  
  # Number of data to be sampled
  n = c(100)
  
  
  for (d in n) {
    
    elapsedTime <- system.time({bestResult = lapply(seq(1,N),recoverDriftMatrix,
                                                    p = s, n = d,
                                                    method = method)})
    
    if (method != "tLasso"){
      accMax = max(unlist(lapply(bestResult,function(x){x$acc})))
      tprAvg = mean(unlist(lapply(bestResult,function(x){x$tpr})))
      fprAvg = mean(unlist(lapply(bestResult,function(x){x$fpr})))
      f1scoreAvg = mean(unlist(lapply(bestResult,function(x){x$f1score})))
      aucrocAvg = mean(unlist(lapply(bestResult,function(x){x$aucroc})))
      aucprAvg = mean(unlist(lapply(bestResult,function(x){x$aucpr})))
      timeAvg = elapsedTime/N
    }
    
    save(accMax,tprAvg,fprAvg,f1scoreAvg,timeAvg,aucrocAvg,aucprAvg,
         file = paste0(getwd(),"/../Results/",d,"_",s,".RData"))
  
  }
}

stopCluster(numCores)