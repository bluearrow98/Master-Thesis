rm(list=ls())

source("helpFunctions.R")
source("discreteFirstOrderAlgos.R")
library(bestsubset)
library(glmnet)
library(ROCR)
library(parallel)
library(doParallel)

# Parallel settings
numCores = makeCluster(detectCores())
registerDoParallel(numCores)
clusterExport(numCores,c("bs"),envir = environment())
clusterEvalQ(numCores,source("helpFunctions.R"))
clusterEvalQ(numCores,source("discreteFirstOrderAlgos.R"))
clusterEvalQ(numCores,library(glmnet))
clusterEvalQ(numCores,library(gurobi))


# Signal size
p = 5

# Number of signals 
N = 100

# Number of data to be sampled
n = 100000

set.seed(99)

elapsedTime <- system.time({bestResult = lapply(seq(1,N),recoverDriftMatrix,
                                                p = p, n = n,numCores = numCores,
                                                method = "bs")})


accMax = max(unlist(lapply(bestResult,function(x){x$acc})))
tprAvg = mean(unlist(lapply(bestResult,function(x){x$tpr})))
fprAvg = mean(unlist(lapply(bestResult,function(x){x$fpr})))
f1scoreAvg = mean(unlist(lapply(bestResult,function(x){x$f1score})))
timeAvg = elapsedTime/N

stopCluster(numCores)