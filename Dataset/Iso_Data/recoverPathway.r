rm(list = ls())

# Import required packages
library(parallel)
library(doParallel)
library(bestsubset)
library(igraph)


mainCodePath <- paste0(getwd(), "/../../Codes/")
source(paste0(mainCodePath, "helpFunctions.R"))
source(paste0(mainCodePath, "modelSelection.R"))
source(paste0(mainCodePath, "metrics.R"))
source(paste0(mainCodePath, "discreteFirstOrderAlgos.R"))
source("pathwayRecoveryUtils.R")

# Parallel settings
numCores = makeCluster(detectCores())
registerDoParallel(numCores)
clusterExport(numCores, c("bs", "mainCodePath"))
clusterEvalQ(numCores, source(paste0(mainCodePath, "helpFunctions.R")))
clusterEvalQ(numCores, source(paste0(mainCodePath, "modelSelection.R")))
clusterEvalQ(numCores, source(paste0(mainCodePath, "metrics.R")))
clusterEvalQ(numCores, source(paste0(mainCodePath, "discreteFirstOrderAlgos.R")))
clusterEvalQ(numCores, source("pathwayRecoveryUtils.R"))
clusterEvalQ(numCores, library(glmnet))
clusterEvalQ(numCores, library(gurobi))


path <- 3
data <- divideIsoData(path)
signal <- list(list(list(setPathwayProb(data))))
n <- dim(data)[1]
p <- dim(data)[2]


split_signal <- train_test_split(signal)

method <- "lasso"
methodDir <- ""

for (s in 1:length(p)) {

  for (d in 1:length(n)) {

        elapsedTime <- system.time(results <-
            recoverDriftMatrix(split_signal$train[[s]][[d]], p, method, paras = lapply(1: length(split_signal$train[[s]][[d]]), function(x){
                                  floor(seq(p, p ** 2, length.out = 10))})))

        bestResult <- saveResults(split_signal$train[[s]][[d]], results,
            p, n, elapsedTime, methodDir)

        paras <- lapply(bestResult, function(x){return(x$bestPara)})
        
        testTime <- system.time(tesResults <- recoverDriftMatrix(split_signal$test[[s]][[d]], p, method, paras, test = TRUE))
  }
}