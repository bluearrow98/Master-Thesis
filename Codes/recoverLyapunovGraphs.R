rm(list = ls())

# Import required packages
library(parallel)
library(doParallel)
library(bestsubset)
source("helpFunctions.R")
source("modelSelection.R")
source("metrics.R")
source("discreteFirstOrderAlgos.R")

# Parallelization settings
numCores = makeCluster(detectCores())
registerDoParallel(numCores)
clusterExport(numCores, c("bs"))
clusterEvalQ(numCores, source("helpFunctions.R"))
clusterEvalQ(numCores, source("modelSelection.R"))
clusterEvalQ(numCores, source("metrics.R"))
clusterEvalQ(numCores, source("discreteFirstOrderAlgos.R"))
clusterEvalQ(numCores, library(glmnet))
clusterEvalQ(numCores, library(gurobi))

# Signal size
p <- c(5, 10, 15, 20, 25)

# Sample sizes
n <- c(100, 200, 500, 1000, 5000, 10000, 100000, Inf)

# Number of signals
N <- 100

# Method (bs or lasso)
method <- "lasso"
methodDir <- "ResultsLasso" # Directory where the results need to be saved

set.seed(99)

# Get all signals
signals <- lapply(p, function(x) {
  lapply(seq(1, N), getSignals, x)
})

signals <- lapply(1 : length(signals), function(x) {
  return(lapply(n, function(y) {
    return(setupProblems(signals[[x]], p[x], y))
  }))
})

split_signals <- train_test_split(signals)

for (s in 1:length(p)) {

  for (d in 1:length(n)) {

    elapsedTime <- system.time(results <-
      recoverDriftMatrix(split_signals$train[[s]][[d]], p[s], method))

    # Evaluate Pre model selection metrics
    getMetricsAlongRegPath(split_signals$train[[s]][[d]], results, p[s], n[d], elapsedTime, methodDir)

    bestResult <- saveResults(split_signals$train[[s]][[d]], results,
      p[s], n[d], elapsedTime, methodDir)

    getMetrics(split_signals$test[[s]][[d]], bestResult, p[s], n[d], method, methodDir)

  }
}

stopCluster(numCores)