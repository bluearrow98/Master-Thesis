rm(list = ls())

# Import required packages
library(parallel)
library(doParallel)
library(bestsubset)
source("helpFunctions.R")
source("modelSelection.R")
source("metrics.R")
source("discreteFirstOrderAlgos.R")

# Parallel settings
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
p <- c(5)

# Sample sizes
n <- c(1000)

# Number of signals
N <- 100

# Method for best subset
method <- "bs"

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

for (s in 1:length(p)) {

  for (d in 1:length(n)) {

    elapsedTime <- system.time(results <-
      recoverDriftMatrix(signals[[s]][[d]], p[s], n[d], method))

    bestResult <- saveResults(signals[[s]][[d]], results,
      p[s], n[d], elapsedTime)

    getMetrics(signals[[s]][[d]], bestResult, p[s], n[d])

  }
}

stopCluster(numCores)