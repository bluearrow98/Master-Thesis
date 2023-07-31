library("Rlab")
library(MASS)
library(Matrix)
library(parallel)
library(doParallel)
library(gurobi)

# Generate stable matrices
stabSignal <- function(p){
  # edge probability
  k = 0.25 * p
  M = matrix(rnorm(p**2) * rbern(p ** 2, prob = k / p), nrow = p)
  nonDiagAbsSum = sum(abs(M)) - sum(diag(abs(M)))
  diag(M) = - nonDiagAbsSum - abs(rnorm(p))
  return(M)
}

# Obtain the population covariance matrix by solving the Lyapunov equation
getPopCovar <- function(M, C){

  p <- dim(M)[1]
  AM <- kronecker(diag(p), M) + kronecker(M, diag(p))
  E <- eigen(AM)
  Sigma <- solve(AM, -c(C))
  Sigma <- matrix(Sigma, nrow = p)

  return(Sigma)
}

getSignals <- function(i, p) {

  set.seed(i)

  # True signal and Volatility matrix
  M_star <- stabSignal(p)
  trueSparsity <- (M_star != 0)
  C <- diag(runif(p, min = 0.2))

  # True number of subsets
  k0 <- sum(trueSparsity)

  # Population covariance
  Sigma_star <- getPopCovar(M_star, C)

  return(list(M = M_star, C = C, popSigma = Sigma_star, gt = trueSparsity, k0 = k0))
}
 
 # Inverse covariance estimation
invCov <- function(Mvec,Cvec){
  # Get number of variables
  p <- sqrt(length(Mvec))

  # Get the drift matrix and volatility matrix
  M <- matrix(Mvec, nrow = p)
  C <- matrix(Cvec, nrow = p)

  # Get the recovered covariance matrix
  Sigma <- getPopCovar(M, C)
  
  # Find the inverse using decomposition
  E <- eigen(Sigma)
  Lambda <- diag(E$values^(-1))
  Gamma <- E$vectors
  invCov <- Gamma %*% Lambda %*% t(Gamma)

  return(invCov)
}

# Get observations
sampleX <- function(n,Sigma){
  X = mvrnorm(n, mu = rep(0, dim(Sigma)[1]), Sigma)
  return(X)
}

# Compute the sample covariance matrix
getDataCovar <- function(X){
  n <- dim(X)[1]
  Sigma <- t(X) %*% X / n
  
  return(Sigma)
}


# Choose the original BS or custom BS
applyBS <- function(i, j, signals) {

  ASigma <- signals[[i]]$ASigma
  vecC <- signals[[i]]$vecC
  sampleSize <- signals[[i]]$n

  n <- dim(ASigma)[1]
  p <- length(vecC)
  # Modified MIO formulation
  result <- customBS(ASigma, -vecC, j, sampleSize, standardize = TRUE)

  # Original MIO formulation by Hastie et al. (By default is commented. Uncomment to activate)
  # result <- bs(ASigma, -vecC, j, intercept = F,
  #             params = list(NonConvex = 2), form = 2)

  # Store objVal and constructed signal
  if (p >= n) {
    result$beta <- result$x[(n + 1):(n + p)]
    }
  else {
    result$beta <- result$x[1:p]
    }

  # Unscale the coefficients  
  result$beta <- result$beta / result$scale_fac

  return(result)

}

# Custom MIO formulation for drift matrix recovery
customBS <- function(X, y, k, sampleSize, polish = FALSE, standardize = FALSE) {

  n <- dim(X)[1]
  p <- dim(X)[2]

  if (standardize){
    # Get the covariance matrix before doing standardization (needed for glasso initialization)
    Sigma <- getSigma(X)
    scale_fac <- apply(X, 2, function(x) {sqrt(var(x) * (n - 1) / n)})
    X <- apply(X, 2, function(x) (x  - mean(x))/sqrt(var(x) * (n - 1) / n))
  }

  # Discrete First order method (uncomment based on requirement)
  # algo2time <- system.time(algo2Sol <- proj_grad(X, y, k, polish = polish))
  # algo2time <- system.time(algo2Sol <- proj_grad_oneEdge(X, y, k, polish = polish))
  # algo2time <- system.time(algo2Sol <- proj_grad_glasso(X, y, k, Sigma = Sigma, N = sampleSize, polish = polish))
  projTime <- system.time(projSol <- applyLasso(X, y, k))

  M_p <- projSol$bm

  ## Choosing bounds
  eps <- 1e-5
  d_ind <- seq(1, p, by = sqrt(p) + 1)
  Mu_d <- 2 * norm(as.matrix(M_p[d_ind]), type = "i")
  if (min(abs(M_p[d_ind])) < 1e-8) Mu_d_upper <- eps else Mu_d_upper <- min(abs(M_p[d_ind])) - eps
  Mu_nd <- 2 * norm(as.matrix(M_p[-d_ind]), type = "i")
  if (Mu_nd == 0) Mu_nd <- Mu_d

  ## Initial guess (Not used for cold start)
  z0 <- rep(0,p)
  z0[M_p != 0] <- 1

  diagA <- rep(-Mu_nd, p)
  diagA[d_ind] <- -Mu_d

  model <- list()
  params <- list()
  if (p < n) {
    # Warm-start model creation (Eq 2.5)
    model$A <- spMatrix(nrow = 1 + 2 * p, ncol = 2 * p, i = c(rep(1, p), rep(seq(2, 1 + p), 2), rep(seq(2 + p, 1 + 2 * p), 2)),
                       j = c(seq(p + 1, 2 * p), rep(seq(1, 2 * p), 2)),
                       x = c(rep(1, p), rep(1, p), diagA, rep(-1, p), diagA))
    model$Q <- spMatrix(nrow = 2 * p,ncol = 2 * p,
                        i = c(rep(seq(1, p), each = p)),
                        j = c(rep(seq(1, p), p)), x = c(0.5 * t(X) %*% X))
    model$obj <- c(-t(X) %*% y, rep(0, p))
    model$ub <- c(rep(-Mu_d_upper, p), rep(1, p))
    model$ub[-c(d_ind, seq(p + 1, 2 * p))] <- Mu_nd
    model$lb <- c(rep(-Mu_d, p), rep(0, p))
    model$lb[-c(d_ind, seq(p + 1, 2 * p))] <- -Mu_nd
    model$rhs <- c(k, rep(0, 2 * p))
    model$sense <- c("<=", rep("<=", 2 * p))
    model$vtype <- c(rep("C", p), rep("B", p))
    model$start <- c(M_p, z0)

    ## Creating SOS constraints (not used anymore. SOS constraints are included as linear inequalities)
    # sos <- list()
    # for (j in 1:p){
    #   sos[[j]] <- list()
    #   sos[[j]]$type <- 1
    #   sos[[j]]$index <- c(j, (p + j))
    #   sos[[j]]$weights <- c(1, 2)
    # }
    # model$sos <- sos

  }else{
    # Warm-start model creation (Eq 2.6)
    M_zeta <- max(colSums(apply(abs(X), 1, sort,
                    decreasing=TRUE)[, 1:k, drop = FALSE])) * Mu_d
    model$A <- spMatrix(nrow = n + 1 + 2 * p, ncol = n + 2 * p,
                       i = c(seq(1, n), rep(seq(1, n), each = p), rep(n + 1, p), rep(seq(n + 2, n + 1 + p), 2), rep(seq(n + 2 + p, n + 1 + 2 * p),2)),
                       j = c(seq(1, n), rep(seq(n + 1, n+p), n), seq(n + p + 1, n + 2 * p), rep(seq(n + 1, n + 2 * p), 2)),
                       x = c(rep(1, n), c(-t(X)), rep(1, p), rep(1, p), diagA, rep(-1, p), diagA))
    model$Q <- spMatrix(nrow = n + 2 * p, ncol = n + 2 * p,
                       i = seq(1, n), j = seq(1, n), x = rep(0.5, n))
    model$obj <- c(rep(0, n), -t(X) %*% y, rep(0, p))
    model$ub <- c(rep(M_zeta, n), rep(-Mu_d_upper, p), rep(1, p))
    model$ub[-c(seq(1, n), n + d_ind, seq(n + p + 1, n + 2 * p))] <-
      c(rep(Mu_nd, p - sqrt(p)))
    model$lb <- c(rep(-M_zeta, n), rep(-Mu_d, p), rep(0, p))
    model$lb[-c(seq(1, n), n + d_ind, seq(n + p + 1, n + 2 * p))] <-
      c(rep(-Mu_nd, p - sqrt(p)))
    model$rhs <- c(rep(0, n), k, rep(0, 2 * p))
    model$sense <- c(rep("=", n), "<=", rep("<=", 2 * p))
    model$vtype <- c(rep("C", n), rep("C", p), rep("B", p))
    model$start <- c(X %*% M_p, M_p, z0)

    ## Creating SOS constraints (not used anymore. SOS constraints are included as linear inequalities)
    # sos <- list()
    # for (j in 1:p){
    #   sos[[j]] <- list()
    #   sos[[j]]$type <- 1
    #   sos[[j]]$index <- c((n + j), (n + p + j))
    #   sos[[j]]$weights <- c(1, 2)
    # }
    # model$sos <- sos

  }

  ## Parameters for Gurobi optimization
  params$OutputFlag <- 1
  params$timeLimit <- 180
  params$NonConvex <- 2
  params$MIPGap <- 0.01

  result <- gurobi(model, params)
  
  result$projTime <- projTime
  result$projSol <- projSol
  result$scale_fac <- scale_fac # Only used when standardization is turned on. It is typically always turned on

  return(result)
}

recoverDriftMatrix <- function(signals, p, method,
                                paras = lapply(1: length(signals), function(x){
                                  floor(seq(p, p ** 2, length.out = 10))}),
                                  test = FALSE) {

  if (method == "bs") {
    k <- paras
    # Choose the original formulation or custom formulation (uncomment as required)
    # resultBS <- lapply(k, applyBS, ASigma, vecC)
    resultBS <- foreach(i = 1:length(signals)) %:%
                  foreach(j = k[[i]]) %dopar% {
                    applyBS(i, j, signals)
                  }
                  
    if (!test){

      # Check if all solutions obtained are optimal or sub-optimal with less than 5 % MIP gap
      feasibleResults <- lapply(1:length(resultBS), function(x){
        checkMIPGap(resultBS[[x]], k[[x]])
      })

      resultBS <- lapply(feasibleResults, function(z) {
            return(z$results)
      })
      feask <- lapply(feasibleResults, function(z) {
            return(z$k)
      })
    }else{
      feask <- k
    }
    # Get solutions and objective values
    betas <- lapply(resultBS, sapply, function(x) {
          return(x$beta)
          })

    objvals <- lapply(resultBS, sapply, function(x) {
          return(x$objval)
          })

    parameters <- feask


  }else if (method == "lasso") {
      # Lasso

      penalty <- rep(1, p ** 2)
      penalty[seq(1, p ** 2, by = p + 1)] <- 0

      if (!test){

        initGrid <- seq(10, 10^-5, length.out = 100)
        coarseLasso <- foreach(i = 1:length(signals)) %dopar% {
                          glmnet(signals[[i]]$ASigma, -signals[[i]]$vecC, 
                            intercept = FALSE, alpha = 1,
                            standardize = TRUE, penalty.factor = penalty,
                            lambda = initGrid)
                        }

        # Find min lambda such that M is diagonal
        lambdas <- lapply(coarseLasso, function(z) {
              return(min(initGrid[(colSums(penalty * (z$beta != 0)) == 0 &
            colSums((1 - penalty) * (z$beta != 0)) == p)]))
              })
        fineGrid <- lapply(lambdas, function(z) {
              return(seq(z, z / 10^4, length.out = 100))
              })
        fineLasso <- foreach(i = 1:length(signals)) %dopar% {
                          glmnet(signals[[i]]$ASigma, -signals[[i]]$vecC, 
                            intercept = FALSE, alpha = 1,
                            standardize = TRUE, penalty.factor = penalty,
                            lambda = fineGrid[[i]])
                        }
        # Get solution                    }
        betas <- lapply(fineLasso, function(z) {
              return(z$beta)
              })

        # Get objective values, if necessary
        objvals <- lapply(1:length(betas), function(i) {
              return(apply(betas[[i]], 2, function(z){
                norm(signals[[i]]$ASigma %*% z + signals[[i]]$vecC, type = "2")}))
              })

        parameters <- fineGrid
      }else{
        # Perform Lasso on test dataset
        testLasso <- foreach(i = 1:length(signals)) %dopar% {
                  glmnet(signals[[i]]$ASigma, -signals[[i]]$vecC, 
                            intercept = FALSE, alpha = 1,
                            standardize = TRUE, penalty.factor = penalty,
                            lambda = paras[[i]])
                  }
         # Get solution                    }
        betas <- lapply(testLasso, function(z) {
              return(z$beta)
              })
      }

    }

    if (length(paras[[1]]) > 1){
      return(list(M = betas, para = parameters))
    }else {
       return(list(M = betas))
    }
}

problem <- function(signal, p, n){

  # Get the population covariance matrix and the Volatility matrix
  C <- signal$C
  Sigma_star <- signal$popSigma


  # Sampled data and covariance
  if (n == Inf) {
    Sigma <- Sigma_star
  }else {
    X <- sampleX(n, Sigma_star)
    signal$X <- scale(X)
  }

  signal$p <- p
  signal$n <- n

  return(signal)

}

setupRegressionPrb <- function(signal, X){

  p <- signal$p
  C <- signal$C

  # Compute the covariance matrix
  if (signal$n != Inf){
    Sigma <- getDataCovar(X)
  }else{
    Sigma <- signal$popSigma
  }

  # Commutaion matrix that transforms column vetorization of
  # a matrix to its row vectorization
  KPermute <- spMatrix(nrow = p ** 2, ncol = p ** 2, i = c(seq(1, p ** 2)),
                      j = c(sapply(seq(1, p),
                            function(x) {seq(x, p ** 2, by = p)})),
                      x = c(rep(1, p ** 2)))
  # Data matrix for regression
  ASigma <- as.matrix(kronecker(Sigma, diag(p)) +
                      kronecker(diag(p), Sigma) %*% KPermute)

  # Vectorize volatility matrix for regression
  vecC <- c(C)

  return(list(Sigma = Sigma, ASigma = ASigma, vecC = vecC))
}

# Split data to training and test 
train_test_split <- function(signals, split = 0.8){

  if (split == 1){
    signals <- lapply(signals, function(x){
      lapply(x, function(y){
        lapply(y, function(z){
          if (z$n == Inf){
            return(c(setupRegressionPrb(z), z))
          }else{
            data <- z$X
            z["X"] <- NULL
            return(c(setupRegressionPrb(z, data), z))
          }
        })
      })
    })
  }else{
    train_signals <- lapply(signals, function(x){
      lapply(x, function(y){
        lapply(y, function(z){
          if (z$n == Inf) {
            return(c(setupRegressionPrb(z), z))
          }else{
            train_ind <- sample(seq(1, z$n), size = split * z$n)
            train_data <- z$X[train_ind, ]
            test_data <- z$X[-train_ind, ]
            return(c(setupRegressionPrb(z, train_data), 
                     list(train = train_data, test = test_data), z))
          }
        })
      })
    })
    
    test_signals <- lapply(train_signals, function(x){
      lapply(x, function(y){
        lapply(y, function(z){
          # Delete entries added for training data
          z[c("Sigma", "ASigma", "vecC")] <- NULL
          if (z$n == Inf){
            return(c(setupRegressionPrb(z), z))
          }else{
            data <- z$test
            z[c("X", "train", "test")] <- NULL
            return(c(setupRegressionPrb(z, data), z))
          }
        })
      })
    })
    train_signals <- lapply(train_signals, function(x){
      lapply(x, function(y){
        lapply(y, function(z){
          if (z$n != Inf) {
            z[c("X", "train", "test")] <- NULL
          }
          return(z)
        })
      })
    })
  }
  
  if (split == 1){
    return(list(train = signals, test = signals))
  }else{
    return(list(train = train_signals, test = test_signals))
  }
}

setupProblems <- function(signals, p, n) {

  signals <- lapply(signals, problem, p, n)

  return(signals)
}

# Accept results only with an MIP gap of less than 5 % 
checkMIPGap <- function(results, k){
  #  Get the indices at which we get a feasible solution
  fInd <- which(sapply(results,
                function(x){if(!is.null(x$mipgap)) {
                  return(x$mipgap <= 0.05)
                }else{
                  return(FALSE)
                }}))

  # Keep only the feasible results
  resultsFeasible <- results[fInd]

  feasibleK <- k[fInd]

  return(list(fInd = fInd, results = resultsFeasible, k = feasibleK))
}

# Fix singular matrices using nearPD (not used)
fixInvertibility <- function(beta){

  p <- sqrt(length(beta))
  M <- matrix(beta, nrow = p)
  E <- eigen(M)
  if(all(abs(E$values)) < 1e-3){
    print(TRUE)
    if(all(E$values) >= 0){
      M <- nearPD(M, base.matrix = TRUE, conv.tol = 1e-4)$mat
    }else{
      M <- -nearPD(-M, base.matrix = TRUE, conv.tol = 1e-4)$mat
    }
  }

  return(c(M))
}

getMetricsAlongRegPath <- function(signals, results, p, n, elapsedTime, methodDir) {

  betas <- results$M

  # Compute metrics for only off-diagonal elements
  diag <- seq(1, p ** 2, by = p + 1)
  # Evaluate Metrics
  results$max_acc <- mean(sapply(1:length(betas), function(x){
    max(apply(betas[[x]], 2, function(y){
      acc(y[-diag], signals[[x]]$gt[-diag])
    }))
  }))
  results$max_f1 <- mean(sapply(1:length(betas), function(x){
    max(apply(betas[[x]], 2, function(y){
      score <- f1score(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    }))
  }))
  results$tpr <- lapply(1:length(betas), function(x){
    apply(betas[[x]], 2, function(y){
      score <- tpr(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    })
  })
  results$fpr <- lapply(1:length(betas), function(x){
    apply(betas[[x]], 2, function(y){
      score <- fpr(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    })
  })
  results$precision <- lapply(1:length(betas), function(x){
    apply(betas[[x]], 2, function(y){
      score <- precision(y[-diag], signals[[x]]$gt[-diag])
      if (is.nan(score)) score <- 1
      return(score)
    })
  })

  # Area under the Curve along the regularization path
  results$aucroc <- mean(sapply(1:length(betas), function(x){
    AUCROC_path(results$tpr[[x]], results$fpr[[x]])
  }))
  results$aucpr <- mean(sapply(1:length(betas), function(x){
    AUCPR_path(results$precision[[x]], results$tpr[[x]])
  }))

  results$elaspsedTime <- elapsedTime
  # Save the results in desired location
    # save(results,
    #   file = paste0(getwd(), "/../Results/", methodDir, "/allModels_edgeProb5/", n, "_", p, "_res.RData"))
}

saveResults <- function(signals, results, p, n, elapsedTime, methodDir) {
  
      betas <- results$M
      paras <- results$para

      # Choose the best drift matrix based on posterior model probabilities
      Siginv_comp <- lapply(1:length(betas), function(x) {
        return(apply(betas[[x]], 2, invCov, signals[[x]]$vecC))
        })
      bic_scores <- lapply(1:length(betas), function(x) {
        return(apply(Siginv_comp[[x]], 2, ebic_score, signals[[x]]$Sigma, signals[[x]]$n, 1))
        }) 
      posterior_prob <- lapply(bic_scores, postprb)
      bestM <- lapply(1:length(betas), function(x) {
        return(b <- betas[[x]][, which.max(posterior_prob[[x]])])
      })

      # Choose the best drift matrix based on objective values
      # bestM <- betas[, which.min(objvals)]
      bestPara <- lapply(1:length(betas), function(x) {
        return(b <- paras[[x]][which.max(posterior_prob[[x]])])
      })
      bestk <- lapply(bestM, function(x) {return(sum(x != 0))})
      bestResult <- lapply(1:length(signals), function(x) {
        return(list(M_star = signals[[x]]$M, M_comp = bestM[[x]],
        bestPara = bestPara[[x]], bestk = bestk[[x]], k = signals[[x]]$k,
        t = elapsedTime))
      })

      # Save the results
      # save(bestResult,
      # file = paste0(getwd(), "/../Results/train/", n, "_", p, "_res.RData"))

      # Return, if required for further calculation of metrics
      return(bestResult)
}

getMetrics <- function(signals, bestResult, p, n, method, methodDir){

    paras <- lapply(bestResult, function(x){return(x$bestPara)})

    # Testing (skipped when population covariance matrix is used)
    if (n != Inf){
      testTime <- system.time(results <- recoverDriftMatrix(signals, p, method, paras, test = TRUE))
    }else{
      results$M <- lapply(bestResult, function(x){
        x$M_comp
      })
      testTime <-  0
    }

    # Evaluate Metrics
    diag <- seq(1, p ** 2, by = p + 1)
    bestResult <- apply(as.matrix(1:length(signals)), 1,  function(x) {

      if (!is.list(results$M[[x]])){
        bestResult[[x]]$M_testComp <- results$M[[x]]


        bestResult[[x]]$acc <- acc(bestResult[[x]]$M_testComp[-diag],
          signals[[x]]$gt[-diag])
        bestResult[[x]]$tpr <- tpr(bestResult[[x]]$M_testComp[-diag],
          signals[[x]]$gt[-diag])
        if (is.nan(bestResult[[x]]$tpr)) bestResult[[x]]$tpr <- 1
        bestResult[[x]]$fpr <- fpr(bestResult[[x]]$M_testComp[-diag],
          signals[[x]]$gt[-diag])
        if (is.nan(bestResult[[x]]$fpr)) bestResult[[x]]$fpr <- 1
        bestResult[[x]]$f1score <- f1score(bestResult[[x]]$M_testComp[-diag],
          signals[[x]]$gt[-diag])
        if (is.nan(bestResult[[x]]$f1score)) bestResult[[x]]$f1score <- 1
        bestResult[[x]]$aucroc <- AUCROC(bestResult[[x]]$M_testComp[-diag],
          signals[[x]]$gt[-diag])
        bestResult[[x]]$aucpr <- AUCPR(bestResult[[x]]$M_testComp[-diag],
          signals[[x]]$gt[-diag])
        return(bestResult[[x]])
      }else {
         return(NULL)
      }
    })

    res <- bestResult[which(sapply(bestResult, function(x){!is.null(x)}))]
    bestResult$InfeasibleResults <- length(which(sapply(bestResult, is.null)))

    # Get average over all signals
    bestResult$accAvg <- mean(unlist(lapply(res, function(x) {x$acc})))
    bestResult$tprAvg <- mean(unlist(lapply(res, function(x) {x$tpr})))
    bestResult$fprAvg <- mean(unlist(lapply(res, function(x) {x$fpr})))
    bestResult$f1scoreAvg <- mean(unlist(lapply(res, function(x) {x$f1score})))
    bestResult$aucrocAvg <- mean(unlist(lapply(res, function(x) {x$aucroc})))
    bestResult$aucprAvg <- mean(unlist(lapply(res, function(x) {x$aucpr})))
    bestResult$trainTimeAvg <- res[[1]]$t / length(signals)
    bestResult$testTimeAvg <- testTime / length(signals)

    # Save the results in desired location
    # save(bestResult,
      # file = paste0(getwd(), "/../Results/test/", n, "_", p, "_res.RData"))

}