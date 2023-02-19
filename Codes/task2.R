rm(list=ls())

# Import required packages
library(parallel)
library(doParallel)
library(bestsubset)
source("helpFunctions.R")
source("discreteFirstOrderAlgos.R")

# Parallel settings
# numCores = makeCluster(detectCores())
# registerDoParallel(numCores)
# clusterExport(numCores,c("bs"))
# clusterEvalQ(numCores,source("helpFunctions.R"))
# clusterEvalQ(numCores,source("discreteFirstOrderAlgos.R"))
# clusterEvalQ(numCores,library(glmnet))
# clusterEvalQ(numCores,library(gurobi))

recoverDriftMatrix <- function(i, p, n, method) {

set.seed(i)

# True signal and Volatility matrix
M_star <- stabSignal(p) 
trueSparsity <- (M_star != 0)
C <- diag(runif(p))
# Population covariance
Sigma_star <- getPopCovar(M_star,C)  
# Sampled data and covariance
if (n == Inf){
  Sigma <- Sigma_star
}else{
  X <- sampleX(n, Sigma_star)
  Sigma <- getDataCovar(X)
}
# Data matrix of regression
KPermute <- spMatrix(nrow = p**2, ncol = p**2, i = c(seq(1,p**2)),
                    j = c(sapply(seq(1,p),function(x){seq(x,p**2,by=p)})),
                    x = c(rep(1,p**2)))
ASigma <- as.matrix(kronecker(Sigma,diag(p)) + kronecker(diag(p),Sigma)%*%KPermute)
vecC <- c(C)

# True number of best subset
k0 <- sum(M_star != 0)

# Storing the best result after hyperparameter search
bestResult <- list()

if (method != "tLasso") {

  if (method == "bs") {
    # Number of best subset to select
    # k = seq(k0-k0%/%2,k0 + k0%/%2)
    k <- floor(seq(p,p**2,length.out = 10))


    resultBS = lapply(k, applyBS, ASigma = ASigma,vecC = vecC)
    # foreach(i = k) %dopar% {applyBS(i,ASigma,vecC)} -> resultBS

    bestM <- resultBS[[which.min(unlist(lapply(resultBS,function(x){x$objval})))]]$beta
    bestResult <- resultBS[[which.min(unlist(lapply(resultBS,function(x){x$objval})))]]
  }
  else if (method == "lasso") {
    # Lasso
    penalty <- rep(1,p**2)
    penalty[seq(1,p**2,by=p+1)] <- 0

    initGrid <- seq(10, 10^-5, length.out = 100)
    coarseLasso <- glmnet(ASigma, -vecC, intercept = FALSE, alpha = 1, 
                            standardize = T, penalty.factor = penalty,
                            lambda = initGrid)
    # Find min lambda such that M is diagonal
    lambda <- min(initGrid[(colSums(penalty*(coarseLasso$beta!=0))==0 & 
                              colSums((1-penalty)*(coarseLasso$beta!=0)) == p)])
    fineGrid <- seq(lambda, lambda/10^4,length.out = 100)
    fineLasso <- glmnet(ASigma,-vecC, intercept = FALSE, alpha = 1, 
                        standardize = T, penalty.factor = penalty,
                        lambda = fineGrid)
    Siginv_comp <- apply(fineLasso$beta, 2, invCov, vecC)
    bic_scores <- apply(Siginv_comp, 2, bic_score, Sigma, n)
    posterior_prob <- postprb_bic(bic_scores)
    bestM <- fineLasso$beta[, which.max(posterior_prob)]
  }

  bestk <- sum(bestM!=0)
  bestResult$M_star <- M_star
  bestResult$M_comp <- bestM
  bestResult$bestk <- bestk
  bestResult$k <- k0

  # Evaluate Metrics  

  bestResult$acc <- acc(bestM, trueSparsity)
  bestResult$tpr <- tpr(bestM, trueSparsity)
  bestResult$fpr <- fpr(bestM, trueSparsity)
  bestResult$f1score <- f1score(bestM, trueSparsity)
  bestResult$aucroc <- AUCROC(bestM, trueSparsity)
  bestResult$aucpr <- AUCPR(bestM, trueSparsity)

  return(bestResult)

}else {
  # Save Dataset to be later loaded in Matlab 
writeMat(con = paste0("D:/TUM/Master-Thesis/Dataset/",n,"_",p,"/X_",n,"_",p,"_",i,".mat"),X = ASigma)
writeMat(con = paste0("D:/TUM/Master-Thesis/Dataset/",n,"_",p,"/y_",n,"_",p,"_",i,".mat"),y = -vecC)
writeMat(con = paste0("D:/TUM/Master-Thesis/Dataset/",n,"_",p,"/M_",n,"_",p,"_",i,".mat"),M = M_star)
}
}

# Signal size
p <- 5

# Number of signals 
N <- 100


# Method for best subset
method <- "lasso"

set.seed(99)


for (s in p) {

  # Number of data to be sampled
  n <- c(100,200,500,1000,10000,100000,Inf)


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
         file = paste0(getwd(),"/../Results/LassoResults/",d,"_",s,"_res.RData"))

  }
}

# stopCluster(numCores)