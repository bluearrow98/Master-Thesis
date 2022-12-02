library("Rlab")
library(MASS)
library(Matrix)
library(parallel)

stabSignal <- function(p){
  k = sample(c(1,2,3,4),size = 1)
  M = matrix(rnorm(p**2)*rbern(p**2,prob = 1/p), nrow = p)
  nonDiagAbsSum = sum(abs(M)) - sum(diag(abs(M)))
  diag(M) = - nonDiagAbsSum - abs(rnorm(p))
  return (M)
}

getPopCovar <- function(M,C){
  p = dim(M)[1]
  Sigma = solve(kronecker(diag(p),M) + kronecker(M,diag(p)),-c(C))
  Sigma = matrix(Sigma,nrow = p)
  
  return (Sigma)
}

sampleX <- function(n,Sigma){
  X = mvrnorm(n,mu = rep(0,dim(Sigma)[1]), Sigma)
  return (X)
}

getDataCovar <- function(X){
  n = dim(X)[1]
  Sigma = t(X)%*%X/n
  
  return (Sigma)
}

standardize <- function(X){
  n = dim(X)[1]
  p = dim(X)[2]
  mean = matrix(rep(t(apply(X,2,function(x){mean(x)})),n) ,ncol = p, byrow = T)
  X_shift = X - mean
  X_norm = X_shift%*%diag(1/apply(X_shift, 2, function(x){norm(x,type="2")}))
  return (X_shift)
}


customBS <- function(X,y,k,polish = T){
  
  p = dim(X)[2]
  
  data = standardize(cbind(X,y))
  X = data[,1:p]
  y = data[,p+1]
  
  # Discrete First order method
  algo2time = system.time(algo2Sol <- proj_grad(X,y,k,polish = polish))
  
  M_p = algo2Sol$bm
  
  #MIO
  Mu = 1.5*norm(as.matrix(M_p),type="i")
  
  
  ## Initial guess (Not used for cold start)
  z0 = rep(0,p)
  z0[M_p!=0] = 1
  x0 = c(M_p,1-z0,norm(as.matrix(M_p),type="i"))
  
  model = list()
  params = list()
  
  # Warm-start model creation (Eq 2.5)
  model$A = spMatrix(nrow = 2, ncol = 2*p+ 1, i = c(rep(1,p),2),
                     j = c(seq(p+1,2*p+1)),
                     x = c(rep(1,p+1)))
  # model$A = spMatrix(nrow = 2*p+1, ncol = 2*p,
  #                    i= c(rep(seq(1,p),2),rep(seq(p+1,2*p),2),rep(2*p+1,p)),
  #                    j=c(rep(seq(1,2*p),2),seq(p+1,2*p)),
  #                    x=c(rep(-1,p),rep(-Mu,p),rep(1,p),rep(-Mu,p),rep(1,p)))
  model$Q = spMatrix(nrow=2*p+1,ncol=2*p+1,i = c(rep(seq(1,p),each = p)),
                     j = c(rep(seq(1,p),p)), x = c(0.5*t(X)%*%X))
  model$obj = c(-t(X)%*%y,rep(0,p+1))
  model$objcon = 0.5*t(y)%*%y
  model$ub = c(rep(Mu,p),rep(1,p),Inf)
  model$lb = c(rep(-Mu,p),rep(0,p),0)
  model$rhs = c(p-k,Mu)
  model$sense = c('>=','<=')
  model$vtype = c(rep('C',p),rep('B',p),'C')
  model$start = x0
  model$genconnorm = list(list(resvar = 2*p + 1,
                               vars = c(seq(1,p)),
                               which = Inf))

  ## Creating SOS constraints
  sos = list()
  for (j in 1:p){
    sos[[j]] = list()
    sos[[j]]$type = 1
    sos[[j]]$index = c(j,p+j)
    sos[[j]]$weights = c(1,2)
  }
  model$sos = sos
  
  ## Parameters for Gurobi optimization
  params$OutputFlag = 0
  params$timeLimit = 100
  params$NonConvex = 2
  
  result = gurobi(model,params)
  
  result$algo2Time = algo2time
  result$algo2Sol = algo2Sol
  
  return(result)
}

applyBS <- function(j,ASigma,vecC) {
  
  # p = length(vecC)
  # result = customBS(ASigma,-vecC,j)
  
  result = bs(ASigma,-vecC,j,intercept = F,
              params = list(NonConvex = 2))
  
  # Store objVal and constructed signal
  # result$beta = result$x[1:p]
  result$objVals = 0.5*norm(ASigma%*%as.matrix(result$beta) + vecC,
                            type = "2")**2
  
  return (result)
  
}

recoverDriftMatrix <- function(i,p,n,numCores,method) {
  
  set.seed(i)
  
  # True signal and Volatility matrix
  M_star = stabSignal(p)
  trueSparsity = (M_star != 0)
  C = diag(runif(p))
  # Population covariance
  Sigma_star = getPopCovar(M_star,C)  
  # Sampled data and covariance
  if (n == Inf){
    Sigma = Sigma_star
  }else{
    X = sampleX(n, Sigma_star)
    Sigma = getDataCovar(X)
  }
  # Data matrix of regression
  KPermute = spMatrix(nrow = p**2, ncol = p**2, i = c(seq(1,p**2)),
                      j = c(sapply(seq(1,p),function(x){seq(x,p**2,by=p)})),
                      x = c(rep(1,p**2)))
  ASigma = as.matrix(kronecker(Sigma,diag(p)) + kronecker(diag(p),Sigma)%*%KPermute)
  vecC = c(C)
  
  # True number of best subset
  k0 = sum(M_star != 0)
  
  
  # Storing the best result after hyperparameter search
  bestResult = list()
  
  if (method == "bs") {
    # Number of best subset to select
    k = seq(k0-k0%/%2,k0 + k0%/%2)
    
    # Pre-allocation 
    resultBS = list()
    
    # resultBS = parLapply(numCores,k, applyBS, ASigma = ASigma,vecC = vecC)
    foreach(i = k) %dopar% {applyBS(i,ASigma,vecC)} -> resultBS
    
    # Evaluate Metrics 
    # resultBS = lapply(resultBS,function(x){x$tpr = tpr(x$beta,M_star)
    # x$fpr = fpr(x$beta,M_star)
    # x$acc = acc(x$beta,M_star)
    # x$f1score = f1score(x$beta,M_star)
    # return (x)})
    
    bestM = resultBS[[which.min(unlist(lapply(resultBS,function(x){x$objVals})))]]$beta
    bestResult = resultBS[[which.min(unlist(lapply(resultBS,function(x){x$objVals})))]]
  }
  else if (method == "lasso") {
    # Lasso
    penalty = rep(1,p**2)
    penalty[seq(1,p**2,by=p)] = 0
    l1 = glmnet(ASigma,-vecC,intercept = FALSE,alpha = 1,standardize = T,penalty.factor = penalty)
    l1.cv = cv.glmnet(ASigma,-vecC,intercept = FALSE,alpha = 1,standardize = T, penalty.factor = penalty)
    bestM = coef(l1,s=l1.cv$lambda.min)[2:(p**2+1)]
    
  }
  
  bestk = sum(bestM!=0)
  bestResult$M_star = M_star
  bestResult$M_comp = bestM
  bestResult$bestk = bestk
  bestResult$k = k0
  
  # Evaluate Metrics  
  tp = length(intersect(which(bestM!=0),which(M_star!=0)))
  fp = length(setdiff(which(bestM!=0),which(M_star!=0)))
  tn = length(intersect(which(bestM==0),which(M_star==0)))
  fn = length(setdiff(which(bestM==0),which(M_star==0)))
  
  
  bestResult$acc = (tp + tn)/(tp + tn + fp + fn)
  bestResult$tpr = tp/(tp + fp)
  bestResult$fpr = fp/(fp + tn)
  bestResult$f1score = tp/(tp + 0.5*(fp + fn))
  
  
  return (bestResult)
}


tp <- function(comp,gt){ return (length(intersect(which(comp!=0),which(gt!=0))))}
fp <- function(comp,gt){return (length(setdiff(which(comp!=0),which(gt!=0))))}
tn <- function(comp,gt){return (length(intersect(which(comp==0),which(gt==0))))}
fn <- function(comp,gt){return (length(setdiff(which(comp==0),which(gt==0))))}
acc <- function(comp,gt){
  
  accuracy <- (tp(comp,gt) + tn(comp,gt)) /(tp(comp,gt) + tn(comp,gt) 
                                            + fp(comp,gt) + fn(comp,gt))
  return (accuracy)}
tpr <- function(comp,gt){
  return (tp(comp,gt)/(tp(comp,gt) + fp(comp,gt)))
}
fpr <- function(comp,gt){
  return (fp(comp,gt)/(fp(comp,gt) + tn(comp,gt)))
}
f1score <- function(comp,gt){
  return (tp(comp,gt)/(tp(comp,gt) + 0.5*(fp(comp,gt) + fn(comp,gt))))
}

sigmoid <- function(x){
  return (exp(abs(x))/(1+exp(abs(x))))
}