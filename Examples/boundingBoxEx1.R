rm(list=ls())

# Load required packages
library(gurobi)
library(MASS)
library(glmnet)

# Load data
load("ex1data.Rdata")
load("ex1coef.Rdata")
n = dim(data)[1]; p = dim(data)[2] - 1; k0 = sum(b)
X = data[,1:p]
y = data[,p+1]
true = which(b!=0)
k = k0

# Discrete First-order Algorithm 2 Implementation 
algorithm2 <- function(X,y,k,inits = 50,p=dim(X)[2]) {
  bm = matrix(rep(0,p*inits),ncol = p)
  gmPlus = rep(0,inits)
  L = norm(t(X)%*%X,type = "2")
  for (i in 1:inits) {
    randk = sample(1:k,1)
    randIndex = sample(1:p,randk)
    bm[i,randIndex] = rnorm(randk,rep(0,randk),2*rep(1,randk))
    etam = c(rep(0,p))
    H = bm[i,] - t(X)%*%(X%*%bm[i,] - y)/L
    etam[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])] = 
      H[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])]
    gm = 0
    gmPlus[i] = 0.5*norm(X%*%etam - y,type = "2")
    iter = 1
    while ((abs(gm - gmPlus[i])) >= 1e-4) {
      # Solution to min g(lambda*eta_m + (1-lambda)*beta_m)
      lambda = solve(t(etam - bm[i,])%*%t(X)%*%X%*%(etam-bm[i,]),t(etam - bm[i,])%*%t(X)%*%(y - X%*%bm[i,]))
      bmPlus = c(lambda)*etam + c((1-lambda))*bm[i,]
      H = bmPlus - t(X)%*%(X%*%bmPlus - y)/L
      etam[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])] = 
        H[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])]
      gm = gmPlus[i]
      gmPlus[i] = 0.5*(norm(X%*%etam - y, type = "2"))**2
      bm[i,] = bmPlus
      iter = iter + 1
    }
  }
  # Take the coefficients with the best objective
  bm = bm[which.min(gmPlus),]
  
  # Send to Algorithm 1 
  gm = 0
  gmPlus = norm(X%*%bm - y,type = "2")**2/2
  iter = 1 
  while (abs(gm - gmPlus) >= 1e-4) {
    bmPlus = c(rep(0,p))
    H = bm - t(X)%*%(X%*%bm - y)/L
    bmPlus[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])] = 
      H[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])]
    gm = gmPlus
    gmPlus = 0.5*norm(X%*%bmPlus - y, type = "2")**2
    bm = bmPlus
    iter = iter + 1 
  }
  
  return (list("bm" = bm, "gm" = min(gmPlus)))
}

algo2time = system.time(algo2Sol <- algorithm2(X,y,k))


# Discrete First-order Algorithm 1 Implementation 
algorithm1 <-function(X,y,k,inits = 100, p = dim(X)[2]) {
  bm = matrix(rep(0,p*inits),ncol = p)
  gmPlus = rep(0,inits)
  L = norm(t(X)%*%X,type = "2")
  for (i in 1:inits) {
    randk = sample(1:k,1)
    randIndex = sample(1:p,randk)
    bm[i,randIndex] = rnorm(randk,rep(0,randk),2*rep(1,randk))
    gmPlus[i] = 0.5*norm(X%*%bm[i,] - y, type = "2")**2
    gm = 0
    k = 1
    while (abs(gm - gmPlus[i]) >= 1e-4) {
      bmPlus = c(rep(0,p))
      H = bm[i,] - t(X)%*%(X%*%bm[i,] - y)/L
      bmPlus[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])] = 
        H[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])]
      gm = gmPlus[i]
      gmPlus[i] = 0.5*norm(X%*%bmPlus - y, type = "2")**2
      bm[i,] = bmPlus
      k = k +1 
    }
  }
  # Take the coefficients with the best objective
  bm = bm[which.min(gmPlus),]
  
  return (list("bm" = bm, "gm" = min(gmPlus)))
}

#algo1time = system.time(algo1Sol <- algorithm1(X,y,k))


# Polishing active coefficients 
Xj = X[,which(algo2Sol$bm!=0)]
lm.fit = lm(y ~ 0 + Xj, data = data.frame(Xj,y))
bm_p = rep(0,p)
bm_p[which(algo2Sol$bm!=0)] = lm.fit$coefficients

# MIO

# Bounding box formulation (Eq 5.3) for high dimensional regime
# Bounds
Mu = 2*norm(as.matrix(bm_p),type="o")
Ml = k*Mu

model = list()
params = list()

## Creating equality terms
eqConst = rep(0,n*3*p)
for (i in 1:p) {
  eqConst[(1 + (i-1)*3*p):(p + (i-1)*3*p)] = rep(1,p)
  eqConst[(p + 1 + (i-1)*3*p):(2*p + (i-1)*3*p)] = -X[i,]
}