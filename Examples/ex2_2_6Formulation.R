rm(list=ls())

# Load required packages
library(gurobi)
library(MASS)
library(glmnet)

# Load data
load("ex2data.Rdata")
load("ex2coef.Rdata")
n = dim(data)[1]; p = dim(data)[2] - 1; k0 = sum(b)
X = data[,1:p]
y = data[,p+1]
true = which(b!=0)
k = 5

# Discrete First-order Algorithm 2 Implementation 
algorithm2 <- function(X,y,k,inits = 50,p=dim(X)[2]) {
  bm = matrix(rep(0,p*inits),ncol = p)
  gmPlus = rep(0,inits)
  L = norm(t(X)%*%X,type = "2")
  for (i in 1:inits) {
    randk = k
    randIndex = sample(1:p,randk)
    bm[i,randIndex] = min(i-1,1)*rnorm(randk,rep(0,randk),2*rep(1,randk))
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
algorithm1 <-function(X,y,k,inits = 150, p = dim(X)[2]) {
  bm = matrix(rep(0,p*inits),ncol = p)
  gmPlus = rep(0,inits)
  L = norm(t(X)%*%X,type = "i")
  for (i in 1:inits) {
    randIndex = sample(1:p,k)
    bm[i,randIndex] = min(i-1,1)*rnorm(k,rep(0,k),2*rep(1,k))
    gmPlus[i] = 0.5*norm(X%*%bm[i,] - y, type = "2")**2
    gm = 0
    iter = 1
    while (abs(gm - gmPlus[i]) >= 1e-4) {
      bmPlus = c(rep(0,p))
      H = bm[i,] - t(X)%*%(X%*%bm[i,] - y)/L
      bmPlus[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])] = 
        H[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])]
      gm = gmPlus[i]
      gmPlus[i] = 0.5*norm(X%*%bmPlus - y, type = "2")**2
      bm[i,] = bmPlus
      iter = iter + 1 
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


#MIO
Mu = 1.5*norm(as.matrix(bm_p),type="i")
Mu_zeta = max(apply(X,1,function(x){sum(as.matrix(sort(abs(x),decreasing = TRUE))[1:k])}))*Mu


## Creating quadratic terms of objective
Q = spMatrix(nrow=n + 2*p + 2,ncol=n + 2*p + 2,i = c(seq(1,n)),
             j = c(seq(1,n)), x = c(rep(0.5,n)))

## Initial guess (Not used for cold start)
z0 = rep(0,p)
z0[bm_p!=0] = 1
x0 = c(X%*%bm_p,bm_p,1-z0,Mu,Mu_zeta)

model = list()
params = list()

# Gurobi Model creation 

# Warm-start model creation (Eq 2.6)
model$A = spMatrix(nrow = n + 3, ncol = n + 2*p + 2, 
                   i = c(seq(1,n),rep(seq(1,n),each = p),rep(n+1,p),n+2,n+3),
                   j = c(seq(1,n),rep(seq(n+1,n+p),n),seq(n + p+1,n + 2*p),n + 2*p+
                           seq(1,2)),
                   x = c(rep(1,n),as.vector(-t(X)),rep(1,p+2)))
model$Q = Q
model$obj = c(rep(0,n),-t(X)%*%y,rep(0,p+2))
model$rhs = c(rep(0,n),p-k,Mu,Mu_zeta)
model$ub = c(rep(Inf,n),rep(Mu,p),rep(1,p),rep(Inf,2))
model$lb = c(rep(-Inf,n),rep(-Mu,p),rep(0,p),rep(0,2))
model$sense = c(rep('=',n),'>=',rep('<=',2))
model$vtype = c(rep('C',n),rep('C',p),rep('B',p),rep('C',2))
model$start = x0
model$genconnorm = list(list(resvar = n + 2*p + 1,
                             vars = c(seq(n+1,n + p)),
                             which = Inf),
                        list(resvar = n + 2*p + 2,
                             vars = c(seq(1,n)),
                             which = Inf))

## Creating SOS constraints
sos = list()
for (i in 1:p){
  sos[[i]] = list()
  sos[[i]]$type = 1
  sos[[i]]$index = c(n+i,n+p+i)
  sos[[i]]$weights = c(1,2)
}
model$sos = sos

## Parameters for Gurobi optimization
params$OutputFlag = 1
params$timeLimit = 100

result = gurobi(model,params)