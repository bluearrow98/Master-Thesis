rm(list=ls())

# Load required packages
library(gurobi)
library(MASS)
library(glmnet)
source("discreteFirstOrderAlgos.R")

# Load data
load("ex1data.Rdata")
load("ex1coef.Rdata")
n = dim(data)[1]; p = dim(data)[2] - 1; k0 = sum(b)
X = data[,1:p]
y = data[,p+1]
true = which(b!=0)
k = k0



algo2time = system.time(algo2Sol <- algorithm2(X,y,k))

#algo1time = system.time(algo1Sol <- algorithm1(X,y,k))


# Polishing active coefficients 
Xj = X[,which(algo2Sol$bm!=0)]
lm.fit = lm(y ~ 0 + Xj, data = data.frame(Xj,y))
bm_p = rep(0,p)
bm_p[which(algo2Sol$bm!=0)] = lm.fit$coefficients

# MIO 
Ml = 2*norm(as.matrix(bm_p),type="o")
Mu = norm(as.matrix(bm_p),type="i")

## Creating quadratic terms of objective
Q = matrix(data = 0, nrow=2*p,ncol=2*p)
Q[1:p,1:p] = 0.5*t(X)%*%X

## Creating inequality terms for max norm
ineqConst = rep(0,2*p*p)
mat = diag(p)
for (i in 1:p) {
  ineqConst[(2*(i-1)*p + 1):(p*(2*i-1))] = mat[i,]
}

## Initial guess (Not used for cold start)
z0 = rep(1,p)
z0[bm_p!=0] = 0
x0 = c(bm_p,z0)

model = list()
params = list()

## Gurobi Model creation 

# Warm-start model creation (Eq 2.5)
# model$A = matrix(c(rep(0,p),rep(1,p),rep(1,p),rep(0,p),ineqConst,ineqConst),ncol=2*p,byrow=T)
# model$Q = Q
# model$obj = c(-t(X)%*%y,rep(0,p))
# model$rhs = c(p-k0,Ml,rep(Mu,p),rep(-Mu,p))
# model$sense = c('>','<',rep('<',p),rep('>',p))
# model$vtype = c(rep('C',p),rep('B',p))
# model$start = x0

# Cold-start model creation (Eq 2.4)
model$A = matrix(c(rep(0,p),rep(1,p)),ncol=2*p)
model$Q = Q
model$obj = c(-t(X)%*%y,rep(0,p))
model$rhs = c(p-k)
model$sense = c('>')
model$vtype = c(rep('C',p),rep('B',p))
#model$start = x0


## Creating SOS constraints
sos = list()
for (i in 1:p){
  sos[[i]] = list()
  sos[[i]]$type = 1
  sos[[i]]$index = c(i,p+i)
  sos[[i]]$weights = c(1,2)
}
model$sos = sos

## Parameters for Gurobi optimization
params$OutputFlag = 1
params$timeLimit = 500

result = gurobi(model,params)

# Print number of non-zeros in the solution obtained and calculate the predictive error of the solution
print(sum(result$x[1:p]!=0))
mioError = norm(X%*%result$x[1:p] - X%*%b,type="2")**2/norm(X%*%b,type="2")**2

# Lasso 
l1 = glmnet(X,y,intercept = FALSE)
l1.cv = cv.glmnet(X,y,intercept = FALSE)
bhat = coef(l1,s=l1.cv$lambda.min)
lassoError = norm(X%*%as.matrix(bhat)[2:101] - X%*%b,type="2")**2/norm(X%*%b,type="2")**2
