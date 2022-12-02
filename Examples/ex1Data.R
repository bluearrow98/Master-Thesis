rm(list=ls())

# Load required packages
library(gurobi)
library(MASS)
library(glmnet)

# Model generation (Example 1 - Synthetic dataset)
set.seed(8)
n = 500; p = 100; k0 = 10; rho = 0.5
Sigma = matrix(nrow=p,ncol=p)
for (i in 1:p) {
  for (j in 1:p) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = mvrnorm(n,mu = rep(0,p),Sigma = Sigma)
mean = matrix(rep(t(apply(X,2,function(x){mean(x)})),n),ncol = p, byrow = T)
X_shift = X - mean
X_norm = X_shift%*%diag(1/apply(X_shift, 2, function(x){norm(x,type="2")}))

# True coefficients 
b = rep(0,p)
b[c(10,20,30,40,50,60,70,80,90,100)] = 1
true = which(b!=0)

# Noise 
sigma = 0.07
epsilon = rnorm(n,0,sigma)

# Response
y = X_norm%*%b + epsilon
SNR = var(X_norm%*%b)/sigma**2

# Saving data
data = matrix(c(X_norm,y),ncol = p+1)
save(data,file = "ex1data.Rdata")
save(b,file= "ex1coef.Rdata")