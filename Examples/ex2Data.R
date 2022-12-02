rm(list=ls())

# Load required packages
library(gurobi)
library(MASS)
library(glmnet)

# Model generation (Example 2 - Synthetic dataset)
set.seed(8)
n = 30; p = 2000; k0 = 5
X = mvrnorm(n,mu = rep(0,p),Sigma = diag(p))
mean = matrix(rep(t(apply(X,2,function(x){mean(x)})),n),ncol = p, byrow = T)
X_shift = X - mean
X_norm = X_shift%*%diag(1/apply(X_shift, 2, function(x){norm(x,type="2")}))


# True coefficients
b = rep(0,p)
b[1:k0] = 1

# Noise 
SNR = 4
sigma = sqrt(var(X_norm%*%b)/SNR)
epsilon = rnorm(n,0,sigma)

# Response
y = X_norm%*%b + epsilon

# Saving data
data = matrix(c(X_norm,y),ncol = p+1)
save(data,file = "ex2data.Rdata")
save(b,file= "ex2coef.Rdata")