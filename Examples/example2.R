rm(list=ls())

# Load required packages
library(gurobi)
library(MASS)
library(glmnet)
source("discreteFirstOrderAlgos.R")
source("helpFunctions.R")

# Load data
load("ex2data.Rdata")
load("ex2coef.Rdata")
n = dim(data)[1]; p = dim(data)[2] - 1; k0 = sum(b)
X = data[,1:p]
y = data[,p+1]
true = which(b!=0)
k = 5

result = customBS(X,y,k)

# Print number of non-zeros in the solution obtained and calculate the predictive error of the solution
print(sum(result$x[1:p]!=0))
mioError = norm(X%*%result$x[1:p] - X%*%b,type="2")**2/norm(X%*%b,type="2")**2

# Lasso 
l1 = glmnet(X,y,intercept = FALSE)
l1.cv = cv.glmnet(X,y,intercept = FALSE)
bhat = coef(l1,s=l1.cv$lambda.min)


