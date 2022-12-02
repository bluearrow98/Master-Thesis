# Discrete First-order Algorithm 2 Implementation 

algorithm2 <- function(X,y,k,inits = 50,p=dim(X)[2],polish = T,maxIter = 1000) {
  
  bm = matrix(rep(0,p*inits),ncol = p)
  gmPlus = rep(0,inits)
  L = norm(t(X)%*%X,type = "2")

  for (i in 1:inits) {
    randk = sample(1:k,1)
    randIndex = sample(1:p,randk)
    bm[i,randIndex] = rnorm(randk,rep(0,randk),2*rep(1,randk))
    etam = rep(0,p)
    H = bm[i,] - t(X)%*%(X%*%bm[i,] - y)/L
    etam[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])] = 
      H[which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])]
    gm = 0
    gmPlus[i] = 0.5*norm(X%*%etam - y, type = "2")**2
    iter = 1
    while ((gm - gmPlus[i])/abs(gmPlus[i]) >= 1e-4)  {
      # Solution to min g(lambda*eta_m + (1-lambda)*beta_m)
      lambda = solve(t(etam - bm[i,])%*%t(X)%*%X%*%(etam-bm[i,]),t(etam - bm[i,])%*%t(X)%*%(y - X%*%bm[i,]))
      bm[i,] = c(lambda)*etam + c((1-lambda))*bm[i,]
      H = bm[i,] - t(X)%*%(X%*%bm[i,] - y)/L
      etam = H
      etam[-which(abs(H) %in% sort(abs(H),decreasing=T)[1:k])] = 0
      gm = gmPlus[i]
      gmPlus[i] = 0.5*norm(X%*%etam - y, type = "2")**2
      iter = iter + 1
      if (iter > maxIter) break
    }
    
    # Send to Algorithm 1 
    gm = 0
    iter = 1 
    while ((gm - gmPlus[i])/abs(gmPlus[i]) >= 1e-4) {
      bm[i,] = bm[i,] - t(X)%*%(X%*%bm[i,] - y)/L
      bm[i,-which(abs(bm[i,]) %in% sort(abs(bm[i,]),decreasing=T)[1:k])] = 0
      gm = gmPlus[i]
      gmPlus[i] = 0.5*norm(X%*%bm[i,] - y, type = "2")**2
      iter = iter + 1 
    }
    
    if (polish == T) {
      Xj = X[,which(bm[i,]!=0)]
      fitCoef = lsfit(Xj,y,intercept = F)$coefficients
      bm[i,which(bm[i,]!=0)] = fitCoef
      gmPlus[i] = 0.5*norm(X%*%bm[i,] - y, type = "2")**2
    }
  }
  # Take the coefficients with the best objective
  bm = bm[which.min(gmPlus),]
  
  
  return (list("bm" = bm, "gm" = min(gmPlus)))
}



# Discrete First-order Algorithm 1 Implementation 
algorithm1 <-function(X,y,k,inits = 50, p = dim(X)[2],polish = T,maxIter = 1000) {
  
  bm = matrix(rep(0,p*inits),ncol = p)
  gmPlus = rep(0,inits)
  L = norm(t(X)%*%X,type = "2")
  
  for (i in 1:inits) {
    randk = sample(1:k,1)
    randIndex = sample(1:p,randk)
    bm[i,randIndex] = rnorm(randk,rep(0,randk),2*rep(1,randk))
    gmPlus[i] = 0.5*norm(X%*%bm[i,] - y, type = "2")**2
    gm = 0
    iter = 1
    while (abs(gm - gmPlus[i])/abs(gmPlus[i]) >= 1e-4) {
      bm[i,] = bm[i,] - t(X)%*%(X%*%bm[i,] - y)/L
      bm[i,-which(abs(bm[i,]) %in% sort(abs(bm[i,]),decreasing=T)[1:k])] = 0
      gm = gmPlus[i]
      gmPlus[i] = 0.5*norm(X%*%bm[i,] - y, type = "2")**2
      iter = iter +1
      if (iter > maxIter) break
    }
    
    if (polish == T) {
      Xj = X[,which(bm[i,]!=0)]
      fitCoef = lsfit(Xj,y,intercept = F)$coefficients
      bm[i,which(bm[i,]!=0)] = fitCoef
      gmPlus[i] = 0.5*norm(X%*%bm[i,] - y, type = "2")**2
    }
    
  }
  # Take the coefficients with the best objective
  bm = bm[which.min(gmPlus),]
  
  return (list("bm" = bm, "gm" = min(gmPlus)))
}


proj_grad = function(X, y, k, inits=50, maxiter=1000,polish=TRUE) {
  
  n = nrow(X)
  p = ncol(X)
  
  # If beta0 is NULL, use thresholded least squares coefficients when p < n,
  # and thresholded marginal regression coefficients when p >= n
  if (p < n) beta0 = lsfit(X,y,int=FALSE)$coef
  else beta0 = t(X)%*%y/colSums(X^2)
  ids = order(abs(beta0), decreasing=TRUE)
  beta0[-ids[1:k]] = 0
  
  # If L is NULL, use the power method to approximate the largest eigenvalue
  # of X^T X, for the step size
  L = norm(t(X)%*%X,type = "2")
  
  beta.beta = beta0
  best.crit = Inf
  beta = beta0
  
  for (r in 1:nruns) {
    for (i in 1:maxiter) {
      beta.old = beta
      
      # Take gradient descent step
      grad = -t(X) %*% (y - X %*% beta)
      beta = beta - grad/L
      
      # Set to zero all but the top k
      ids = order(abs(beta), decreasing=TRUE)
      beta[-ids[1:k]] = 0
      
      # Perform least squares polishing, if we are asked to
      if (polish) beta[ids[1:k]] = lsfit(X[,ids[1:k]],y,int=FALSE)$coef
      
      # Stop if the relative difference in coefficients is small enough
      if (norm(beta - beta.old) / max(norm(beta),1) < tol) break
    }
    
    # Compute the criterion for the current coefficients, compare to the
    # best so far
    cur.crit = sum((y - X%*%beta)^2)
    if (cur.crit < best.crit) {
      best.crit = cur.crit
      best.beta = beta
    }
    
    # Start the next run off at a random spot (particular choice matches
    # Rahul's Matlab code)
    beta = beta0 + 2*runif(p)*max(abs(beta0),1)
  }
  
  obj = 0.5*norm(X%*%best.beta - y, type = "2")**2
  
  return(list("bm" = best.beta, "gm" = obj))
}