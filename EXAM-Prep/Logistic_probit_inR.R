#logisticRegression 

Probit = 1
wineData = data.frame(readRDS("WineData.rds"))
Nobs = dim(wineData)[1]
y = wineData[,dim(wineData)[2]]
y = ifelse(y > 5, yes = 1, no = 0)
Cov = c(1:11)

standardize <- TRUE # If TRUE, covariates/features are standardized to mean 0 and variance 1

lambda <- 1# scaling factor for the prior of beta 

# add intercept to our data
wineData = data.frame(rep(1,Nobs),wineData)

#the covariants - matrix : note x should be a matrix

X = as.matrix(wineData[,Cov])
Xnames <- colnames(X)

#Let's standardize X
#X = apply(X, 2, scale)

if(standardize){
  
  Indx = 2:(length(Cov )-1)
  X[,Indx] = scale(X[,Indx])
}

Npar = dim(X)[2]

# now we can set up the prior

mu = numeric(Npar) # Prior mean vector
Sigma = (1/lambda)*diag(Npar) # Prior covariance matrix


# now the Posterior Functions

# Functions that returns the log posterior for the logistic and probit regression.
# First input argument of this function must be the parameters we optimize on, 
# i.e. the regression coefficients beta.

logPostLogist <- function(beta,y,X,mu,Sigma){
  
  library(mvtnorm)
  linPred = X %*% beta
  loglike = sum((y*linPred) - log(1 + exp(linPred)))
  logPrior = dmvnorm(beta,mu,Sigma,log = TRUE)
  return( loglike + logPrior)
  
}

LogPostProbit <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  SmallVal <- .Machine$double.xmin
  logLik <- sum(y*log(pnorm(linPred)+SmallVal) + (1-y)*log(1-pnorm(linPred)+SmallVal) )
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  return(logLik + logPrior)
}

#now the optimization and finding the best Beta

#first define Init
InitVal = matrix(rep(0,Npar))

if (Probit==1){
  logPost = LogPostProbit;
} else{
  logPost = logPostLogist;
}

# The argument control is a list of options to the optimizer optim, where fnscale=-1 means that 
# we minimize the negative log posterior. Hence, we maximize the log posterior.  

OptimRes = optim(InitVal,fn = logPost,y,X,mu,Sigma,method = c("BFGS"),
                 control = list(fnscale=-1),gr = NULL,hessian = TRUE)
names(OptimRes$par) <- Xnames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(-solve(OptimRes$hessian))) # Computing approximate standard deviations.
names(approxPostStd) <- Xnames # Naming the coefficient by covariates
print('The posterior mode is:')
print(OptimRes$par)
print('The approximate posterior standard deviation is:')
approxPostStd <- sqrt(diag(-solve(OptimRes$hessian)))
print(approxPostStd)