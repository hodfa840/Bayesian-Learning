TempLink = read.table("TempLinkoping.txt",header = TRUE)
attach(TempLink)

dim(TempLink)

head(TempLink)

# setting the initial values
mu.0 = c(-10,100,10)
omega.0 = 0.01*diag(3)
nu.0 = 4
sigma2.0 = 1



tau2<- function(data,mu,n){
    sum((log(data)-mu)^2)/n
}
# Random generation from a scaled inverse chisquare
rinvchisq <- function(draws, n, tau) {
chi_square <- rchisq(draws, n)
return( tau*(n-1)/chi_square )
}
# Density of a scaled inverse chisquare
dinvchisq <- function(data, df, tau) {
return( (tau2*df/2)^(df/2)/gamma(df/2) * exp(-df*tau2/(2*data)) / data^(1+df/2) )
}


lmTemp = lm( temp ~ time + I(time^2), data = TempLink)

summary(lmTemp)

sqrt(26.7)

sum(lmTemp$residuals**2)/length(lmTemp$residuals)

plot(y=temp,x=time,col='deeppink',pch=19,lwd=3)
lines(y=lmTemp$fitted.values,x=time,col='blue',lwd=3)

library(mvtnorm )


sigma2.prior <- function(){
    rinvchisq(draws = 1,n = nu.0,tau = sigma2.0)
}

Beta.prior <- function(sigma2){
    rmvnorm(mean =mu.0,n=1,sigma = sigma2*solve(omega.0))
}

# create empty structure for sigma,Beta and error
NDraws = 200
ErrorTerm = numeric(NDraws)
sigma2 = numeric(NDraws)
BetaList = matrix(,NDraws,3)
colnames(BetaList) = c('B0','B1','B2')


for(i in 1:NDraws){
    sigma2[i] = sigma2.prior()
    BetaList[i,1] = Beta.prior(sigma2[i])[1]
    BetaList[i,2] = Beta.prior(sigma2[i])[2]
    BetaList[i,3] = Beta.prior(sigma2[i])[3]
    ErrorTerm[i] = rnorm(1,mean = 0,sd = sqrt(sigma2))
}

Bayes.Regressor = matrix(,length(time),NDraws)

for(i in 1:NDraws){
    Bayes.Regressor[,i]= BetaList[i,1] +BetaList[i,2]*time + BetaList[i,3]*(time^2) + ErrorTerm[i]
}
colnames(Bayes.Regressor)=paste0('model',1:NDraws)

head(data.frame(Bayes.Regressor),1)

TempLink1=data.frame(TempLink)
models=data.frame(TempLink1,Bayes.Regressor)

# PLOT for start value of parameters
library(ggplot2)
plot_start_param= ggplot(models , aes(y=temp,x = time)) +
  labs(title =expression(paste("Linkoping Temperature"," " ,mu[0],',',sigma[0],",",
                               Omega[0],',',nu[0])),x = "Time", y="Tempreture") + theme_minimal()

for(i in names(models)[-c(1,2)]){
  plot_start_param = plot_start_param +
    geom_line(aes_string(y = i), color="blue", alpha=0.2)
}

plot_start_param = plot_start_param +
  geom_point(aes(y = temp), alpha=0.5,color='deeppink')
plot_start_param

lmTemp$coefficients

mu.0 = c(-11, 103,-95)
omega.0 = 0.01*diag(3)
nu.0 = 10
sigma2.0 =  0.03
NDraws=100
Bayes.Regressor2 = matrix(,length(time),NDraws)
ErrorTerm2 = numeric(NDraws)
sigma22 = numeric(NDraws)
BetaList2 = matrix(,NDraws,3)
colnames(BetaList2) = c('B0','B1','B2')
for(i in 1:NDraws){
    sigma22[i] = sigma2.prior()
    BetaList2[i,1] = Beta.prior(sigma22[i])[1]
    BetaList2[i,2] = Beta.prior(sigma22[i])[2]
    BetaList2[i,3] = Beta.prior(sigma22[i])[3]
    ErrorTerm[i] = rnorm(1,mean = 0,sd = sqrt(sigma22))
}
for(i in 1:NDraws){
    Bayes.Regressor2[,i]= BetaList2[i,1] +BetaList2[i,2]*time + BetaList2[i,3]*(time^2) + ErrorTerm2[i]
}
TempLink2=data.frame(TempLink)
models2=data.frame(TempLink2,Bayes.Regressor2)
plot_new_param= ggplot(models2 , aes(y=temp,x = time)) +
  labs(title =expression(paste("Linkoping Temperature revised value for "," " ,mu[0],',',sigma[0],",",
                               Omega[0],',',nu[0])),x = "Time", y="Tempreture") + theme_minimal()

for(i in names(models2)[-c(1,2)]){
  plot_new_param = plot_new_param +
    geom_line(aes_string(y = i), color="blue", alpha=0.2)
}

plot_new_param = plot_new_param +
  geom_point(aes(y = temp), alpha=0.5,color='deeppink')
plot_new_param

X=model.matrix(lmTemp)
y = temp
# setting the initial values
mu.0 =  c(-11, 103,-95)
omega.0 = 0.01*diag(3)
nu.0 = 4
sigma2.0 = 1
n = dim(X)[1]
NDraws = 1000

omega.n = t(X) %*% X + omega.0
nu.n = nu.0 + n -3 
betaHat = solve(t(X) %*% X) %*% t(X) %*% y
mu.n = solve(t(X) %*% X + omega.0) %*% (t(X) %*% X %*% betaHat + omega.0 %*% mu.0)
sigma2.n = (nu.0 * sigma2.0 + (t(y) %*% y + t(mu.0) %*% omega.0%*% mu.0 - t(mu.n) %*% omega.n %*% mu.n)) / nu.n

sigma2n.pos = numeric(NDraws)
BetaList2n.pos = matrix(,NDraws,3)
colnames(BetaList2n.pos)=c('B0','B1','B2')
pos.sigma2 <- function(nu.n,sigma2.n){
    rinvchisq(1,n = nu.n,tau = sigma2.n)
    
}
pos.Beta <- function(sigma2_n,mu_n,omega_n){
    rmvnorm(1, mu_n, solve(omega_n)*as.numeric(sigma2_n))
}

for (i in 1:NDraws ){
    sigma2n.pos[i] = pos.sigma2(nu.n,sigma2.n)
    BetaList2n.pos[i,] = pos.Beta(sigma2_n = sigma2.n,mu_n = mu.n,omega_n = omega.n)
}

par(mfrow=c(2,2))
p1 =hist(BetaList2n.pos[,1],breaks = 20,probability = TRUE,
         xlab = expression(beta[0]),col='khaki1',main=expression(paste('Histogram of'," " ,beta[0])))
lines(density(BetaList2n.pos[,1]),lwd=3,col='blue')
p2 = hist(BetaList2n.pos[,2],breaks = 20,probability = TRUE,
          xlab = expression(beta[2]),col='khaki1',main=expression(paste('Histogram of'," " ,beta[1])))
lines(density(BetaList2n.pos[,2]),lwd=3,col='blue')

p3 = hist(BetaList2n.pos[,3],breaks = 20,probability = TRUE,
          xlab = expression(beta[3]),col='khaki1',main=expression(paste('Histogram of'," " ,beta[2])))
lines(density(BetaList2n.pos[,3]),lwd=3,col='blue')

p4= hist(sigma2n.pos,probability=TRUE,col='khaki1',
         xlab = expression(sigma[n]),breaks=20,,main=expression(paste('Histogram of'," " ,sigma[n])))
lines(density(sigma2n.pos),lwd=3,col='blue')

Beta.median = apply ( BetaList2n.pos , 2, median )

f.time.median = Beta.median %*% t(X)

length(f.time.median)

# Estimation of the whole dataset with 1000 different beta parameters
ypost = BetaList2n.pos %*% t(X)

CI <- matrix(, n, 2)
colnames(CI) <- c("lower","upper")

# calculate the 95% credible interval
for(i in 1:n){
CI[i,] <- quantile( ypost[,i], probs = c(0.025,0.975))
}

plot(y=temp,x=time,lwd=2,pch=19,col = 'darkblue')
lines ( y = CI[ ,2] ,x= time , col= " brown ",lwd=3,lty = 3)
lines ( y = CI[ ,1] ,x= time , col= " brown ",lwd=3, lty = 3)
lines (y=f.time.median  ,x= time ,col ="darkolivegreen",lwd=2)
legend("topright",c('median','Credible interval'),lty=c(1:3),lwd=c(2,3),col=c('darkolivegreen','brown'),bg='lightgrey')

dim(ypost)

highest.Temp = numeric(n)
highest.Temp = apply(ypost,2,max)

X.tilda <- (-1*BetaList2n.pos[,2])/(2*BetaList2n.pos[,3])
hist(X.tilda,probability = TRUE,breaks = 20,col='lightblue')
lines(density(X.tilda),lwd=3,lty=3,col='darkviolet')

library(mvtnorm)
ww_data <- read.table('WomenWork.dat', header = TRUE)
rows <- nrow(ww_data)
cols <- ncol(ww_data)
y <- as.matrix(ww_data[1])
X <- as.matrix(ww_data[,2:cols])

params <- dim(X)[2]
mu <- as.matrix(rep(0,params))
tau = 10
Sigma = (tau^2)*diag(params)

LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum(linPred*y - log(1 + exp(linPred)))
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE)
  
  return(logLik + logPrior)
}

initValue <- matrix(0,params,1)

# Optimum beta(coefficient) are calculated
OptimRes <- optim(initValue,
                  LogPostLogistic, gr=NULL, y, X, mu, Sigma, method=c("BFGS"),
                  control=list(fnscale=-1), hessian=TRUE)

print('Posterior Mode: ')
print(OptimRes$par)
print('Inverse of hessian matrix')
inversehessian <- solve(OptimRes$hessian)
print(inversehessian)

print('Estimtes obtained using GLM model:')
print(glm(Work ~ 0 + ., data = ww_data, family = binomial)$coefficients)

beta_post <- OptimRes$par
names(beta_post) <- colnames(ww_data[,2:cols])
approx_par_NSC <- rnorm(1000,beta_post['NSmallChild'],-inversehessian[7,7])
lowerInterval <- quantile(approx_par_NSC,0.05)
upperInterval <- quantile(approx_par_NSC, 0.95)

plot(density(approx_par_NSC),lwd = 3,main = '95% Posterior Probability Interval of NSmallChild variable')
polygon(density(approx_par_NSC), col = 'lightblue2')
abline(v = lowerInterval, col = 'red', lwd = 3)
abline(v = upperInterval, col = 'red', lwd = 3)
arrows(lowerInterval,0.3,upperInterval,0.3,length = 0.1,col = 'black',lwd = 3)
arrows(upperInterval,0.3,lowerInterval,0.3,length = 0.1,col = 'black',lwd = 3)
text(-1.4,0.5, '95% CI', lwd = 3,cex = 1.3)

posteriorPredictive <- function(x, beta_post, inversehessian){
  post_sample <- rmvnorm(1, beta_post,-inversehessian)
  logist_prob <- (exp(x %*% t(post_sample)))/(1 + exp(x %*% t(post_sample)))
  return(logist_prob)
}

x <- c(1,13, 8, 11, (11/10)^2, 37, 2, 0)
nsamples = 1000
post_predict <- c(rep(0,nsamples))
for(i in 1:nsamples){
  post_predict[i] <- posteriorPredictive(x, beta_post, inversehessian)
}
#print(post_predict)
hist(post_predict, breaks = 100, 
     col = 'lightblue2',lwd = 3,
     xlab = 'Pr(y=1|X)',
     main = 'Posterior Predictive Plot')
lines(density(post_predict), col = 'red',lwd = 3)
polygon(density(post_predict), col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5))

posteriorPredictiveBinomial <- function(x, beta_post, inversehessian){
  post_sample <- rmvnorm(1, beta_post,-inversehessian)
  logist_prob1 <- (exp(x %*% t(post_sample)))/(1 + exp(x %*% t(post_sample)))
  logist_prob <- sum(rbinom(1,8,logist_prob1))
  return(logist_prob)
}

test_data <- matrix(x,nrow=8, ncol=8,byrow = TRUE)
nsamples = 1000
post_predict <- c(rep(0,nsamples))
for(i in 1:nsamples){
  post_predict[i] <- posteriorPredictiveBinomial(test_data, beta_post, inversehessian)
}
hist(post_predict,breaks = 100,col = 'lightblue2',
     xlab = 'Number of women who works',
     main = 'Histogram of number of women who works ')
