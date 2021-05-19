## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)


## -----------------------------------------------------------------------------
rainfall_data <- read.table('rainfall.dat')


## -----------------------------------------------------------------------------
# Data and prior Initialization
log_y <- as.matrix(log(rainfall_data))
mu0 <- 0
tau0_2 <- 1
v0 <- 1
si0_2 <- 1

# Initial value
si_2 <- 1
n_iterations <- 1000
mu_vec <- rep(NA,n_iterations)
si_vec <- rep(NA,n_iterations)
logy_post <- rep(NA,n_iterations)
n <- length(log_y)

# Pre calculation
log_y_mean <- mean(log_y)

# Function for getting samples from inverse chisquare distributino
rinvchisquare <- function(num_draws, n, tau_sq){
  #set.seed(1234)
  x <- rchisq(num_draws,df = n)
  x_inv <- ((n)*tau_sq)/x
  return(x_inv)
}

# Gibbs sampling iterations
for(k in 1:n_iterations){
  # Sampling mu
  w <- (n/si_2)/((n/si_2) + (1/tau0_2))
  taun_2 <- w/(n/si_2)
  mun <- (w * log_y_mean) + ((1 - w) * mu0)
  mu <- rnorm(1, mun, sqrt(taun_2))

  # Sampling si_2
  vn <- v0+n
  sin_2 <- (v0*si0_2 + sum((log_y - mu)^2))/vn
  si_2 <- rinvchisquare(1, vn, sin_2)
  
  # sample from posterior predictive
  logy_post[k] <- rnorm(1, mu, sqrt(si_2))
  
  mu_vec[k] <- mu
  si_vec[k] <- si_2
}



## -----------------------------------------------------------------------------
library(MASS)

rho_mu <- acf(mu_vec)
rho_si <- acf(si_vec)

IF_mu <- 1 + 2*sum(rho_mu$acf)
IF_si <- 1 + 2*sum(rho_si$acf)

cat('Inefficiency factor for parameter mu : ',IF_mu,'\n')
cat('Inefficiency factor for parameter mu : ',IF_si)

par(mfrow = c(1,2))
plot(x = 1:n_iterations, y = mu_vec,type = 'l', col = 'lightblue3',
     xlab = 'Iterations', ylab = 'mu')
plot(x = 1:n_iterations, y = si_vec,type = 'l', col = 'lightblue3',
     xlab = 'Iterations', ylab = 'si')


## -----------------------------------------------------------------------------
y_post <- exp(logy_post)
y_prior <- as.matrix(rainfall_data)

plot(density(y_prior),lwd = 4, main = 'Daily precipitation vs Posterior Predictive')
lines(density(y_post),col = 'orange',lwd = 2)
legend(300, 0.020,
       legend = c('Prior','Posterior Predictive'),
       col = c('black','orange'),
       lty=c(1,1),lwd = 3)



## -----------------------------------------------------------------------------
ebay_data <- read.table('ebayNumberOfBidderData.dat',header = TRUE)
rows <- nrow(ebay_data)
cols <- ncol(ebay_data)
y <- as.matrix(ebay_data[1])
X <- as.matrix(ebay_data[,2:cols])


## -----------------------------------------------------------------------------
x_glm <- X[,2:ncol(X)]
model <- glm(y ~ x_glm,family = poisson())
print(model$coefficients)


## -----------------------------------------------------------------------------
library(mvtnorm)
params <- dim(X)[2]
mu <- as.matrix(rep(0,params))
Sigma = 100*solve(t(X)%*%X)

LogPostPoisson <- function(betas,y,X,mu,Sigma){
  n <- nrow(X)
  lamda <- exp(X%*%betas)
  logLik <- (sum((X%*%betas)*y) - sum(lamda))
  Prior <- dmvnorm(betas, mu, Sigma, log=TRUE)
  return(logLik + Prior)
}

initValue <- matrix(0,params,1)

# Optimum beta(coefficient) are calculated
OptimRes <- optim(initValue,
                  LogPostPoisson, gr=NULL, y,X, mu, Sigma, method=c("BFGS"),
                  control=list(fnscale=-1), hessian=TRUE)


## -----------------------------------------------------------------------------
print('Posterior Mode: ')
print(OptimRes$par)
print('Inverse of hessian matrix')
inversehessian <- solve(-OptimRes$hessian)
print(inversehessian)


## -----------------------------------------------------------------------------
logPosteriorPoisson <- function(betas){
  params <- ncol(X)
  mu <- as.matrix(rep(0,params))
  Sigma = 100*solve(t(X)%*%X)
  
  n <- nrow(X)
  lamda <- exp(X%*%betas)
  logLik <- (sum((X%*%betas)*y) - sum(lamda))
  Prior <- dmvnorm(t(betas), mu, Sigma, log=TRUE)
  return(logLik + Prior)
}

sampleFromPosterior <- function(num_iterations, theta, c, postDensityFun){
  theta_current <- theta
  samples <- matrix(theta_current, num_iterations, nrow(theta))
  accept <- c()
  accept[1] <- 1
  for(i in 2:num_iterations){
    #shift <- as.vector(rmvnorm(1, mean = samples[i-1], sigma = c*inversehessian))
    #theta_new <- samples[i-1,] + shift
    theta_new <- as.vector(rmvnorm(1, mean = samples[i-1,], sigma = c*inversehessian))
    
    p_log_target_val <- postDensityFun(theta_new)
    log_target_val <- postDensityFun(samples[i-1,])
    
    alpha = min(1, exp(p_log_target_val - log_target_val))
    if(runif(1) <= alpha){
      #theta_current <- theta_new
      samples[i,] <- theta_new
      accept[i]<-1
    }
    else{
      samples[i,] <- samples[i-1,]
      accept[i]<-0
    }
  }
  print(sum(accept)/n_iterations)
  return(samples)
}
set.seed(1234)
beta_init <- matrix(0,dim(X)[2],1)
iter = 1000
sample_bet <- sampleFromPosterior(iter, beta_init, 0.8, logPosteriorPoisson)

for(i in 1:ncol(sample_bet)){
  plot(c(1:iter),sample_bet[,i],'l',
       ylab = paste('beta',i), xlab = 'Num Samples')
}


## -----------------------------------------------------------------------------
input <- c(1,1,1,1,0,1,0,1,0.7)
pred <- c()
n_samp <- 1000
for(i in 1:n_samp){
  pred[i] <- rpois(1,lambda = input*colMeans(sample_bet))
}
hist(pred, breaks = 100, xlab = 'Number of Bids',
     main = 'Posterior Predictive Distribution')

zero_bid <- length(pred[pred==0])
zero_bid_prob <- zero_bid/n_samp
print(paste('Probability of no bidders in the new auction is:',zero_bid_prob))


## -----------------------------------------------------------------------------
library(rstan)

mu = 20
sigma2 = 4
nT = 200 

AR.1 <- function(phi,mu=20,sigma2=4,nT=200) {
  
  x=numeric(nT)
  x[1]=mu + rnorm(n = 1,sd = sqrt(sigma2),mean = 0)
  for (i in 2:nT) {
    x[i]=mu+phi*(x[i-1]-mu)+ rnorm(1,0,sqrt(sigma2))
  }
  return(x)
}

phi <- c(-0.75, -.4,0, 0.4, 0.75,1)
set.seed(12345)
simulation = data.frame(sapply(phi, AR.1))

par(mfrow=c(2,3))

for(i in 1:length(simulation)){
  
  plot(simulation[,i],x=1:nT,type='l',col=i+1,
       xlab='Time',ylab='X_t',
       main=paste(expression(phi),"=",phi[i]))
}





## -----------------------------------------------------------------------------
#b

x <- AR.1(phi=0.3)
y <- AR.1(phi=0.9)

StanModel='
data {
  int<lower=0> nT; 
  vector[nT] x;
}
parameters {
  real mu;
  real<lower=0> sigma2;
  real<lower=-1, upper=1> phi;
}
model {
  mu ~ normal(0,1000);
  sigma2 ~ scaled_inv_chi_square(1000,2); //sigma2 initial value is 4
  x[2:nT] ~ normal( mu + phi*(x[1:(nT-1)]-mu), sqrt(sigma2)); //
}'

fit.x <- stan(
  model_code = StanModel,  
  data = list(N = length(x), d = x),    # named list of data
  chains = 1,
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  refresh = 0
)
fit.y <- stan(
  
  model_code = StanModel,  
  data = list(N = length(y), d = y),    # named list of data
  chains = 1,
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  refresh = 0
)

fit.x
fit.y

xDraws = extract(fit.x)
yDraws = extract(fit.y)


## -----------------------------------------------------------------------------
plot(fit.x)
plot(fit.y)


## -----------------------------------------------------------------------------
trace1 <- traceplot(fit.x,  pars = c("mu", "phi"), 
                    inc_warmup = TRUE,col='deeppink',include = TRUE)

trace2 <- traceplot(fit.y, pars = c("mu", "phi"), 
                    inc_warmup = TRUE,col='darkblue',include = TRUE)
trace1
trace2



## -----------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(x=xDraws$mu, y=xDraws$phi,
     xlab=expression(mu),
     ylab=expression(phi),
     main=expression(paste(phi,"="," 0.3")),
     col='darkolivegreen',pch=20)

plot(x=yDraws$mu, y=yDraws$phi,
     xlab=expression(mu), ylab=expression(phi),
     main=expression(paste(phi,"="," 0.9")),
     col='darkolivegreen',pch=20)


