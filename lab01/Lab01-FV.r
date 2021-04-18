s = 8
n = 24
f = 16
alpha_0 = beta_0 = 3


NDraws = 10000
res = matrix(,NDraws,3)
colnames(res) = c("iter","mean","sd")
alpha_n=alpha_0+s
beta_n = beta_0+f
TrueMean = (alpha_n)/(alpha_n+beta_n)
TrueVar = sqrt((alpha_n*beta_n)/((alpha_n+beta_n)^2*(alpha_n+beta_n+1)))
TrueMean
TrueVar

PostB <- function(alpha= alpha_n,beta =beta_n ,n=NDraws){
    set.seed(43)

    for(i in 1:n){
        postB = rbeta(i,alpha,beta)
        res[i,] = c(i,mean(postB),sd(postB))

    }
    
    return(list(post=postB,mean_sd=res))
}
Posterior = PostB()

Posterior$mean_sd[1,3]=0
head(Posterior$mean_sd)

x <-seq( from =-5, to =10 , by =0.001)
hist(Posterior$post,col='lightyellow',freq = FALSE,breaks=50,xlab='posterior',main='Histogram of Posterior')
curve(dbeta(x,alpha_n,beta_n),add=TRUE,col='blue', lwd =2,lty=3)


plot(x=1:NDraws,y=Posterior$mean_sd[,2],col='lightblue',lwd =1,type='l',xlab="Number of Draws",ylab=expression(mu))
abline(h=TrueMean,col='deeppink3',lwd=3)

plot(x=1:NDraws,y=Posterior$mean_sd[,3],col='lightblue',lwd =1,type='l',xlab="Number of Draws",ylab=expression(sigma))
abline(h=TrueVar,col='deeppink3',lwd=3)


exact_prb =1- pbeta(0.4,alpha_n,beta_n)
pos_prob = length(Posterior$post[Posterior$post > 0.4])/length(Posterior$post)
prb = data.frame(exact_prb,pos_prob)
colnames(prb) = c('excact_prob',"simulated_value")
prb


phi = log(Posterior$post/(1- Posterior$post))
hist(phi,breaks=100,probability = TRUE,col='cyan',xlab=expression(phi))
lines(density(phi),col='deeppink3',lty=2,lwd=3)

num_draws = 10000 #sample size
mu = 3.8
observ = c(38,20, 49, 58, 31, 70, 18, 56, 25,78)
n_obser = length(observ)


tau2 = function(data,mu,n_obser){
    sum((log(data) -mu)^2 )/n_obser
}

tau2(observ,mu,n_obser )

rinvchisq <- function(num_draws, n_obser, tau_sq){
  set.seed(1234)
  x <- rchisq(num_draws,df = n_obser-1)
  x_inv <- ((n_obser-1)*tau_sq)/x
  return(x_inv)
}

post_sigma2 <- function(m){
    set.seed(12345)
    rinvchisq(num_draws = num_draws,n_obser = n_obser,tau_sq = tau2(observ,mu,n_obser=n_obser))
}

dinvchisq <- function(x,n_obser,tau_sq){
  res <- (((tau_sq*(n_obser-1))/2)^(n_obser-1)/2 * exp(-(tau_sq*(n_obser-1))/(2*x)))/((x^(1+(n_obser-1)/2)) * gamma((n_obser-1)/2))
  return(res)
}

x <-seq( from =0, to =10 , by =0.001)
hist(post_sigma2(m = num_draws ),probability = TRUE,col='pink',breaks = 100,
     main=expression(paste('Simulated and theoratical posterior of ', sigma^2) ),
    ,xlab = expression(paste(sigma ^2)))
    curve(dinvchisq(x ,n_obser= n_obser,tau_sq = tau2(observ,mu,n_obser )),add=TRUE,col='blue',lwd=2)

sigma2 = post_sigma2(m)
g_sigma = sqrt(sigma2)/sqrt(2) 
gpdf =2* pnorm(q = g_sigma,mean = 0,sd = 1)-1
hist(gpdf,probability = TRUE,col='yellow',breaks=50,main='the posterior
distribution of the Gini coefficient',xlab=expression(paste('Gini ', sigma)))
lines(density(gpdf),col='darkblue',lty=3,lwd=2)

alpha_conf = 0.1
q_lower = quantile(gpdf,alpha_conf/2)
q_upper = quantile(gpdf,1-alpha_conf/2)
c(q_lower, q_upper)
true_mean = mean(gpdf)


true_mean

KerDenG=density(gpdf)
xx = KerDenG$x
yy = sort(KerDenG$y,decreasing = TRUE,index.return=TRUE)


CumSum = cumsum(yy$x)
HPD_D = CumSum[length(CumSum)]*0.9
temp = which(CumSum<HPD_D)
HPD_interval = range(xx[which(CumSum<HPD_D)])

df = data.frame('HPD_interval'=HPD_interval,'credible_interval'=c(q_lower, q_upper))


y_temp=density(gpdf,n=10000)
plot(x=y_temp$x,y=y_temp$y,type = 'l', lwd = 2, col = 'violet',main= 'distribution of the Gini coefficient'
     ,xlab=expression(paste('Gini ', sigma)),ylab='Density')
polygon(x=y_temp$x,y=y_temp$y, col = 'lightgrey',border = 'darkblue')
abline(v = true_mean, col="deeppink", lwd=3, lty=2)

segments(x0 =q_lower,y0 =0.5,x1 =q_upper,y1 =  0.5,col='deepskyblue1' ,lwd=3)
abline(v=q_lower,col='deepskyblue1' ,lwd=3)
abline(v=q_upper,col='deepskyblue1' ,lwd=3)
segments(x0=HPD_interval[1],y0=0.8,x1=HPD_interval[2],y=0.8,col='darkolivegreen',lwd=3)
abline(v=HPD_interval[1],col='darkolivegreen',lwd=3)
abline(v=HPD_interval[2],col='darkolivegreen',lwd=3)
arrows(x0 =q_lower,y0 =0.5,x1 =q_upper,y1 =  0.5,col='deepskyblue1' ,lwd=3,code = 3)
arrows(x0=HPD_interval[1],y0=0.8,x1=HPD_interval[2],y=0.8,col='darkolivegreen',lwd=3,code=3)
legend("topright", 
  legend = c(" Credible_interval", "HPD",expression(mu)), 
  col = c('deepskyblue1','darkolivegreen','deeppink'), lty=c(1,1,2), cex=0.8,lwd=3,bg='lightyellow')
                

y_data =c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

n =length(y_data)
mu = 2.39
n

# a function to compute prior for k (exponential with lambda = 1)
prior_k <- function(k){
    dexp(x = k,rate = 1)
}

likelikood_k <- function(y,mu ,n,k){
    
    exp(k*sum(cos(y-mu)))/((2*pi*besselI(k,nu = 0))^n)
}

k= seq(from = 0,to = 20,by = 0.01)

#
posterior_k <- function(y=y_data,mu=mu ,n=n,k=k){
    likelikood_k(y,mu ,n,k )*prior_k(k)
}


post=posterior_k(y=y_data,mu=mu ,n=n,k=k)

plot(y=post,x=k,col='blue',type='l',lwd=3)


k[which.max(post)]

hist(y_data,probability = TRUE,breaks = 50,col='brown')

y_data[which(y_data<0)]
