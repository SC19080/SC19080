## ------------------------------------------------------------------------
ad_data = read.csv("advertising.csv")  ##Import the data
lg = lm(sales~TV,data = ad_data)  ## Linear Regression
summary(lg)$coef

## ------------------------------------------------------------------------
knitr::kable(head(ad_data))    ##Transform to .Rmd

## ------------------------------------------------------------------------
par(mfrow=c(1,1))   
plot(lg)            

## ------------------------------------------------------------------------
ge_rayleigh = function(n, sigma){ 
  #generate n random samples from Rayleigh(sigma)
  
  V = runif(n)
  return(sigma * sqrt(-2*log(V)) )
}


## ------------------------------------------------------------------------
set.seed(12)

hist_rayleigh = function(n, sigma, y){
  # Draw a histogram to check
  # y is a sequence
  
  x = ge_rayleigh(n, sigma)                      # Generate the random sample
  
  hist(x, prob = T,main="Rayleigh Distribution")
  
  lines(y, y*exp(-y^2/(2*sigma^2))/(sigma^2))
}

hist_rayleigh(100000,1,seq(0,5,0.05))         # Histogram when sigma = 1

hist_rayleigh(100000,3,seq(0,15,0.15))        # Histogram when sigma = 3

hist_rayleigh(100000,0.2,seq(0,1,0.01))       # Histogram when sigma = 0.2




## ------------------------------------------------------------------------
set.seed(27)

ge_NLM = function(n, p1){
    # Generate a random smaple of size n from a normal location mixture, with p_1 = p1.  
    # I applied the (first) algorithm on Statistical Computing with R, Page 78.
  
    U = runif(n)                                     
    
    X = rep(0,n)
    
    for(i in 1:n){
      
      if(U[i] <= p1){
        
        X[i] = rnorm(1)              # Generate X from N(0,1)
      
      }
      else{
        
        X[i] = 3 + rnorm(1)          # Generate X from N(3,1)
      
      }
    }
    return(X)
}

hist_NLM = function(n, p1){
  # Generate the sample and graph the histogram
  
  X = ge_NLM(n, p1)
  
  hist(X, probability = T, main = "Normal Location Mixture")             # Graph the histogram    
  
  y = seq(-4,6,0.1)
  
  lines(y, (p1*exp(-1/2 * y^2) + (1-p1)*exp(-1/2 * (y-3)^2))/sqrt(2*pi))      
  # graph the probability density function
  
}

hist_NLM(1000, 0.75)


for(p in seq(0.1,0.9,0.1)){
  hist_NLM(1000,p)                       # Repeat with different values
}



## ------------------------------------------------------------------------
library(matrixcalc)

ge_Wishart = function(n, Sigma){
  # Generate a random sample from a W_d(sigma,n) distribution.
  
  if(is.numeric(n) == F | length(n) > 1){      # Check that n and sigma satiesfy the assumptions needed.
    return("n is not a number")
  }
  
  n = as.integer(n)
  
  if(is.matrix(Sigma) == F){
    return("Sigma is not a matrix")
  }
  if(isSymmetric(Sigma) == F){
    return("Sigma is not symmetric")
  }
  if(matrixcalc::is.positive.definite(Sigma) == F){
    return("Sigma is not positive definite")
  }
  
  d = nrow(Sigma)
  
  if(n <= d+1){
    return("error in dimension")                # Check the assumption mentioned in the question.
  }
  
  L = chol(Sigma)                  
  # Obtain the Choleski factorization. Note that L is upper triangular.
  
  TT = matrix(0, nrow = d, ncol = d)
  # Create the matrix
  
  for(i in 1:d){                                
  #applied the algorithm on Statistical Computing with R, Page 80
    
    TT[i,i] = sqrt(rchisq(1, n-i+1))
    
    if(i > 1){
      
      for(j in 1:(i-1)){
        
        TT[i,j] = rnorm(1)
      } 
    }
  }
  
  return(t(L) %*% TT %*% t(TT) %*% L)
}



## ------------------------------------------------------------------------

Sigma = cbind(c(4,12,-16),c(12,37,-43),c(-16,-43,98))       # An example of Sigma

ge_Wishart(5,Sigma)       # Generate a random sample

## ------------------------------------------------------------------------
set.seed(1011)

U = runif(10000, min = 0, max = pi/3)     # generate 10000 random samples of U

pi / 3 * mean(sin(U))                     # Calculate the estimation of the integral
pi / 3 * mean(sin(U)) - 1/2               # calculate the bias
var(pi /3 * sin(U))                       # Calculate the sample variance

## ------------------------------------------------------------------------
MC.antithetic = function(R = 1000, antithetic = T){   
  # Modified from the codes on page 131
  
  u = runif(R/2)
  
  if(!antithetic){
    v = runif(R/2)
  }
    else{
    v = 1 - u  
  } 
  
  u = c(u,v)
  
  E = exp(-u)/(1+u^2)
  
  return(mean(E))
}


MC.antithetic(antithetic = T)   # estimate with antithetic sampling


m=1000

MC1 = rep(0,1000)
MC2 = rep(0,1000)

for(i in 1:1000){
  MC1[i] = MC.antithetic(antithetic = F)
  MC2[i] = MC.antithetic(antithetic = T)
}

sd(MC1)                         # sample variance without antithetic sampling

sd(MC2)                         # sample variance with antithetic sampling

(var(MC1)-var(MC2))/var(MC1)    # reduction in variance


## ------------------------------------------------------------------------
m <- 10000 

g <- function(x){ 
  exp(-x - log(1+x^2)) * (x > 0) * (x < 1) 
}

u <- runif(m)                            #f3, inverse transform method 
x <- - log(1 - u * (1 - exp(-1))) 
fg <- g(x) / (exp(-x) / (1 - exp(-1))) 

theta.hat <- mean(fg)                    #estimate of importance sampling
se <- sd(fg)                             #estimated standard error

## ------------------------------------------------------------------------
m = 2000                                  #replicate 2000 times on each strata 

Quantile=c(0,0.1,0.2,0.3,0.4,1)

# the following codes are modified from Example 5.12

fg_im = matrix(0,nrow = 5, ncol = m)              #initialization

for(j in 1:5){
  u = runif(m,(j-1)/5,j/5)       #f3, inverse transform method
  x = - log(1 - u * (1 - exp(-1))) 
  fg_im[j,] = g(x) / (exp(-x) / (1 - exp(-1))) 
}
  
fg2 = colMeans(fg_im)

theta.hat2 = mean(fg2) 
se2 = sd(fg2)

result = rbind(c(theta.hat,se),c(theta.hat2,se2))
rownames(result) = c("Importance Sampling","Stratified Importance Sampling")
colnames(result) = c("estimate","estimated standard error")

result

(var(fg)-var(fg2))/var(fg)          # reduction in variance

## ------------------------------------------------------------------------
set.seed(1018)

m = 1e5                    #repeat 1e5 times
n = 20                     #sample size

estimate = rep(0,m)
estimate_2 = rep(0,m)

for(i in 1:m){
  
  X = rchisq(n, df = 2)       # generate random samples from  chi_2^2 
  
  estimate[i] = mean(X)- qt(0.975,n-1)*var(X)/sqrt(n) <= 2 & mean(X) + qt(0.975,n-1)*var(X)/sqrt(n) >= 2
  # notice that the mean of the chi_2^2 is 2
  # calculate whether the mean of chi_2^2 is 2 is covered by the t-interval in the ith simulation
  
  estimate_2[i] = (n-1)*var(X)/qchisq(0.05, df = n-1) >= 4
  # the variance of the chi_2^2 is 4
  # calculate whether the variance of chi_2^2 is covered by the confidence interval in Example 6.4
  
  }

mean(estimate)              # estimated coverage probability of the t-interval for mean
mean(estimate_2)            # estimated coverage probability of the confidence interval for variance in Example 6.4

## ------------------------------------------------------------------------
q = c(0.025,0.05,0.95,0.975)

sk = function(x){ 
  #computes the sample skewness coeff. 
  xbar = mean(x) 
  m3 = mean((x - xbar)^3) 
  m2 = mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) 
} 


Quantile_skew = function(mean = 0, sd = 1, n = 500, m = 1e4, q){
  # Generate the table of the answer
  
  # Input: mean and sd are mean and standard error of the normal distribution.
  # n is the sample size. m is the repeating times.
  # q is the sequence of the quantiles.
  
  skew = rep(0,m)
  
  for(i in 1:m){
  X = rnorm(n, mean = mean, sd = sd)          
  skew[i] = sk(X)                         # Calculating the sample skewness coefficients
  }
  
  a = quantile(skew, probs = q)           
  # the samples quantiles of the sample skewness coefficients          
  
  b_1 = qnorm(q, sd = sqrt(6/n) )         
  # theoretical quantiles of the sample skewness coefficients, with the first variance sqrt{6/n} applied
  
  b_2 = qnorm(c(0.025,0.05,0.95,0.975), sd = sqrt(6*(n-2)/(n+1)/(n+3)))
  # theoretical quantiles of the sample skewness coefficients, with the second variance applied
  
  c_1 = sqrt( q*(1-q)/n/( pnorm(b_1, sd = sqrt(6/n)) )^2 )
  # the standard error of estimates, with the first variance sqrt{6/n} applied
  
  c_2 = sqrt( q*(1-q)/n/( pnorm(b_2, sd = sqrt(6/n)) )^2 )
  # the standard error of estimates, with the second variance applied
  
  result = cbind(a,b_1,b_2,c_1,c_2)
  colnames(result) = c("sample quantile","theoretical quantile 1","theoretical quantile 2","standard error 1", "standard error 2")
  
  knitr::kable(result)
  
}

set.seed(123)
Quantile_skew(q = q)
Quantile_skew(n = 20, q=q)

## ------------------------------------------------------------------------
sk = function(x) { 
  #computes the sample skewness coeff. 
  xbar = mean(x) 
  m3 = mean((x - xbar)^3) 
  m2 = mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) 
} 


Pwr_Skew_NB = function(alpha, sl = 0.1, n = 30, m = 2500  ){
  # return the estimated power
  # alpha is the parameter in Beta distribution
  # sl is the significance level
  # n is the sample size
  # m is the replication times
  
  N = length(alpha) 

  pwr = numeric(N) 
 
  cv = qnorm(1-sl/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))      
  #critical value for the skewness test

  for (j in 1:N) {              #for each epsilon 
  
    shape = alpha[j] 
  
    sktests= numeric(m) 
  
    for (i in 1:m) {            #for each replicate 
    
      x = rbeta(n, shape1 = shape, shape2 = shape) 
    
      sktests[i] = as.integer(abs(sk(x)) >= cv) 
    
    } 
  
  pwr[j] = mean(sktests) 
  }
  
  return(pwr)
}

alpha = seq(0,2,0.1)
#plot power vs alpha
plot(alpha, Pwr_Skew_NB(alpha = alpha), type = "b", xlab = bquote(alpha), ylab = "power", ylim = c(0,0.1)) 
abline(h = .1, lty = 3)


alpha = seq(2,20,1)
#plot power vs alpha
plot(alpha, Pwr_Skew_NB(alpha = alpha), type = "b", xlab = bquote(alpha), ylab = "power", ylim = c(0,0.1)) 
abline(h = .1, lty = 3)


## ------------------------------------------------------------------------
Pwr_Skew_NT = function(alpha, sl = 0.1, n = 30, m = 2500  ){
  # return the estimated power
  # alpha is the parameter in t distribution
  # sl is the significance level
  # n is the sample size
  # m is the replication times
  
  N = length(alpha) 

  pwr = numeric(N) 
 
  cv = qnorm(1-sl/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))      
  #critical value for the skewness test

  for (j in 1:N) {              #for each epsilon 
  
    freedom = alpha[j] 
  
    sktests= numeric(m) 
  
    for (i in 1:m) {            #for each replicate 
    
      x = rt(n, df = freedom) 
    
      sktests[i] = as.integer(abs(sk(x)) >= cv) 
    
    } 
  
  pwr[j] = mean(sktests) 
  }
  
  return(pwr)
}

alpha=1:20
#plot power vs alpha
plot(alpha, Pwr_Skew_NT(alpha = alpha), type = "b", xlab = bquote(alpha), ylab = "power", ylim = c(0,1)) 
abline(h = .1, lty = 3) 

## ------------------------------------------------------------------------

Type_1_chisq = function(sl, df = 1, m = 1e4, n = 50){
  # reutrn the empirical Type I error rate of chi square distribution
  # sl is the significance level(sl can be a sequence)
  # df is the degree of freedom
  # m is the number of the replication times
  # n is the sample size
  
  N = length(sl)
  result = rep(0,N)
  
  for(j in 1:N){
    
    esl = rep(0,n) 

    for(i in 1:m){
      x = rchisq(n, df = 1)                                               # generate the random sample
      xbar = mean(x)
      q = qt(1-sl[j]/2, df = n-1)                                       # calculate the quantile
      S = var(x)
      esl[i] = mean(xbar-q*S/sqrt(n) > 1 || xbar+q*S/sqrt(n) < 1 )        
      # calculate empirical Type I error rate
      # the mean of chi^2(1) is 1
    }
    result[j] = mean(esl)
  }
  return(result)
}

sl = c(0.05,0.1,0.5)
M = rbind(sl, Type_1_chisq(sl))
rownames(M) = c("significance level","empirical type I error rate")
M

## ------------------------------------------------------------------------

Type_1_unif = function(sl, minimum = 0, maximum = 2, m = 1e4, n = 30){
  # reutrn the empirical Type I error rate of uniform distribution
  # sl is the significance level(sl can be a sequence)
  # minimum and maximum is the minimum and maximum of the uniform distribution
  # m is the number of the replication times
  # n is the sample size
  
  N = length(sl)
  result = rep(0,N)
  
  for(j in 1:N){
    
    esl = rep(0,n) 

    for(i in 1:m){
      x = runif(n, min = minimum, max = maximum)                                  # generate the random sample
      xbar = mean(x)
      q = qt(1-sl[j]/2, df = n-1)                                       # calculate the quantile
      S = var(x)
      esl[i] = mean(xbar-q*S/sqrt(n) > 1 || xbar+q*S/sqrt(n) < 1 )        
      # calculate empirical Type I error rate
      # the mean of chi^2(1) is 1
    }
    result[j] = mean(esl)
  }
  return(result)
}

sl = c(0.05,0.1,0.5)
M = rbind(sl, Type_1_unif(sl))
rownames(M) = c("significance level","empirical type I error rate")
M

## ------------------------------------------------------------------------
Type_1_exp = function(sl, rate = 1, m = 1e4, n = 30){
  # reutrn the empirical Type I error rate of uniform distribution
  # sl is the significance level(sl can be a sequence)
  # rate is the parameter in exponential distribution
  # m is the number of the replication times
  # n is the sample size
  
  N = length(sl)
  result = rep(0,N)
  
  for(j in 1:N){
    
    esl = rep(0,n) 

    for(i in 1:m){
      x = rexp(n, rate = rate)                                 # generate the random sample
      xbar = mean(x)
      q = qt(1-sl[j]/2, df = n-1)                                       # calculate the quantile
      S = var(x)
      esl[i] = mean(xbar-q*S/sqrt(n) > 1 || xbar+q*S/sqrt(n) < 1 )        
      # calculate empirical Type I error rate
      # the mean of chi^2(1) is 1
    }
    result[j] = mean(esl)
  }
  return(result)
}

sl = c(0.05,0.1,0.5)
M = rbind(sl, Type_1_exp(sl))
rownames(M) = c("significance level","empirical type I error rate")
M

## ------------------------------------------------------------------------
library(bootstrap)
k = ncol(bootstrap::scor)
names = colnames(bootstrap::scor)
par(mfrow = c(2, 3))
for(i in 1:(k-1)){
  for(j in (i+1):k){
    plot(bootstrap::scor[,i],bootstrap::scor[,j],sub = paste("correlation","=",cor(bootstrap::scor[,i],bootstrap::scor[,j])), xlab=paste(names[i]),ylab=paste(names[j]) )
  }
}

## ------------------------------------------------------------------------
par(mfrow = c(1, 1))


rho = function(x,i){
  return(cor(x[i,1],x[i,2]))
}

library(boot)

#boostrap estimate of cor(mec,vec)
sd(boot::boot(data = scor[,c(1,2)], statistic = rho, R=2000)$t)

#boostrap estimate of cor(alg,ana)            
sd(boot::boot(data = scor[,c(3,4)], statistic = rho, R=2000)$t)

#boostrap estimate of cor(alg,sta)
sd(boot::boot(data = scor[,c(3,5)], statistic = rho, R=2000)$t)

#boostrap estimate of cor(ana,sta)
sd(boot::boot(data = scor[,c(4,5)], statistic = rho, R=2000)$t)

## ------------------------------------------------------------------------
sk = function(x) { 
  #computes the sample skewness coeff. 
  xbar = mean(x) 
  m3 = mean((x - xbar)^3) 
  m2 = mean((x - xbar)^2) 
  return( m3 / m2^1.5 ) 
}

sk_boot= function(x,i){
  #function used for bootstrap 
  return( sk(x[i,]) )
}


sk_cp_norm = function(n = 30, m = 1000, R = 200, conf = 0.95){
  # return the four estimated coverage probabilities of the four intervals for nomal populations
  # m is the repeating times of this estimation
  # R is the number of bootstrap replicates
  # conf is the confidence level
  
  cp_norm = rep(0,m)
  cp_perc = rep(0,m)
  cp_basic = rep(0,m)
  
  
  for(i in 1:m){
    
    x = as.matrix(rnorm(n))
    
    obj = boot.ci(boot(data = x, statistic = sk_boot, R =R), conf = conf, type = c("norm", "basic", "perc"))
    
    # the skewness of a normal distribution is 0
    cp_norm[i] = obj$norm[2] <= 0 && obj$norm[3] >= 0
    cp_perc[i] = obj$perc[4] <= 0 && obj$perc[5] >= 0
    cp_basic[i] = obj$basic[4] <= 0 && obj$basic[5] >= 0
    
  }
  
  result =c(mean(cp_norm),mean(cp_perc),mean(cp_basic))
  return(result) 
}

sk_cp_chisq = function(n = 30, m = 1000, R = 200, conf = 0.95){
  
  cp_norm = rep(0,m)
  cp_perc = rep(0,m)
  cp_basic = rep(0,m)
  
  
  for(i in 1:m){
    
    x = as.matrix(rchisq(n, df = 5))
    
    obj = boot.ci(boot(data = x, statistic = sk_boot, R =R), conf = conf, type = c("norm", "basic", "perc"))
    
    # the skewness of chi square distribution  is 2*sqrt(2/degrees of freedom)
    cp_norm[i] = obj$norm[2] <= 2*sqrt(2/5) && obj$norm[3] >= 2*sqrt(2/5)
    cp_perc[i] = obj$perc[4] <= 2*sqrt(2/5) && obj$perc[5] >= 2*sqrt(2/5)
    cp_basic[i] = obj$basic[4] <= 2*sqrt(2/5) && obj$basic[5] >= 2*sqrt(2/5)
    
  }
  
  result = c(mean(cp_norm),mean(cp_perc),mean(cp_basic))
  return(result) 
}

compare = function(n = 30, m = 1000, R = 200, conf = 0.95){
  # compare the results and construct a table
  r1 = sk_cp_norm(n=n, m=m, R=R, conf = conf)
  r2 = sk_cp_chisq(n=n, m=m, R=R, conf = conf)
  result = rbind(r1,r2)
  rownames(result) = c("normal","chi sqaure")
  colnames(result) = c("normal","precentile","basic")
  knitr::kable(result)
}

compare(n = 30, m=1000, R = 200, conf = 0.95)
compare(n = 30, m=1000, R = 200, conf = 0.9)
compare(n = 30, m=1000, R = 200, conf = 0.5)

## ------------------------------------------------------------------------
theta2 = function(xdata){                     # function used to calculate theta.hat
  vec = eigen(cov(xdata))$values
  return(vec[1]/sum(vec))
}

library(bootstrap)
xdata = as.matrix(scor)
n = nrow(xdata)

theta.hat = theta2(xdata)

theta.jack=numeric(n)
for(i in 1:n){
  theta.jack[i] = theta2(xdata[-i,])         
}

(n-1) * (mean(theta.jack) - theta.hat)          #bias of Jacknife

sqrt((n-1)/n * sum((theta.jack - mean(theta.jack))^2))        #standard error of Jacknife

## ------------------------------------------------------------------------
library(lattice)
library(DAAG) 
attach(ironslag)

n = length(magnetic) #in DAAG ironslag 
e1 = e2 = e3 = e4 = numeric(n)

# for n-fold cross validation # fit models on leave-one-out samples 
for (k in 1:n) { 
  y = magnetic[-k] 
  x = chemical[-k]
  
  J1 = lm(y ~ x) 
  yhat1 = J1$coef[1] + J1$coef[2] * chemical[k] 
  e1[k] = magnetic[k] - yhat1

  J2 = lm(y ~ x + I(x^2)) 
  yhat2 = J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2 
  e2[k] = magnetic[k] - yhat2

  J3 = lm(log(y) ~ x) 
  logyhat3 = J3$coef[1] + J3$coef[2] * chemical[k] 
  yhat3 = exp(logyhat3) 
  e3[k] = magnetic[k] - yhat3
  
  # cubic polynomial model
  J4 = lm(y ~ poly(x,3,raw=T)) 
  yhat4 = J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] = magnetic[k] - yhat4

} 

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)) 

## ------------------------------------------------------------------------
y = magnetic
x = chemical
r = rep(0,4)

L1 = lm(y~x)
r[1] = summary(L1)$adj.r.squared # extract adjusted R^2

L2 = lm(y ~ x + I(x^2))
r[2] = summary(L2)$adj.r.squared

L3 = lm(log(y) ~ x) 
r[3] = summary(L3)$adj.r.squared

L4 = lm(y ~ poly(x,3,raw=T)) 
r[4] = summary(L4)$adj.r.squared

r

## ------------------------------------------------------------------------
maxout <- function(x, y) { 
  # copied from Page 177
  X <- x - mean(x) 
  Y <- y - mean(y) 
  outx <- sum(X > max(Y)) + sum(X < min(Y)) 
  outy <- sum(Y > max(X)) + sum(Y < min(X)) 
  return(max(c(outx, outy))) 
}

per_count5test = function(x,y,N=999, alpha = 0.05){
  # N is the number of the permutations
  # alpha is the significance level
  # return 0 if we accept that X and Y are from same distribution(accept H_0), and return 1 otherwise(reject H_0) 
  
  if(length(x)<=length(y)){
    m = length(x)
    n = length(y)
    z = c(x,y)
    theta = maxout(x,y)
  }
  else{     
    #if X has more samples than Y, we exchange them
    m = length(y)
    n = length(x)
    z = c(y,x)
    theta = maxout(y,x)
  }
  
  theta_hat = rep(0,N)
  
  for(i in 1:N){
    k = sample(1:(m+n), size = 2*m, replace = F)
    z1 = z[k]
    x1 = z1[1:m]
    y1 = z[(m+1):(2*m)]
    
    # Calculate the maxout for every permutation
    theta_hat[i] = maxout(x1, y1)
  }
  
  # calculate the empirical p value, see Page 217 for details
  p = (1+sum(theta_hat>=theta))/(N+1)
  
  # alpha = 0.05 is the desired significance level
  return(as.integer(p <= alpha))
}

## ------------------------------------------------------------------------
n1 = 20
n2 = 30
m = 1e4
D = rep(0,m)
for(i in 1:m){
  D[i] = per_count5test(rnorm(n1),rnorm(n2))
}
mean(D)

## ------------------------------------------------------------------------
n1 = 20
n2 = 30
m = 1e4
D = rep(0,m)
for(i in 1:m){
  D[i] = per_count5test(rnorm(n1),rchisq(n2, df=1))
}
# estimate the power
mean(D)

## ------------------------------------------------------------------------
library(Ball)
library(energy)


pow_comparison = function(n, N=1e3, alpha = 0.05){
  # return the estimated power of ball test and distance correlation test
  # n is the sample size
  # N is the replication times
  # alpha is the significance level
  
  p_bd = rep(0,n)
  p_dcor = rep(0,n)
  
  for(i in 1:N){
    X = matrix(rnorm(2*n), nrow = n, ncol = 2)
    e = matrix(rnorm(2*n), nrow = n, ncol = 2)
    Y1 = X/4 + e
    Y2 = X/4 * e
    p_bd[i] = bd.test(x = Y1, y = Y2, seed = i*n^2, R = 99)$p.value
    p_dcor[i] = dcor.test(x = Y1, y = Y2, R = 99)$p.value
  }
  return(c(mean(p_bd < alpha), mean(p_dcor < alpha)))
}


n = seq(5, 100, 5)
pow = matrix(0, nrow = length(n), ncol = 2)
for(i in 1:length(n)){
  pow[i,] = pow_comparison(n[i])
}
pow = cbind(n,pow)
colnames(pow) = c("sample size", "ball test", "dcor test")
knitr::kable(pow)

## ------------------------------------------------------------------------
plot(n, pow[,2], type="l", col = "blue", lty=5, ylim = c(0,1), ylab="power")
lines(n, pow[,3], type = "l", col = "red", lty =3)
legend("bottomright", c("ball","dcor"), fill = c("blue","red"))

## ------------------------------------------------------------------------
f = function(x){
  # return the density of standard Laplace distribution
  return(exp(-abs(x))/2)
}


Metropolis = function(m = 10000, sd = 1){
  # m is the number of replications
  # sd is the standard error of the normal distribution
  # k is the number of the rejected candidate points
  # x is the Markov Chain we obtain 
  
  x = numeric(m)
  
  #simulate from normal distribution
  x[1] = rnorm(1, mean = 0, sd = sd) 
  k = 0 
  u = runif(m)
  for (i in 2:m) { 
    xt = x[i-1] 
    y = rnorm(1, mean = xt, sd = sd)
    num = f(y)    # here it is the Metropolis sampler, see P253
    den = f(xt) 
    if (u[i] <= num/den){
      x[i] = y
    } 
    else { 
      x[i] = xt 
      k = k+1 #y is rejected 
    } 
  }
  return(list(x = x, k = k, sd = sd, m=m))
}

set.seed(1129)
#sd = 1
result = Metropolis(m = 10000, sd = 1)
1 - result$k / result$m # acceptance rate
index = 5000:5500 
y1 = result$x[index] 
plot(index, y1, type="l", main="sd=1", ylab="x") 

#sd = 10
result = Metropolis(m = 10000, sd = 10)
1 - result$k / result$m # acceptance rate
index = 5000:5500 
y1 = result$x[index] 
plot(index, y1, type="l", main="sd=10", ylab="x") 

#sd = 0.1
result = Metropolis(m = 10000, sd = 0.1)
1 - result$k / result$m # acceptance rate
index = 5000:5500 
y1 = result$x[index] 
plot(index, y1, type="l", main="sd=0.1", ylab="x") 

## ------------------------------------------------------------------------
x = rep(T,10)
for(i in 1:10){
  x[i] = log(exp(i))-exp(log(i))==0
}
x

## ------------------------------------------------------------------------
x = rep(T,10)
for(i in 1:10){
  x[i] = all.equal(log(exp(i))-exp(log(i)), 0)
}
x

## ------------------------------------------------------------------------
SS = function(a, k){
  # the function S_k in Exercise 11.4
  1 - pt(sqrt((a^2 *k)/(k+1-a^2)), df = k)
}


AA = function(k){
  # calculate the points A(k) in Exercise 11.4
  # k CANNOT be a vector
  
  # I found that the point is usually between 1 and 2... 
  uniroot(function(x) SS(x,k-1)-SS(x,k) , c(1,min(sqrt(k)-0.01,2)))$root    
}

## ------------------------------------------------------------------------
dt_new = function(x, df){
  # return the integrand in RHS of the equation in Exercise 11.5
  # return the same value as 2*dt(x, df)

  2*(1+x^2 /df)^(-(df+1)/2) /sqrt(pi*df) *exp(lgamma((df+1)/2)-lgamma(df/2))
}

III = function(a,k){
  # calculate the integral in the equation in Exercise 11.5
  
  c = sqrt(a^2 *k/(k+1-a^2))
  res = integrate(function(x) dt_new(x,k), lower = 0, upper = c)$value
  return(res)         
}

BB = function(k){
  # calculate the solution in Exercise 11.5
  # k CANNOT be a vector
  
  uniroot(function(x) III(x,k-1)-III(x,k) , c(1,min(sqrt(k)-0.01,2) ))$root
}

## ------------------------------------------------------------------------
K = c(c(4:25))
compare = rep(0,length(K))
for(i in 1:length(K)){
  compare[i] = all.equal(AA(K[i]),BB(K[i]))
}
max(compare) == 1 #return T iff A(k)=B(k) for all the k above

## ------------------------------------------------------------------------
ll = function(theta){
  p = theta[1]
  q = theta[2]
  -(28*log(2-p-2*q) +24*log(2-q-2*p) + 92*log(1-p-q) + 98*log(p) + 94*log(q))
}

optim(c(0.3,0.3), ll)

## ------------------------------------------------------------------------
l = function(theta, p_0,q_0){
  p = theta[1]
  q = theta[2]
  r = 1-p-q
  nAA = 28*p_0/(2-p_0-2*q_0)
  nBB = 24*q_0/(2-q_0-2*p_0)
  -((98+nAA)*log(p) + (94+nBB)*log(q) + (134-nAA-nBB)*log(r) )
}

# initialization
N = 10
Theta = matrix(0, ncol = 2, nrow = N+1)
Theta[1,] = c(0.3,0.3)

# EM algorithm
for(i in 1:N){
  Theta[i+1,] = optim(c(0.4,0.4), function(theta) l(theta, p_0 = Theta[i,1], q_0 = Theta[i,2]))$par
}

# result
Theta[N+1,]
ll(Theta[N+1,])

x = 1:N
y = rep(0,N)
for(i in 1:N){
  y[i] = -ll(Theta[i,])
}
plot(x,y, type="l", col="red")


## ----eval=FALSE----------------------------------------------------------
#  
#  #include <cmath>
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  //[[Rcpp::export]]
#  double f(double x) {
#    return exp(-abs(x));
#  }
#  
#  //[[Rcpp::export]]
#  NumericVector Metropolis_C (int N, double sd) {
#    NumericVector x(N);
#    NumericVector initial = rnorm(1,0,sd);
#    x[0] = initial[0];
#    NumericVector u = runif(N);
#    for (int i = 1; i < N;i++ ) {
#      NumericVector y = rnorm(1, x[i-1], sd);
#      if (u[i] <= (f(y[0]) / f(x[i-1]))){
#        x[i] = y[0];
#      }
#      else {
#        x[i] = x[i-1];
#      }
#    }
#    return(x);
#  }
#  

## ------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

# R
f_R = function(x){
    # return the density of standard Laplace distribution
    return(exp(-abs(x))/2)
}


Metropolis_R = function(m = 10000, sd = 1){
  # m is the number of replications
  # sd is the standard error of the normal distribution
  # k is the number of the rejected candidate points
  # x is the Markov Chain we obtain 
  
  x = numeric(m)
  
  #simulate from normal distribution
  x[1] = rnorm(1, mean = 0, sd = sd) 
  k = 0 
  u = runif(m)
  for (i in 2:m) { 
    xt = x[i-1] 
    y = rnorm(1, mean = xt, sd = sd)
    num = f_R(y)    # here it is the Metropolis sampler, see P253
    den = f_R(xt) 
    if (u[i] <= num/den){
      x[i] = y
    } 
    else { 
      x[i] = xt 
      k = k+1 #y is rejected 
    } 
  }
  return(list(x = x, k = k, sd = sd, m=m))
}

    dir_cpp = 'D:/Statistical_Computing/Homework/Homework_12/'
    sourceCpp(paste0(dir_cpp,"Metropolis_C.cpp"))
    N = 10000
    sd = 2
    (time = microbenchmark(rwR=Metropolis_R(N,sd),rwC=Metropolis_C(N,sd)))

## ------------------------------------------------------------------------

set.seed(12345)
rwR = Metropolis_R(N,sd)$x[-(1:1000)]
rwC = Metropolis_C(N,sd)[-(1:1000)]
qqplot(rwR,rwC)
abline(a=0,b=1,col='black')

