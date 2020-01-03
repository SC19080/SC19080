## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include=FALSE------------------------------------------------------
L1 = function(){
  
  itemp <- function(y, offset, parms, wt) {
    if (is.matrix(y) && ncol(y) > 1)
       stop("Matrix response not allowed")
    if (!missing(parms) && length(parms) > 0)
        warning("parameter argument ignored")
    if (length(offset)) y <- y - offset
    sfun <- function(yval, dev, wt, ylevel, digits ) {
		  paste("  mean=", format(signif(yval, digits)),
			", MSE=" , format(signif(dev/wt, digits)),
			sep = '')
    }
    environment(sfun) <- .GlobalEnv
    list(y = c(y), parms = NULL, numresp = 1, numy = 1, summary = sfun)
  }

  etemp <- function(y, wt, parms) {
    wmean <- median(y)
    rss <- sum(abs(y-wmean))
    list(label = wmean, deviance = rss)
  }

  stemp <- function(y, wt, x, parms, continuous){
    
    n <- length(y)
    
    if (continuous) {
        
        goodness = rep(0,n-1)
        direction = rep(0,n-1)
        
        for(i in 1:(n-1)){
            lmedian = median(y[1:i])
            rmedian = median(y[(i+1):n])
            goodness[i] = 1/log(sum(abs(y[1:i]-lmedian)) + sum(abs(y[(i+1):n]-rmedian)) )
            direction[i] = sign(lmedian - rmedian)
        }
        
        list(goodness = goodness, direction = direction)
        
    } 
    else {
        
        ux = sort(unique(x))
        ymean  <- tapply(y,x,mean)
        ymean_data = as.data.frame(ymean)
        ord <- order(ymean)
        m <- length(ord)
        
        goodness = rep(0,m-1)
        direction = rep(0,m-1)
        
        INDEX = rep(0,n)

        for(i in 1:n){
            INDEX[i] = which(rownames(ymean_data) == x[i])
        }
        

        for(i in 1:m-1){
            index_left = which(INDEX %in% ord[1:i])
            y_left = y[index_left]
            y_right = y[-index_left]
            lmedian = median(y_left)
            rmedian = median(y_right)
            goodness[i] = 1/log(sum(abs(y_left-lmedian)) + sum(abs(y_right-rmedian)) )
        }
        
        
        list(goodness= goodness, direction = ux[ord])
    }
  }

  list(eval = etemp, split = stemp, init = itemp)
}

library(kedd)

Max_abs = function(X){
    p = ncol(X)
    Max=0
    for(i in 1:p){
      for(j in 2:p){
        if(j>i){
          if(abs(X[i,j])>Max){
            Max = abs(X[i,j])
          }
        }
      }
    }
    return(Max)
}

Max_abs_cor= function(n = 30, p = 200, replicate = 500){

  n = n
  p = p
  m = replicate
  x = rep(0,m)

  for(i in 1:m){
    X = matrix(rnorm(p*n), nrow = n, ncol = p)
    x[i] = SC19080::Max_abs(cor(X))
  }

  kedd::dkde(x, deriv.order = 0)
}

## ------------------------------------------------------------------------
library(rpart)
library(rpart.plot)
mystate <- data.frame(state.x77, region=state.region)
names(mystate) <- casefold(names(mystate))
fit1 <- rpart(murder ~ population + illiteracy + income + life.exp +hs.grad + frost + region, data = mystate, method = L1(), minsplit = 10)
rpart.plot(fit1)

## ------------------------------------------------------------------------
data = matrix(rnorm(5*4), nrow = 5, ncol = 4)
cor(data)
Max_abs(cor(data))

## ------------------------------------------------------------------------
set.seed(12345)
result = Max_abs_cor(n = 30, p = 200, replicate = 500)
plot(result)

