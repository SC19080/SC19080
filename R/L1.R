#' @title L1 splitting function
#' @description This function provides the L1 splitting function used in rpart. For more details, please refer to the description of the package rpart and its vignettes.
#' @return Returns a list of functions used in the argument method in rpart
#' @examples
#' \dontrun{
#' library(rpart)
#' library(rpart.plot)
#' mystate <- data.frame(state.x77, region=state.region)
#' names(mystate) <- casefold(names(mystate))
#' fit1 <- rpart(murder ~ population + illiteracy + income + life.exp +hs.grad + frost + region, data = mystate, method = L1())
#' rpart.plot(fit1, type =3)
#' }
#' @export
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

  return(list(eval = etemp, split = stemp, init = itemp))
}
