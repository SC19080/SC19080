#' @title Maximum of the Absolute Value of A Correlation Matrix
#' @description This function provides calculates the maximum of the absolute value of a correlation matrix except the diagonal elements.
#' @param X the dimensions of the standard normal distribution
#' @return the maximum, namely, a numeric element
#' @examples
#' \dontrun{
#' set.seed(12345)
#' data = cor(matrix(rnorm(50),nrow = 5, ncol = 10))
#' Max_abs(data)
#' }
#' @export
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

