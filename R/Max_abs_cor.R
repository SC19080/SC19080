#' @title Estimation of Distributions of Maximum Absolute Sample Correlation Coefficient
#' @description This function provides an estimation of the maximum absolute sample correlation coefficient of standard normal distribution in high dimensions, which is done by generating samples and an Gaussian kernel estimation.
#' @param n the number of samples
#' @param p the dimensions of the standard normal distribution
#' @param replicate the number of replication times, i.e., the number of samples for the kernel density estimation
#' @return the same object as the function dkde does
#' @examples
#' \dontrun{
#' set.seed(12345)
#' result = Max_abs_cor(n = 30, p = 200, replicate = 500)
#' plot(result)
#' }
#' @export
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
