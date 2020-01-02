#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

double f(double x) {
  return exp(-abs(x));
}

//' @title A Metropolis sampler using Rcpp
//' @description A Metropolis sampler using Rcpp used in homework 12
//' @param N the number of samples
//' @param sd the standard error of the normal distribution
//' @return a random sample of size \code{n}
//' @export
//[[Rcpp::export]]
NumericVector Metropolis_C (int N, double sd) {
  NumericVector x(N);
  NumericVector initial = rnorm(1,0,sd);
  x[0] = initial[0];
  NumericVector u = runif(N);
  for (int i = 1; i < N;i++ ) {
    NumericVector y = rnorm(1, x[i-1], sd);
    if (u[i] <= (f(y[0]) / f(x[i-1]))){
      x[i] = y[0];
    }
    else {
      x[i] = x[i-1];
    }
  }
  return(x);
}
