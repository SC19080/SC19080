% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Max_abs.R
\name{Max_abs}
\alias{Max_abs}
\title{Maximum of the Absolute Value of A Correlation Matrix}
\usage{
Max_abs(X)
}
\arguments{
\item{X}{the dimensions of the standard normal distribution}
}
\value{
the maximum, namely, a numeric element
}
\description{
This function provides calculates the maximum of the absolute value of a correlation matrix except the diagonal elements.
}
\examples{
\dontrun{
set.seed(12345)
data = cor(matrix(rnorm(50),nrow = 5, ncol = 10))
Max_abs(data)
}
}
