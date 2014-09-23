\name{snht-package}
\alias{snht-package}
\docType{package}
\title{
Robust and Non-Robust Standard Normal Homogeneity Test
}
\description{
Computes test statistics for the SNHT and robust SNHT test.  Additionally, users may supply a custom function for estimating the mean and standard deviation, and this function will be used for computing the test statistic.
}
\details{
\tabular{ll}{
Package: \tab snht\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-23\cr
License: \tab What license is it under?\cr
}
The main function is snht, which then calls the other functions in this package.  However, users may also wish to call robustSNHT which allows for a custom estimator function.
}
\author{
Josh Browning

Maintainer: Josh Browning <jbrownin@mines.edu>
}
\references{
L. Haimberger. Homogenization of radiosonde temperature time series using innovation
statistics. Journal of Climate, 20(7): 1377-1403, 2007.
}
\keyword{ package }
