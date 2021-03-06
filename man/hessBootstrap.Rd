% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trendFit.R
\name{hessBootstrap}
\alias{hessBootstrap}
\title{Compute bootstrap confidence intervals based on sampling from the variance-covariance matrix.}
\usage{
hessBootstrap(trend, nBoot = 500)
}
\arguments{
\item{trend}{An object of class trend.}

\item{nBoot}{The number of bootstrap samples to draw.}
}
\value{
A trend object with the bootstrapped trend estimates appended.
}
\description{
Draws bootstrap samples using the estimated variance matrix of the fitted gam model.
}
\details{
This function is used by \link{ptrend} and would typically not be called directly.
Bootstrap samples are drawn using the estimated coefficients and covariance matrix \link[mgcv]{vcov.gam} 
of the fitted gam model. The default values of \link[mgcv]{vcov.gam} which gives the Bayesian posterior
covariance matrix.

Bootstrapped samples computed in this way do not account for any uncertainty in the selection of degrees
of freedom.
}
\author{
Jonas Knape
}
