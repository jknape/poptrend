#' Analyze population trends from survey count data.
#'
#' The package provides functions for fitting and analysing trend models of data obtained 
#' from population count surveys.
#' 
#' The package provides functions for estimating smooth trends with generalized additive mixed models, 
#' as well as linear trends and population indices.
#' It is intended as a simple interface to basic trend estimation, allowing estimation of
#' trends accounting for effects of covariates in the form of both smooth terms and random effects.
#' The model fitting engine is the function \code{\link[mgcv]{gam}} of package \link[mgcv]{mgcv}.
#' Background for the package is given in Knape (2016).
#' 
#' @name poptrend
#' @docType package
#' @references Knape, J. 2016. Decomposing trends in Swedish bird populations using generalized additive mixed models. 
#'             Journal of Applied Ecology, 53:1852-1861. DOI:10.1111/1365-2664.12720.
#' @import stats mgcv
#' @importFrom graphics grid lines par plot points polygon segments
#' @importFrom grDevices adjustcolor
NULL
