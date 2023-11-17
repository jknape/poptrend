##' The function estimates a trend from count survey data.
##' 
##' The function estimates smooth or loglinear population trends, or indexes from simple design count survey data. 
##' It is essentially a wrapper around a call to \code{\link[mgcv]{gam}}, processing its output using \code{\link[mgcv]{predict.gam}}
##' to produce a trend estimate.
##' For smooth trends, cubic regression splines for the temporal variable are set up by the term \code{s(var, k = k, fx = fx , bs = "cr")}
##' where \code{var} is the first argument to \code{\link[poptrend]{trend}} in the formula.
##' For loglinear trends, the identity of \code{var} is used, and for index models a factor variable is constructed from \code{var}.
##'
##' Temporal random effects are set up by converting the temporal variable supplied to \code{\link{trend}} to a factor variable
##' and adding this factor variable to the data supplied to \code{\link[mgcv]{gam}}.
##' 
##' Bootstrap confidence intervals are computed by drawing normally distributed random variable with means equal to the
##' estimated coefficients and covariance matrix equal to the Bayesian posterior covariance matrix (see \link[mgcv]{vcov.gam}). 
##' 
##' @title Fit a smooth or linear trend to count survey data.
##' @param formula A trend formula. This is a GAM formula with an extra term \code{\link[poptrend]{trend}} describing the 
##'                time variable and properties of the trend. All terms except the trend term are treated as covariates. 
##'                Effect of temporal variation in these covariates are not included in the calculation of the trend. 
##' @param data A data frame containing response variables and covariates.
##' @param family The distributional family of the response. The family most use a log-link, it defaults to a quasi-Poisson.
##' @param nGrid The number of grid points over which to compute the trend.
##'                  If the length of the argument is one, an equally spaced grid over the survey period of length nGrid is set up.
##'                  nGrid can also be a vector of length 3, in which case the first element is the number of grid points and the
##'                  second and third elements give, respectively, the start and endpoints of the grid.
##' @param nBoot The number of bootstrap samples to draw. 
##' @param bootType Only one method, "hessian", currently implemented. Type "hessian", draws bootstrap samples using the Bayesian
##'                 covariance matrix of the parameters (see \code{\link[mgcv]{vcov.gam}}).
##' @param gamModel If true, the fit of the underlying gam model is saved in the output. May be set to FALSE to save memory,
##'                 but with the side effect that the fit of the gam model cannot be checked.
##' @param engine   If 'gam', the default, model fitting is done via \code{\link[mgcv]{gam}}. If 'bam', model fitting is done via 
##'                 \code{\link[mgcv]{bam}}, which is less memory hungry and can be faster for large data sets.          
##' @param ... Further arguments passed to \code{\link[mgcv]{gam}}.
##' @return An object of class trend.
##' @examples
##' ## Simulate a data set with 15 sites and 25 years
##' data = simTrend(15, 25)
##' ## Fit a smooth trend with fixed site effects, random time effects,
##' ## and automatic selection of degrees of freedom
##' trFit = ptrend(count ~ trend(year, tempRE = TRUE, type = "smooth") + site, data = data)
##' ## Check the model fit
##' checkFit(trFit)
##' ## Plot the trend
##' plot(trFit)
##' summary(trFit)
##' ## Check the estimated percent change from year 8 to 25
##' change(trFit, 8, 25)
##' 
##' ## Fit a loglinear trend model with random site effects and random time effects 
##' ## to the same data set.
##' trLin = ptrend(count ~ trend(year, tempRE = TRUE, type = "loglinear") +
##'                  s(site, bs = "re"), data = data)
##' plot(trLin)
##' summary(trLin)
##' 
##' ## Fit an index model with fixed site effects and an (unrelated) continous covariate 
##' ## as a smooth effect.
##' # Simulate mock covariate unrelated to data.
##' cov = rnorm(nrow(data))
##' trInd = ptrend(count ~ trend(year, type = "index") + site + s(cov), data = data)
##' plot(trInd)
##' summary(trInd)
##' @export
##' @author Jonas Knape
ptrend = function(formula, data = list(), family = quasipoisson(), nGrid = 500, nBoot = 500, bootType = "hessian", gamModel = TRUE, engine = 'gam', ...) {
  engine = match.arg(engine, c('gam', 'bam'))
  call = match.call()
  if (!inherits(formula, "formula"))
    stop("Argument formula needs to be an object of class formula.")
  if (is.character(family)) 
    family <- eval(parse(text = family))
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) 
    stop("family not recognized")
  if(family$link != "log") {
    stop("Only log links allowed.")
  }
  tf = interpret.trendF(formula)
  timeVar = deparse(tf$tVar, width.cutoff = 500)
  if (!identical(timeVar, tf$predName)) {
    stop(paste("Trend variable", timeVar, "must be simple, expressions are not implemented."))
  }
#  tPoints = unique(eval(as.name(tf$predName), data, environment(formula)))
  tPoints = unique(eval(tf$tVar, data, environment(formula)))
  tPoints = tPoints[!is.na(tPoints)]
  tVarOut = eval(tf$tVar, data, environment(formula))
  if (tf$type %in% c("smooth", "loglinear")) {
    if(!is.numeric(tVarOut) | !is.numeric(tPoints))
      stop(paste("Trend variable needs to be numeric."))
  }  
  timeVarFac = NULL
  if(tf$tempRE | tf$type == "index") {
    timeVarFac = paste0(tf$predName, tf$tVarExt)
    data[[timeVarFac]] = factor(tVarOut)
    contrasts(data[[timeVarFac]]) = contr.treatment(nlevels(data[[timeVarFac]]))
  }
  if (tf$type == "index") {
    trendGrid =tPoints
  } else {
    if (length(nGrid) == 1) {
      minT = min(tPoints, na.rm = TRUE)
      maxT = max(tPoints, na.rm = TRUE)
      trendGrid = c(tPoints, seq(minT - .01 * (maxT - minT), maxT + .01 * (maxT - minT), length.out = nGrid))
    } else if (length(nGrid) == 3) {
      if (nGrid[2] > nGrid[3])
        stop("nGrid[2] larger than nGrid[3]")
      trendGrid = c(tPoints, seq(nGrid[2], nGrid[3], length.out = nGrid[1]))
    } else 
      stop("trendGrid needs to be of length one or three")
  }
  ix = sort.int(trendGrid, index.return = TRUE)$ix
  trendFrame = data.frame(trendGrid[ix])
  names(trendFrame) = tf$predName
  trendFrame$isGridP = c(rep(FALSE, length(tPoints)), rep(TRUE, nGrid[1]))[ix]
  trendFrame[[deparse(tf$tVar)]] = eval(tf$tVar, trendFrame)
  
  if (engine == 'gam') {
    gamFit = mgcv::gam(formula = tf$formula, data = data, family = family, ...) 
  } else { 
    gamFit = mgcv::bam(formula = tf$formula, data = data, family = family, ...) 
  }
  if (tf$tempRE | tf$type == "index") {
    tLevs =  levels(model.frame(gamFit)[[timeVarFac]])
    trendFrame[[timeVarFac]] = factor(sapply(trendFrame[[timeVar]], 
                                             function(time, levs) {levs[which.min(abs(time - as.numeric(levs)))]}, 
                                             levs = tLevs), levels = tLevs)
  }
  mf = get_all_vars(gamFit, data)
  predFNames = attr(terms(gamFit$pred.formula), "term.labels")
  mf = mf[which(colnames(mf) %in% setdiff(predFNames, c(tf$response, tf$predName, timeVarFac)))]
  
  for (i in seqP(1, ncol(mf))) {
    cname = colnames(mf)[i]
    if(is.numeric(mf[[cname]])) {
      trendFrame[[cname]] = mean(mf[[i]], na.rm = TRUE) # dfhead[, i][1]
    } else {
      #browser()
      trendFrame[[cname]] = mf[[i]][which(!is.na(mf[[i]]))[1]]
    }
  }
  pred = predict(gamFit, trendFrame, block.size = nrow(trendFrame), type = "lpmatrix") # What if contrasts set differently
  if (tf$type == "smooth") {
    tCol = grep(extractGamLabel(timeVar, formula(gamFit)), colnames(pred), fixed = TRUE)
  }
  if (tf$type == "loglinear") {
    tCol = which(colnames(pred) == timeVar)
  }
  if (tf$type == "index") {
    tCol = c(grep("(Intercept)", colnames(pred), fixed = TRUE), grep(timeVarFac, colnames(pred), fixed = TRUE))
  }
  if (tf$type != "index") {
    trendFrame$trend = exp(pred[, tCol, drop = FALSE] %*% coef(gamFit)[tCol, drop = FALSE]) # Assumes log-link!
    trendFrame$trendResid = rep(NA, nrow(trendFrame))
  } else {
    trendFrame$trend = NA
    trendFrame$trendResid = exp(pred[, tCol, drop = FALSE] %*% coef(gamFit)[tCol, drop = FALSE])
  }
  if(tf$tempRE) {
    tCol = grep(extractGamLabel(timeVarFac, formula(gamFit)), colnames(pred), fixed = TRUE)
    trendFrame$trendResid = exp(pred[, tCol] %*% coef(gamFit)[tCol])
  }
  out = list(trendFrame = trendFrame,
             formula = formula, family = family, gam = gamFit,
             countVar = tf$response, timeVar = timeVar, timeVarFac = timeVarFac, tPredName = tf$predName, 
             timeRE = tf$tempRE, trendType = tf$type, k = tf$k, fx = tf$fx, bootTrend = NULL, bootResid = NULL, bootType = NULL,
             call = call)
  class(out) = "trend"
  if (nBoot > 0) {
    if (bootType == "hessian") 
      out = hessBootstrap(out, nBoot)
    if (bootType == "semi-parametric") 
      stop("semi-parametric bootstrap method not implemented.")
    if (bootType == "fewster")
      stop("Fewster bootstrap method not implemented.")
  }  
  if (!gamModel)
    out$gam = NULL
  out
}

interpret.trendF = function(formula) {
  tf = terms(formula, specials = "trend")
  response  = attr(tf, "variables")[[attr(tf, "response") + 1]]
  tind1 = attr(tf, "specials")$trend
  if (is.null(tind1))
    stop("No trend term found in formula.")
  if (length(tind1) > 1)
    stop("Only one trend term allowed in formula.")
  tCall = attr(tf, "variables")[[tind1 + 1]]
  #tind2 = grep(deparse(tCall, width.cutoff = 500), attr(tf, "term.labels"), fixed = TRUE)
  #tf = drop.terms(tf, dropx = tind2, keep.response = TRUE)
  trval = eval(tCall)
  if (trval$type != "index") {
    newF = gsub(deparse(tCall, width.cutoff = 500), trval$gTrend, x = deparse(tf, width.cutoff = 500), fixed = TRUE)
  } else {
    newF = gsub(deparse(tCall, width.cutoff = 500), trval$gFac, x = deparse(tf, width.cutoff = 500), fixed = TRUE)
  }
  if (trval$tempRE)
    newF = paste0(newF, " + ", trval$gFac)
  c(list(formula = as.formula(newF, env = environment(formula))), response = response,  trval)
  #update.formula(tf, newF)
}

##' The function is used to set up the trend component used in ptrend formulas.
##'
##' The function extracts information about the trend component of a formula supplied to ptrend. 
##' It returns a list containing variable names, information, and \code{\link[mgcv]{s}} components as strings used in subsequent calls to gam.
##' @title Define a trend component.
##' @param var A numeric time variable over which a trend or index will be computed.
##' @param tempRE If TRUE, this will set up random time effects. The random effects will be constructed by converting the
##'        var argument to a factor. Note that this yields a random effect level for each unique value in var.
##'        If this is not appropriate, an alternative is to set tempRE to false and manually add temporal random
##'        effects in the trend formula (using s(..., bs = "re")). Temporal random effects cannot be used with
##'        index estimation.
##' @param type The type of trend to be estimated. One of "smooth", "loglinear" or "index".
##' @param by Currently ignored.
##' @param k The dimension of the basis for the cubic regression spline of smooth trend fits.
##' @param fx If true, automatic selection of degrees of freedom are used for smooth trends.
##' @return A list containing information to set up the trend.
##' @export
##' @author Jonas Knape
##' @examples 
##' ## Simulate a data set with 15 sites and 25 years
##' data = simTrend(15, 25)
##' ## Fit a smooth trend with fixed site effects, but no random time effects,
##' ## and fixed degrees of freedom
##' trFit = ptrend(count ~ trend(year, tempRE = FALSE, k =  8, fx = FALSE, type = "smooth") +
##'                  site, data = data)
##' plot(trFit)
trend = function(var, tempRE = FALSE, type = "smooth", by = NA, k = -1, fx = FALSE) {
  #browser()
  type = match.arg(type, c("smooth", "loglinear", "index"))
  .tVarExt = "__Fac"
  tVar = substitute(var)
  if (length(all.vars(tVar)) != 1)
    stop("Trend term must have one variable.")
  if (type == "smooth")
    gTrend = paste0("s(", deparse(tVar),", k = ", k, ", fx = ", fx, ", bs = \"cr\")")
  if (type == "loglinear") {
#    tVar = substitute(I(tVar)) 
    gTrend = deparse(tVar)
  }
  if (type == "index") {
    if (tempRE)
      stop("Random effects cannot be added to index model.")
    #tVar = substitute(factor(tVar)) 
    gTrend = ""
  }
  if (tempRE) {
    gFac = paste0("s(", all.vars(tVar), .tVarExt, ", bs = \"re\")")
  } else if (type == "index") {
    gFac =paste0(all.vars(tVar), .tVarExt)
  } else {
    gFac = ""
  }
  list(gTrend = gTrend, gFac = gFac, 
       tVar =  tVar, predName = all.vars(tVar), 
       k = k, fx = fx, tempRE = tempRE, type = type, tVarExt = .tVarExt)
}


#Extract labels for named terms
extractGamLabel = function(variable, formula) {
  # TODO: Need to detect cases where 'variable' is hidden in compound terms, e.g. s(variable, foo)
  smoothTerms = sapply(mgcv::interpret.gam(formula)$smooth, function (x) c(paste(x$term, collapse = ""), x$label))
  paramTerms = labels(terms(mgcv::interpret.gam(formula)$pf))
  tlInd = which(variable == smoothTerms[1, ])
  if (length(tlInd) == 1) {
    if(variable %in% paramTerms) warning(paste("variable", variable, "occurs as both smoothed and parametric term, smooth term used."))
    return(smoothTerms[2, tlInd])
  } 
  if (length(tlInd) > 1) {
    stop(paste("variable", variable, "occurs in multiple smooths."))
  }
  if (length(tlInd) == 0) {
    tlInd = which(variable == paramTerms)
  }
  if (length(tlInd) != 1)
    stop(paste("variable", variable, "not found isolated in gam formula."))
  paramTerms[tlInd]  
}

convertTimeToFactor = function(times, fac) {
  sapply(times, 
         function(time, levs) {levs[which.min(abs(time - as.numeric(levs)))]}, 
         levs = levels(fac))
}


##' Draws bootstrap samples using the estimated variance matrix of the fitted gam model.
##' 
##' This function is used by \link{ptrend} and would typically not be called directly.
##' Bootstrap samples are drawn using the estimated coefficients and covariance matrix \link[mgcv]{vcov.gam} 
##' of the fitted gam model. The default values of \link[mgcv]{vcov.gam} which gives the Bayesian posterior
##' covariance matrix.
##' 
##' Bootstrapped samples computed in this way do not account for any uncertainty in the selection of degrees
##' of freedom.
##' @title Compute bootstrap confidence intervals based on sampling from the variance-covariance matrix.
##' @param trend An object of class trend.
##' @param nBoot The number of bootstrap samples to draw.
##' @return A trend object with the bootstrapped trend estimates appended.
##' @export
##' @author Jonas Knape
hessBootstrap = function(trend, nBoot = 500) {
  if (is.null(trend$bootType)) {
    trend$bootType = "hessian"
  } else if (trend$bootType != "hessian") {
    stop("Can't mix different bootstrap methods.")
  }
  X = predict(trend$gam, newdata = trend$trendFrame, type = "lpmatrix")
  if (trend$trendType == "smooth") {
    tCol = grep(extractGamLabel(trend$timeVar, formula(trend$gam)), colnames(X), fixed = TRUE)
    tFacCol = integer(0)
  }
  if (trend$trendType == "loglinear") {
    tCol = which(colnames(X) == trend$timeVar)
    tFacCol = integer(0)
  }
  if (trend$trendType == "index") {
    tCol = integer(0)
    tFacCol = c(grep("(Intercept)", colnames(X), fixed = TRUE), grep(trend$timeVarFac, colnames(X), fixed = TRUE))
  }
  if(trend$timeRE) {
    tFacCol = grep(extractGamLabel(trend$timeVarFac, formula(trend$gam)), colnames(X), fixed = TRUE)
  } 
  X = X[, c(tCol, tFacCol), drop = FALSE]
  cf = coef(trend$gam)[c(tCol, tFacCol)]
  vc  = vcov(trend$gam)[c(tCol, tFacCol), c(tCol, tFacCol)]
  chvc = t(chol(vc))
  bdt = NULL
  bdtFac = NULL
  cfn = cf + chvc %*% matrix(rnorm(length(cf) * nBoot), ncol = nBoot, nrow = length(cf))
  if (trend$trendType != "index") {
    bdt = cbind(bdt, exp(X[, seqP(1,length(tCol)), drop = FALSE] %*% cfn[seqP(1, length(tCol)), , drop = FALSE]))
  }  else { # Standardize each simulation to avoid numerical problems. Could also be done for non-index trends. TODO.
    XC = X[, seqP(length(tCol) + 1, length(tCol) + length(tFacCol)), drop = FALSE] %*% 
      cfn[seqP(length(tCol) + 1, length(tCol) + length(tFacCol)), , drop = FALSE]
    bdtFac = cbind(bdtFac, exp(XC-kronecker(rep(1,nrow(XC)), t(colSums(XC)/nrow(XC)))))
  }
  if (trend$timeRE) {
    bdtFac = cbind(bdtFac, exp(X[, seqP(length(tCol) + 1, length(tCol) + length(tFacCol)), drop = FALSE] %*% 
                                 cfn[seqP(length(tCol) + 1, length(tCol) + length(tFacCol)), , drop = FALSE]))
  }
  trend$bootTrend = cbind(trend$bootTrend, bdt)
  trend$bootResid = cbind(trend$bootResid, bdtFac)
  trend
}

seqP = function(from, to) {
  if (to >= from) return(from:to)
  return(integer(0))
}



