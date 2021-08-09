##' The function estimates a trend from count survey data.
##' 
##' The function estimates smooth or loglinear population trends, or indexes from simple design count survey data. 
##' It is essentially a wrapper around a call to \code{\link[mgcv]{gam}}, processing its output using \code{\link[mgcv]{predict.gam}}
##' to produce a trend estimate.
##' Smooth trends for the temporal variable are set up by the term \code{s(var, k = k, fx = fx , bs = bs)}
##' where \code{var} is the first argument to \code{\link[poptrend]{trend}} in the formula.
##' For loglinear trends, the identity of \code{var} is used, and for index models a factor variable is constructed from \code{var}
##' with one level for each unique time point in \code{var} with corresponding observations in the data.
##'
##' Temporal random effects are set up by converting the temporal variable supplied to \code{\link{trend}} to a factor variable
##' and adding this factor variable to the data supplied to \code{\link[mgcv]{gam}}.
##' 
##' Bootstrap confidence intervals are computed by drawing normally distributed random variables with means equal to the
##' estimated coefficients and covariance matrix equal to the Bayesian posterior covariance matrix (see \link[mgcv]{vcov.gam}). 
##' 
##' @title Fit a smooth or linear trend to count survey data.
##' @param formula A trend formula. This is a GAM formula with an extra term \code{\link[poptrend]{trend}} describing the 
##'                time variable and properties of the trend. All terms except the trend term are treated as covariates. 
##'                Effects of temporal variation in these covariates are not included in the calculation of the trend. 
##' @param data A data frame containing response variables and covariates.
##' @param family The distributional family of the response. The family must use a log-link. Defaults to a negative binomial \code{\link[mgcv]{nb}}.
##' @param nGrid The number of grid points over which to compute the trend.
##'                  If the length of the argument is one, an equally spaced grid over the survey period of length nGrid is set up.
##'                  nGrid can also be a vector of length 3, in which case the first element is the number of grid points and the
##'                  second and third elements give, respectively, the start and endpoints of the grid.
##' @param nBoot The number of bootstrap samples to draw. 
##' @param engine   If 'gam', the default, model fitting is done via \code{\link[mgcv]{gam}}. If 'bam', model fitting is done via 
##'                 \code{\link[mgcv]{bam}}, which is less memory hungry and can be faster for large data sets.          
##' @param ... Further arguments passed to \code{\link[mgcv]{gam}}.
##' @return An object of class trend that inherits from class gam.
##' @examples
##' ## Simulate a data set with 15 sites and 25 years
##' data = simTrend(15, 25)
##' ## Fit a smooth trend with fixed site effects, random time effects,
##' ## and automatic selection of degrees of freedom
##' trSmo = ptrend(count ~ trend(year, tempRE = TRUE, type = "smooth") + site, data = data)
##' ## Check the model fit
##' gam.check(trSmo)
##' ## Plot the trend
##' plot(trSmo)
##' summary(trSmo)
##' ## Check the estimated percent change from year 8 to 25
##' change(trSmo, 8, 25)
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
##' 
##' ## trend objects inherit from 'gam' so that most methods available for gam objects 
##' ## can be directly used with trend objects. This can be used to extract further information.
##' 
##' summary.gam(trSmo)
##' 
##' AIC(trSmo)
##' 
##' plot(residuals(trSmo))
##' 
##' ## Compare smooth and loglinear models
##' anova(trSmo, trLin)
##' 
##' @export
##' @author Jonas Knape
ptrend = function(formula, data = list(), family = 'nb', nGrid = 500, nBoot = 500, engine = 'gam', ...) {
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
  if (!identical(family$link, "log")) {
    stop("Only log links allowed.")
  }
  ##
  if ('gamModel' %in% names(list(...))) {
    warning('Argument gamModel is deprecated. Will be ignored.')
  }
  ##
  tf = interpret.trendF(formula)
  timeVar = deparse(tf$tVar, width.cutoff = 500)
  if (!identical(timeVar, tf$predName)) {
    stop(paste("Trend variable", timeVar, "must be simple, expressions are not implemented."))
  }
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
    data[[timeVarFac]] = factor(tVarOut, exclude = as.character(tVarOut[!is.finite(tVarOut)]))
    contrasts(data[[timeVarFac]]) = contr.treatment(nlevels(data[[timeVarFac]]))
  }
  if (tf$type != "index") {
    if (length(nGrid) == 1) {
      minT = min(tPoints, na.rm = TRUE)
      maxT = max(tPoints, na.rm = TRUE)
      trendGrid = seq(minT - .01 * (maxT - minT), maxT + .01 * (maxT - minT), length.out = nGrid)
    } else if (length(nGrid) == 3) {
      if (nGrid[2] > nGrid[3])
        stop("nGrid[2] larger than nGrid[3]")
      trendGrid = seq(nGrid[2], nGrid[3], length.out = nGrid[1])
    } else 
      stop("trendGrid needs to be of length one or three")
  }
  
  # Fit model
  args = list(formula = tf$formula, data = data, family = family, ...)
  args$gamModel = NULL # Deprecated.
  if (engine == 'gam')
    gamFit = do.call(mgcv::gam, args)
  else 
    gamFit = do.call(mgcv::bam, args)
  # Construct data frame for plotting trend line
  if (tf$type != "index") {
    ix = sort.int(trendGrid, index.return = TRUE)$ix
    trendFrame = data.frame(trendGrid[ix])
    names(trendFrame) = tf$predName
    trendFrame[[deparse(tf$tVar)]] = eval(tf$tVar, trendFrame)
  } else {trendFrame = NULL}
  
  # Start constructing trend object
  out = c(gamFit, list(trendFormula = formula, 
             countVar = tf$response, timeVar = timeVar, timeVarFac = timeVarFac, tPredName = tf$predName, 
             timeRE = tf$tempRE, trendType = tf$type, k = tf$k, fx = tf$fx, ptrendCall = call))
  
  class(out) = class(gamFit)
  gamFit = NULL
  if (out$trendType != 'index') {
    out$timeVarLab = extractGamLabel(out, timeVar)
  }
  if (out$timeRE | out$trendType == 'index') {
    out$timeVarFacLab = extractGamLabel(out, timeVarFac)
  }
  # Construct data frame for annual indices
  tPoints = sort(tPoints)
  indexFrame = data.frame(tPoints)
  names(indexFrame) = out$tPredName
  indexFrame[[deparse(tf$tVar)]] = eval(tf$tVar, indexFrame)
  if (out$timeRE | out$trendType == "index") {
    tLevs =  levels(model.frame(out)[[timeVarFac]])
    indexFrame[[timeVarFac]] = factor(indexFrame[[timeVar]], levels = tLevs)
    
    if (out$trendType != "index") {
    trendFrame[[timeVarFac]] = factor(sapply(trendFrame[[timeVar]], 
                                             function(time, levs) {levs[which.min(abs(time - as.numeric(levs)))]}, 
                                              levs = tLevs), levels = tLevs)
#      trendFrame[[timeVarFac]] = indexFrame[[timeVarFac]][1]
    }
  }
  
  mf = get_all_vars(out, data)
  predFNames = attr(terms(out$pred.formula), "term.labels")
  mf = mf[which(colnames(mf) %in% setdiff(predFNames, c(out$response, out$tPredName, timeVarFac)))]
  
  # Set non-trend terms to constants for prediction. Note that this includes any offset.
  for (i in seq_len(ncol(mf))) {
    cname = colnames(mf)[i]
#    if(is.numeric(mf[[cname]])) {
#      indexFrame[[cname]] = mean(mf[[i]], na.rm = TRUE)
#    } else {
      indexFrame[[cname]] = mf[[i]][which(!is.na(mf[[i]]))[1]]
#    }
    if (out$trendType != "index") 
      trendFrame[[cname]] = indexFrame[[cname]][1]
  }

  rawT = computeRawTrend(out, indexFrame, estimate = TRUE, resid = TRUE, nBoot = 0)
  
  indexFrame = cbind(indexFrame, rawT)
  if (tf$type != "index") {
    rawT = computeRawTrend(out, trendFrame, estimate = TRUE, resid = FALSE, nBoot = 0)
    trendFrame = cbind(trendFrame, rawT)
  } 
  
  out$trendFrame = trendFrame
  out$indexFrame = indexFrame
  class(out) = c('trend', class(out)) 
  if (nBoot > 0) {
    out = hessBoot(out, nBoot)
  }  
  
  out = updateBaseline(out)
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
##'        If this is not appropriate, an alternative is to set tempRE to FALSE and manually add temporal random
##'        effects in the trend formula (using s(..., bs = "re")). Temporal random effects cannot be used with
##'        index estimation.
##' @param type The type of trend to be estimated. One of "smooth", "loglinear" or "index".
##' @param by Currently ignored.
##' @param k The dimension of the basis for the cubic regression spline of smooth trend fits.
##' @param fx If true, automatic selection of degrees of freedom are used for smooth trends.
##' @param bs The basis to use for smooth trends. Defaults to a cubic regression spline "cr".
##'           Currently available options are 'cr', 'cs' , 'tp', 'ts', and 'gp'. See \code{\link[mgcv]{smooth.terms}}.
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
trend = function(var, tempRE = FALSE, type = "smooth", by = NA, k = -1, fx = FALSE, bs = "cr") {
  type = match.arg(type, c("smooth", "loglinear", "index"))
  bs = match.arg(bs, c('cr', 'cs' , 'tp', 'ts', 'gp'))
  .tVarExt = "__RE"
  tVar = substitute(var)
  if (length(all.vars(tVar)) != 1)
    stop("Trend term must have one variable.")
  if (type == "smooth")
    gTrend = paste0("s(", deparse(tVar),", k = ", k, ", fx = ", fx, ", bs = \"",  bs, "\")")
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
  if (tempRE)
    gFac = paste0("s(", all.vars(tVar), .tVarExt, ", bs = \"re\")")
  else if (type == "index")
    gFac =paste0(all.vars(tVar), .tVarExt)
  else
    gFac = ""
  list(gTrend = gTrend, gFac = gFac, 
       tVar =  tVar, predName = all.vars(tVar), 
       k = k, fx = fx, tempRE = tempRE, type = type, tVarExt = .tVarExt)
}


#Extract labels for named terms
extractGamLabel = function(object, variable) {
  # TODO: Need to detect cases where 'variable' is hidden in compound terms, e.g. s(variable, foo)
  smoothTerms = sapply(object$smooth, function (x) c(paste(x[['term']], collapse = ""), x[['label']]))
  if (length(smoothTerms) == 0) 
    smoothTerms = NULL
  paramTerms = labels(terms(object$pterms))
  if (any(!(labels(terms(object)) %in% c(paramTerms, smoothTerms[1,]))))
    stop() 
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
##' @param object An object of class trend.
##' @param nBoot The number of bootstrap samples to draw.
##' @return A trend object with the bootstrapped trend estimates appended.
##' @export
##' @author Jonas Knape
hessBoot = function(object, nBoot = 500) {
  if (is.null(object$bootType)) {
    object$bootType = "hessian"
  }
  newData = object$indexFrame
  newDataT = NULL
  if (object$trendType != 'index') {
    newDataT = object$trendFrame
  }  
  newDataT$indexRaw = newDataT$resid = newData$indexRaw = newData$resid = NULL
  
  if (!is.null(newData) & !is.null(newDataT) & !identical(colnames(newData), colnames(newDataT)))
    stop('Something is wrong. Please report bug to package maintainer.')
  # Concatenate to use the same bootstrap parameters for indexFrame and trendFrame.
  boot = computeRawTrend(object, rbind(newData, newDataT), estimate = FALSE, resid = FALSE, nBoot)
  iInd = seq_len(nrow(newData))
  object[['boot']] = list(nBoot = boot$nBoot, index = boot$bootIndex[iInd, ])
  if (object$trendType != 'index') {
    tInd = max(iInd) + seq_len(nrow(newDataT))
    object$boot$trend = boot$bootIndex[tInd, ]
  }
  if (object$timeRE) 
    object$boot$residuals = boot$bootResid[iInd, ]

  ## TODO!!
  ## updateBaseline if hessBoot is called directly on previously fitted object. 
    
  object
}

computeRawTrend = function(object, newdata, estimate, resid, nBoot = 0) {
  is.null(newdata) # Just to give an error if newdata is missing. predict.gam swallows it.
  X = mgcv::predict.gam(object, newdata = newdata, type = "lpmatrix")
  if (object$trendType == "smooth") {
    tCol = grep(object$timeVarLab, colnames(X), fixed = TRUE)
    tResCol = integer(0)
  }
  if (object$trendType == "loglinear") {
    tCol = which(object$timeVarLab == colnames(X))
    tResCol = integer(0)
  }
  if (object$trendType == "index") {
    tResCol = integer(0)
    tCol = c(grep("(Intercept)", colnames(X), fixed = TRUE), grep(object$timeVarFac, colnames(X), fixed = TRUE))
  }
  if(object$timeRE) {
    tResCol = grep(object$timeVarFacLab, colnames(X), fixed = TRUE)
  } 
  if (length(intersect(tResCol, tCol)) > 0)
    stop()
  #X = X[, c(tCol, tResCol), drop = FALSE]
  cf = coef(object)
  
  out = list()
  
  if (estimate) {
    out$indexRaw = as.numeric(X[, tCol, drop = FALSE] %*% cf[tCol])
  }
  if (resid & object$timeRE) {
    out$resid = as.numeric(X[, tResCol, drop = FALSE] %*% cf[tResCol])
  }
  
  if (nBoot > 0) {
    # Subset to avoid unnecessary chol factorization of full variance matrix.
    X = X[ , c(tCol, tResCol), drop = FALSE]
    cf = cf[c(tCol, tResCol)]
    vc  = mgcv::vcov.gam(object)[c(tCol, tResCol), c(tCol, tResCol)]
    chvc = t(chol(vc)) # Use mgcv::mroot instead?
    cfn = cf + chvc %*% matrix(rnorm(length(cf) * nBoot), ncol = nBoot, nrow = length(cf))
    
    index = resid = NULL
    index = X[, seq_len(length(tCol)), drop = FALSE] %*% cfn[seq_len(length(tCol)), , drop = FALSE]
    if (object$timeRE) {
      resid = X[, length(tCol) + seq_len(length(tResCol)), drop = FALSE] %*% 
        cfn[length(tCol) + seq_len(length(tResCol)), , drop = FALSE]
    }
    
    out = c(out, list(nBoot = nBoot, bootIndex = index, bootResid = resid))
  }
  out
}




updateBaseline = function(object, baseline = NA, level = .95) {
  timeVar = object$timeVar
  if (length(baseline) == 1) { 
    if (is.na(baseline)) {
      baseline = min(object$indexFrame[[timeVar]])
    }
  }
  if (!all(baseline %in% object$indexFrame[[timeVar]])) { # TODO: Use numeric comparison instead of exact equality.
    stop('baseline can only contain time points with observations')
  }
  bootRefValue = NULL
  if (length(baseline) == 1) {
    trendRefValue = object$indexFrame$indexRaw[object$indexFrame[[timeVar]] == baseline]
    if(!is.null(object[['boot']][['index']])) {
      bootRefValue = object$boot$index[object$indexFrame[[timeVar]] == baseline,]
    }
  } else {
    trendRefValue = log(mean(exp(object$indexFrame$indexRaw[object$indexFrame[[timeVar]] %in% baseline])))
    if(!is.null(object[['boot']][['index']])) {
      bootRefValue = log(apply(exp(object$boot$index[object$indexFrame[[timeVar]] %in% baseline, ]), 2, mean))
    }
  }
  if (object$trendType != 'index') {
    object$trendFrame$index = exp(object$trendFrame$indexRaw - trendRefValue)
    object$trendFrame$d1 = fderiv(object$trendFrame[, 'index', drop = FALSE], object$trendFrame[[object$timeVar]], order = 1)
    object$trendFrame$d2 = fderiv(object$trendFrame[, 'index', drop = FALSE], object$trendFrame[[object$timeVar]], order = 2)
    attr(object$trendFrame, 'baseline') = baseline
    object$indexFrame$index = exp(object$indexFrame$indexRaw - trendRefValue)
    if (object$timeRE) {
      object$indexFrame$index.re = exp(object$indexFrame$indexRaw + object$indexFrame$resid - trendRefValue)
    }
  } else {
    object$indexFrame$index = exp(object$indexFrame$indexRaw - trendRefValue)
  }
  attr(object$indexFrame, 'baseline') = baseline
  
  object$trendRefValue = trendRefValue
  
  object$boot$bootRefValue = bootRefValue
  
  object = updateConfLevel(object, level = level)
  
  object
}

updateConfLevel = function(object, level) {
  if (!is.numeric(level) | level <= 0 | level >= 1)
    stop('level need to be a number between 0 and 1.')
  if (is.null(object[['boot']])) {
    stop('No uncertainty estimate available.')
  }
  if (is.null(object[['boot']][['bootRefValue']])) {
    stop('Baseline not set.')
  }
  bootcheck(object$boot$nBoot, level)
  probs = c(.5 - level / 2, .5 +  level / 2)
  bIT = object$boot$index - kronecker(rep(1,nrow(object$boot$index)), t(object$boot$bootRefValue))
  ciIT = t(apply(exp(bIT), 1, quantile, probs = probs))
  ciII = ciT = NULL
  if (object$timeRE) {
    ciII = t(apply(exp(bIT + object$boot$residuals), 1, quantile, probs = probs))
  }
  if (object$trendType != 'index') {
    bT = exp(object$boot$trend - kronecker(rep(1,nrow(object$boot$trend)), t(object$boot$bootRefValue)))
    ciT = apply(bT, 1, quantile, probs = probs)
    
    grad1 = apply(fderiv(bT, grid = object$trendFrame[[object$timeVar]], order = 1), 1, quantile, probs = probs)
    grad2 = apply(fderiv(bT, grid = object$trendFrame[[object$timeVar]], order = 2), 1, quantile, probs = probs)
    trendCI = data.frame(t(ciT), t(grad1), t(grad2))
    colnames(trendCI) = paste0(rep(c('ci.low.', 'ci.upp.'), 3), rep(c('index', 'd1', 'd2'), each = 2))
    attr(trendCI, 'level') = level
    object$trendCI = trendCI
  }
  if (object$timeRE) {
    ciII = t(apply(exp(bIT + object$boot$residuals), 1, quantile, probs = probs))
    indexCI = data.frame(ciIT,apply(exp(bIT), 1, sd), apply(bIT, 1, sd), ciII)
  } else {
    indexCI = data.frame(ciIT, apply(exp(bIT), 1, sd), apply(bIT, 1, sd))
  }
  if (object$trendType == 'index')
    colnames(indexCI) = c('ci.low.index', 'ci.upp.index', 'SE.index', 'SE.log.index')[1:ncol(indexCI)]
  else
    colnames(indexCI) = c('ci.low.index', 'ci.upp.index', 'SE.index', 'SE.log.index', 'ci.low.re', 'ci.upp.re')[1:ncol(indexCI)]
  attr(indexCI, 'level') = level
  object$indexCI = indexCI
  object

}

# Very simple check that there are 'enough' bootstrap samples.
# Requires at least 25 samples in the tails.
bootcheck =  function(nBoot, level) {
  if (nBoot * (1-level) >= 25) {
    return()
  } else {
    warning('Number of bootstrap samples (', nBoot, ') low for the given confidence level (', level, '). 
            Consider increasing the number of bootstraps to at least ' , ceiling(25/(1-level)), ', or use a lower confidence level.')
  }
}


