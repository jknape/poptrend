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
##' @return An object of class trend.
##' @examples
##' ## Simulate a data set with 15 sites and 25 years
##' data = simTrend(15, 25)
##' ## Fit a smooth trend with fixed site effects, random time effects,
##' ## and automatic selection of degrees of freedom
##' trSmo = ptrend(count ~ trend(year, tempRE = TRUE, type = "smooth") + site, data = data)
##' ## Check the model fit
##' checkFit(trSmo)
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
  
  # Construct data frame for annual indices
  tPoints = sort(tPoints)
  indexFrame = data.frame(tPoints)
  names(indexFrame) = tf$predName
  indexFrame[[deparse(tf$tVar)]] = eval(tf$tVar, indexFrame)
  
  if (tf$tempRE | tf$type == "index") {
    tLevs =  levels(model.frame(gamFit)[[timeVarFac]])
    indexFrame[[timeVarFac]] = factor(indexFrame[[timeVar]], levels = tLevs)
    
    if (tf$type != "index") {
    trendFrame[[timeVarFac]] = factor(sapply(trendFrame[[timeVar]], 
                                             function(time, levs) {levs[which.min(abs(time - as.numeric(levs)))]}, 
                                              levs = tLevs), levels = tLevs)
#      trendFrame[[timeVarFac]] = indexFrame[[timeVarFac]][1]
    }
  }
  
  mf = get_all_vars(gamFit, data)
  predFNames = attr(terms(gamFit$pred.formula), "term.labels")
  mf = mf[which(colnames(mf) %in% setdiff(predFNames, c(tf$response, tf$predName, timeVarFac)))]
  
#  for (i in seqP(1, ncol(mf))) {
  for (i in seq_len(ncol(mf))) {
    cname = colnames(mf)[i]
#    if(is.numeric(mf[[cname]])) {
#      indexFrame[[cname]] = mean(mf[[i]], na.rm = TRUE)
#    } else {
      indexFrame[[cname]] = mf[[i]][which(!is.na(mf[[i]]))[1]]
#    }
    if (tf$type != "index") 
      trendFrame[[cname]] = indexFrame[[cname]][1]
  }
  predIndex = predict(gamFit, indexFrame, block.size = nrow(indexFrame), type = "lpmatrix")
  if (tf$type != "index") {
    predTrend = predict(gamFit, trendFrame, block.size = nrow(trendFrame), type = "lpmatrix") # What if contrasts set differently
    stopifnot(identical(colnames(predTrend), colnames(predIndex))) # DEBUG
  }
  if (tf$type == "smooth") {
    tCol = grep(extractGamLabel(timeVar, formula(gamFit)), colnames(predIndex), fixed = TRUE)
  }
  if (tf$type == "loglinear") {
    tCol = which(colnames(predIndex) == timeVar)
  }
  if (tf$type == "index") {
    stopifnot(any(grepl("(Intercept)", colnames(predIndex), fixed = TRUE))) # DEBUG
    tCol = c(grep("(Intercept)", colnames(predIndex), fixed = TRUE), grep(timeVarFac, colnames(predIndex), fixed = TRUE))
  }
  indexFrame$indexRaw = as.numeric(predIndex[, tCol, drop = FALSE] %*% coef(gamFit)[tCol, drop = FALSE])
  if (tf$type != "index") {
    trendFrame$indexRaw = as.numeric(predTrend[, tCol, drop = FALSE] %*% coef(gamFit)[tCol, drop = FALSE])
  } 
  if(tf$tempRE) {
    tCol = grep(extractGamLabel(timeVarFac, formula(gamFit)), colnames(predIndex), fixed = TRUE)
    indexFrame$resid = as.numeric(predIndex[, tCol] %*% coef(gamFit)[tCol])
  }
  
  
#  out = list(trendFrame = trendFrame, indexFrame = indexFrame,
#             formula = formula, family = family, gam = gamFit,
#             countVar = tf$response, timeVar = timeVar, timeVarFac = timeVarFac, tPredName = tf$predName, 
#             timeRE = tf$tempRE, trendType = tf$type, k = tf$k, fx = tf$fx, bootI = NULL, bootIRes = NULL, bootT = NULL,
#             call = call)
#  class(out) = "trend"

  ## TEST
  out = list(trendFrame = trendFrame, indexFrame = indexFrame,
                          formula = formula, family = family,
                          countVar = tf$response, timeVar = timeVar, timeVarFac = timeVarFac, tPredName = tf$predName, 
                          timeRE = tf$tempRE, trendType = tf$type, k = tf$k, fx = tf$fx,
                          call = call)
  
  
  out = c(gamFit, out)
  class(out) = c('trend', class(gamFit)) 
  
  ## END TEST
  
  if (nBoot > 0) {
    out = hessBootstrap(out, nBoot)
  }  
  out = setBaseline(out)
  
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
  bs = match.arg(bs, c('cr', 'cs' , 'tp', 'ts', 'gp', 'ds'))
  .tVarExt = "__Fac"
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
hessBootstrap = function(object, nBoot = 500, newTime = NULL) {
  if (is.null(object$bootType)) {
    object$bootType = "hessian"
  } else if (trend$bootType != "hessian") {
    stop("Can't mix different bootstrap methods.")
  }
  if (object$trendType != 'index')
    XT = mgcv::predict.gam(object, newdata = object$trendFrame, type = "lpmatrix")
  XI = mgcv::predict.gam(object, newdata = object$indexFrame, type = "lpmatrix")
  if (object$trendType == "smooth") {
    tCol = grep(extractGamLabel(object$timeVar, mgcv::formula.gam(object)), colnames(XI), fixed = TRUE)
    tFacCol = integer(0)
  }
  if (object$trendType == "loglinear") {
    tCol = which(colnames(XI) == object$timeVar)
    tFacCol = integer(0)
  }
  if (object$trendType == "index") {
    tCol = integer(0)
    tFacCol = c(grep("(Intercept)", colnames(XI), fixed = TRUE), grep(object$timeVarFac, colnames(XI), fixed = TRUE))
  }
  if(object$timeRE) {
    tFacCol = grep(extractGamLabel(object$timeVarFac, mgcv::formula.gam(object)), colnames(XI), fixed = TRUE)
  } 
  XI = XI[, c(tCol, tFacCol), drop = FALSE]
  if (object$trendType != 'index') {
    XT = XT[, c(tCol, tFacCol), drop = FALSE]
    stopifnot(identical(colnames(XT), colnames(XI))) #DEBUG
  }
  cf = coef(object)[c(tCol, tFacCol)]
  vc  = mgcv::vcov.gam(object)[c(tCol, tFacCol), c(tCol, tFacCol)]
  chvc = t(chol(vc))
  bootI = bootI = NULL
  bootT = NULL
  bootIRes = NULL
  cfn = cf + chvc %*% matrix(rnorm(length(cf) * nBoot), ncol = nBoot, nrow = length(cf))
  if (object$trendType != "index") {
    bootT = XT[, seq_len(length(tCol)), drop = FALSE] %*% cfn[seq_len(length(tCol)), , drop = FALSE]
    bootI = XI[, seq_len(length(tCol)), drop = FALSE] %*% cfn[seq_len(length(tCol)), , drop = FALSE]
    if (object$timeRE) {
      bootIRes = XI[, length(tCol) + seq_len(length(tFacCol)), drop = FALSE] %*% 
                                     cfn[length(tCol) + seq_len(length(tFacCol)), , drop = FALSE]
    }
    
  }  else { # Standardize each simulation to avoid numerical problems. Could also be done for non-index trends. TODO.
    #XC = XI[, length(tCol) + seq_len(length(tFacCol)), drop = FALSE] %*% 
    #  cfn[length(tCol) + seq_len(length(tFacCol)), , drop = FALSE]
    #bdtFac = cbind(bdtFac, exp(XC-kronecker(rep(1,nrow(XC)), t(colSums(XC)/nrow(XC)))))
    bootI = XI[, length(tCol) + seq_len(length(tFacCol)), drop = FALSE] %*% 
                cfn[length(tCol) + seq_len(length(tFacCol)), , drop = FALSE]
  }
  
  object$boot = list(index = bootI, trend = bootT, residuals = bootIRes)
  object
}

setBaseline = function(object, baseline = NA, level = .95) {
  timeVar = object$timeVar
  if (length(baseline) == 1) { 
    if (is.na(baseline)) {
      baseline = min(object$indexFrame[[timeVar]])
    }
  }
  if (!all(baseline %in% object$indexFrame[[timeVar]])) {
    stop('baseline can only contain time points with observations')
  }
  bootRefValue = NULL
  if (length(baseline) == 1) {
    trendRefValue = object$indexFrame$indexRaw[object$indexFrame[[timeVar]] == baseline]
    if(!is.null(object$boot$index)) {
      bootRefValue = object$boot$index[object$indexFrame[[timeVar]] == baseline,]
    }
  } else {
    trendRefValue = log(mean(exp(object$indexFrame$indexRaw[object$indexFrame[[timeVar]] %in% baseline])))
    if(!is.null(object$boot$index)) {
      bootRefValue = log(apply(exp(object$boot$index[object$indexFrame[[timeVar]] %in% baseline, ]), 2, mean))
    }
  }
  if (object$trendType != 'index') {
    object$trendFrame$index = exp(object$trendFrame$indexRaw - trendRefValue)
    object$indexFrame$index = exp(object$indexFrame$indexRaw - trendRefValue)
    if (object$timeRE) {
      object$indexFrame$index.re = exp(object$indexFrame$indexRaw + object$indexFrame$resid - trendRefValue)
    }
  } else {
    object$indexFrame$index = exp(object$indexFrame$indexRaw - trendRefValue)
  }
  
  object$baseline = baseline
  object$trendRefValue = trendRefValue
  
  object$boot$bootRefValue = bootRefValue
  
  
  object = updateConfLevel(object, level = level)
  
  object
}

updateConfLevel = function(object, level) {
  if (!is.numeric(level) | level <= 0 | level >= 1)
    stop('level need to be a number between 0 and 1.')
  if (is.null(object$boot)) {
    stop('No uncertainty estimate available.')
  }
  if (is.null(object$boot$bootRefValue)) {
    stop('Baseline not set.')
  }
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
    grad1 = getGradient(bT, order = 1)
    grad2 = getGradient(bT, order = 2)
    trendCI = data.frame(t(ciT), rowMeans(grad1 > 0), rowMeans(grad2 > 0))
    colnames(trendCI) = c('ciLow', 'ciUpp', 'grad1p', 'grad2p')
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
  object$indexCI = indexCI
  object$ciLevel = level
  object

}



