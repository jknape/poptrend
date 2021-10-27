##' Simulate population survey data.
##'
##' Simulates simple count survey data with a non-linear trend, and site and temporal random effects. 
##' The logistic function is used to create a trend the reduces the expected population size to half 
##' its initial value over the time period.
##'
##' @param nyear The number of years in the simulated survey.
##' @param nsite The number of sites in the simulated survey
##' @param mu The expected mean of the counts at the start of the survey.
##' @param timeSD Standard deviation (at log-scale) of annual mean deviation from the trend.
##' @param siteSD Standard deviation (at log-scale) of simulated among site variation.
##' @param size The size parameter of the negative binomial distribution. Defaults to Inf in which case the data are Poisson distributed.
##' 
##' @return A data frame containing simulated data.
##' @export
##' @author Jonas Knape
simTrend = function(nyear = 30, nsite = 40, mu = 3, timeSD = 0.1, siteSD = 0.3, size = Inf){
  if (mu <= 0)
    stop("mu must be positive")
  nyear = nyear
  nsite = nsite
  yeareff = log(.5 + .5*plogis(20 * (nyear/2 - 1:nyear) / nyear)) + timeSD * rnorm(nyear)
  siteeff =  siteSD * rnorm(nsite)
  data = data.frame(count = NA, year = rep(1:nyear, nsite), site = factor(rep(paste0("site", 1:nsite), each = nyear)))
  data$count = rnbinom(nyear*nsite, mu  = exp(log(mu) + yeareff[data$year] + siteeff[rep( 1:nsite, each = nyear)]), size = size)
  data
}

##' Plot an estimated trend.
##'
##' The function plots an estimated trend or index, as well as estimates of any temporal random effects included in the
##' trend term.
##'
##' Trends and indexes are relative measures and therefore are compared against some reference value.
##' By default, the first observed time point is used as the reference value.
##'
##' If the estimated trend contains bootstrap samples, confidence intervals are plotted as well.
##' For smooth trend models, time periods where the trend is significantly declining or increasing are marked with
##' a different color (set by arguments \var{decCol} and \var{incCol}). Periods where the second derivative is 
##' significantly positive or negative are marked by coloured boxes at the bottom of the plot. 
##' 
##' 
##' There is an additional option of plotting each of the bootstrapped trends.
##' @param x A fitted object of class trend.
##' @param baseline A single time point, or vector of time points, that are used to define the reference 
##'               level for the index. Can only contain time points for which there are observations in the data.
##'               If baseline is a vector of length larger than one, the estimated mean over all the time 
##'               points is taken as the reference. Defaults to NULL in which case the current baseline of the
##'               x is used (often the first time point in the data).
##' @param ylab The label of the y-axis.
##' @param alpha The alpha level of confidence intervals.
##' @param trendCol The color of the trend line.
##' @param lineCol The color of bootstrapped trend lines, if plotted.
##' @param shadeCol The color of the confidence region.
##' @param incCol The color of regions where the first or second derivative is significantly increasing.
##' @param decCol The color of regions where the first or second derivative is significantly decreasing.
##' @param plotGrid If true, grid lines are plotted.
##' @param gridCol Color of grid lines.
##' @param plotLines If true, the bootstrapped trends are plotted.
##' @param ranefCI If true, confidence intervals for random effects are plotted.
##' @param lwd Width of the trend line for smooth or loglinear trends.
##' @param ... Further arguments passed to \code{\link[graphics]{plot.default}}.
##' @export
##' @author Jonas Knape
plot.trend = function(x, baseline = NULL, ylab = "abundance index", trendCol = "#0072B2", lineCol = adjustcolor("#0072B2", alpha.f = 0.05), 
                      shadeCol = adjustcolor("#0072B2", alpha.f = 0.4), incCol = "#009E73", decCol ="#d90211",
                      plotGrid = TRUE, gridCol = 'gray80', plotLines = FALSE, ranefCI = TRUE, lwd = 3, ...) {
  timeVar = x$timeVar
  if (!is.null(baseline)) {
    if (!identical(baseline, x$baseline))
      x = updateBaseline(x, baseline = baseline)
  }
  if (x$trendType != 'index') {
    tGrid = x$trendFrame[[timeVar]]
    trendEst = x$trendFrame[['index']]
  } else {
    tGrid = x$indexFrame[[timeVar]]
    trendEst = x$indexFrame[['index']]
  }
  
  level = attr(x$indexCI, 'level')
  ci = NULL
  if (!is.null(x[['trendCI']])) {
    pGradInd = which(x$trendCI$ci.low.d1 > 0)
    nGradInd = which(x$trendCI$ci.upp.d1 < 0)
#    if (x$trendType == "smooth" & secDeriv) {
#      pGrad2Ind = getRuns(which(x$trendCI$grad2p > .5 + level/2))
#      nGrad2Ind = getRuns(which(x$trendCI$grad1p < .5 - level/2))
#    } else {
#      pGrad2Ind = nGrad2Ind = NULL
#    }
    ci$low = lowess(tGrid, x$trendCI[['ci.low.index']], f = .03)$y
    ci$upp = lowess(tGrid, x$trendCI[['ci.upp.index']], f = .03)$y
  }
  cip = NULL
  if (x$timeRE) {
    timeVarFac = x$timeVarFac
    cip = rbind(x$indexCI$ci.low.re, x$indexCI$ci.upp.re)
  }
  if(x$trendType == "index") {
    timeVarFac = x$timeVarFac
    cip = rbind(x$indexCI$ci.low.index, x$indexCI$ci.upp.index)
  }   
  ## Start plotting
  plotArgs = list(x=tGrid, y =  trendEst, type = "n", ylab = ylab, ...)
  if (!("xlab" %in% names(plotArgs))) {
    plotArgs$xlab = x$timeVar
  } 
  if (!("ylim" %in% names(plotArgs))) {
    plotArgs$ylim = c(min(trendEst, ci$low, cip[1,], na.rm=TRUE), max(trendEst, ci$upp, cip[2,], na.rm=TRUE))
  }
  do.call(plot, plotArgs)
  if (plotGrid)
    grid(nx = NA, ny = NULL, col = gridCol, lty = 1)
  if (!is.null(x$boot$trend)) {
    bymin = par()$usr[3]
    bymax = bymin + .1 * (par()$usr[4] - bymin)
    if (plotLines) {
      for (i in 1:ncol(x$boot$trend)) {
        points(tGrid, exp(x$boot$trend[,i] - x$boot$bootRefValue[i]), type = "l", 
               col = lineCol)
      }
    }
    polygon(c(tGrid,rev(tGrid)), c(ci$low, rev(ci$upp)), col = shadeCol, border = NA)
  }
  if (x$trendType != "index")
    points(tGrid, trendEst , type = "l", lwd = lwd, col = trendCol)
  if (!is.null(x$trendCI)) {
    segments(tGrid[replace(pGradInd-1,pGradInd==1,1)], trendEst[replace(pGradInd-1,pGradInd==1,1)],
             tGrid[pGradInd], trendEst[pGradInd], lwd = lwd, col = incCol)
    segments(tGrid[pGradInd], trendEst[pGradInd],
             tGrid[pGradInd+1], trendEst[pGradInd+1], lwd = lwd, col = incCol)
    segments(tGrid[replace(nGradInd-1,nGradInd==1,1)], trendEst[replace(nGradInd-1,nGradInd==1,1)],
             tGrid[nGradInd], trendEst[nGradInd], lwd = lwd, col = decCol)
    segments(tGrid[nGradInd], trendEst[nGradInd],
             tGrid[nGradInd+1], trendEst[nGradInd+1], lwd = lwd, col = decCol)
  }
  if (x$timeRE | x$trendType == "index") {
    if(!is.null(x$boot$residuals)) {
      if (plotLines) {
        for (i in 1:ncol(x$boot$residuals)) {
          #if (x$timeRE)
            points(x$indexFrame[[timeVar]], exp(x$boot$index[,i] + x$boot$residuals[, i] - x$boot$bootRefValue[i]), type = "p", pch = 20,cex = .5,
                   col = lineCol)
          #else
          #  points(resGrid, x$bootResid[ind, i] / bDiv[i], type = "p", pch = 20,cex = .5,
           #        col = lineCol)
        }
      }
    }
    if (ranefCI | x$trendType == "index") { # Plot confidence intervals for random effects or index.
        apply(cbind(x$indexFrame[[timeVar]], t(cip )), 1, 
              function(row) lines(x = c(row[1], row[1]), y = row[2:3], lwd = 1, col = trendCol))
    }
    if (x$timeRE)
      points(x$indexFrame[[timeVar]], x$indexFrame[['index.re']] , type = "p", pch = 20, col = trendCol)
    else
      points(x$indexFrame[[timeVar]], x$indexFrame[['index']] , type = "p", pch = 20, col = trendCol)
  }
  invisible(list(x$trendFrame, x$indexFrame))
}


##' Computes the estimated percentual change in the population between two given time points, 
##' and an approximate bootstrap confidence interval for the change.
##'
##' The function computes the estimated change between two chosen time points. 
##' When random effects are present, the change is computed for the underlying linear or
##' smooth trend term. 
##' For index models, the change is estimated from the difference between indices.
##' Changes can only be computed between time points that were included in the \code{trendGrid}
##' argument to \link{ptrend}, if the two time points are not included the nearest points in the grid are chosen.
##' 
##' Confidence intervals are computed using quantiles of the bootstrapped trends.
##'
##' @title Compute the change in the population over a time interval.
##' @param trend A fitted object of class trend.
##' @param start Start time for the comparison.
##' @param end End time for the comparison.
##' @param level Approximate level for confidence intervals.
##' @param nBoot Number of bootstrap samples to draw. For high precision of confidence intervals, a large number is recommended.
##'          Defaults to 2000. 
##' @param type Determines whether the change is computed as a difference at the log scale, as a multiplicative factor, or as the percentage change. 
##'          One of 'factor', 'log' or 'percent'.
##' @param pointwise If TRUE the change between every pair of time points in start and end are computed. 
##'           If FALSE the change between the mean over the time points in start and the mean over the time points in end is computed. 
##'           If both start and end are scalar, this argument is of little consequence.
##' @return A list containing the estimated change, and start and end points.
##' @note If \code{start} or \code{end} are not contained in the \code{trendgrid} argument of the \code{\link{ptrend}} function, 
##' the change is computed between the values in the grid that are closest to these points.
##' @export
##' @author Jonas Knape
##' @examples
##' ## Simulate a data set with 10 sites and 30 years
##' data = simTrend(30, 10)
##' ## Fit a smooth trend with fixed site effects, random time effects,
##' ## and automatic selection of degrees of freedom
##' trFit = ptrend(count ~ trend(year, type = "smooth") + site, data = data)
##' ## Check the estimated percent change from year 2 to 20
##' change(trFit, 10, 20)
change = function(trend, start, end, level = .95, nBoot = 2000, type = 'percent', pointwise = FALSE) {
  type = match.arg(type, c('factor', 'log', 'percent'))
  if (trend$trendType == 'index' & !all(c(start, end) %in% trend$indexFrame[[trend$timeVar]]))
    stop('For index models, all time points in start and end need to correspond exactly to time points in observations.')
  if (nBoot < trend$boot$nBoot) {
    cat('Number of bootraps in object (' ,  trend$boot$nBoot, ') is larger than  nBoot supplied (' , nBoot  ,'). Increasing nBoot to the larger value.')
    nBoot = trend$boot$nBoot
  }
  newData = trend$indexFrame[rep(1, length(c(start, end))),]
  newData$indexRaw = newData$index = newData$resid = newData$index.re = NULL
  if (trend$trendType == 'index') {
    newData[[trend$timeVarFac]] = trend$indexFrame[[trend$timeVarFac]][pmatch(c(start,end), trend$indexFrame[[trend$timeVar]])]
  } else {
    newData[[trend$timeVar]] = c(start, end)
  }
  ## Note: Ignoring random effects!
   
  sInd = seq_along(start)
  eInd = max(sInd) + seq_along(end)
  ct = computeRawTrend(trend, newdata = newData, estimate = TRUE, resid = FALSE, nBoot)
  
  if (pointwise)
    getChange = function(tr) outer(tr[eInd], tr[sInd], FUN = '-')
  else
    getChange = function(tr) log(mean(exp(tr[eInd]))) -log(mean(exp(tr[sInd])))
  ch = getChange(ct$indexRaw)
  bootch = apply(ct$bootIndex, 2, getChange)
  probs = c(1 - level, 1 + level)/2
  
  if (is.matrix(ch) & prod(dim(ch)) > 1) {
    colnames(ch) = start
    rownames(ch) = end
    ci = apply(bootch, 1,quantile, probs = probs)
    pr = rownames(ci)
    est = array(c(ch, ci[1,], ci[2,]),dim = c(dim(ch), 3))
    dimnames(est) = c(dimnames(ch), list(c('log-change', pr)))
  } else {
    ci = quantile(bootch, probs = probs)
    est = c(ch, ci)
  }
  if (type == 'percent') {
    est = 100 * (exp(est) - 1)
  }
  if (type == 'factor') {
    est = exp(est)
  }

  out = list(change = est, level = level, start = start, end = end, type = type, nBoot = nBoot, timeVar = trend$timeVar, pointwise = pointwise)
  class(out) = 'indexchange'
  out
}

print.indexchange = function(object) {
  if (is.vector(object$change)) {
    mc = ifelse(object$pointwise, ' from ' , ' in the mean from ')
    cat("Estimated ", object$type , " change" , mc, object$timeVar, " = ",  paste(object$start, collapse = ', '), " to ", object$timeVar, " = ", 
        paste(object$end, collapse = ', '), ": ", format(object$change[1], digits = 2), " ", 
        paste0("(", format(object$change[2], digits = 2),", ", format(object$change[3], digits =2), ")"), "\n",sep = "")
  } else {
    print(object$change[,,1])
  }
}

summary.indexchange = function(object) {
  object$change
}




# Computes derivatives along rows of mat using finite differences.
# Assumes regular grid, and that grid points are close enough.
fderiv = function(mat, grid, order = 1) {
  nr = nrow(mat)
  h = grid[2]-grid[1]
  if (any(abs(diff(grid) - h) > .Machine$double.base^0.5))
    stop('Irregular grid.')
  fdiff = matrix(NA_real_, nrow = nr, ncol = ncol(mat))
  if (order == 1) {
    fdiff[2:(nr-1), ] =  (mat[3:nr,] - mat[1:(nr-2),])/(2*h)
    fdiff[1, ] = (mat[2,]-mat[1,])/h
    fdiff[nr, ] = (mat[nr,]-mat[nr-1,])/h
  }
  if (order == 2) { 
    fdiff[2:(nr-1), ] = (mat[3:nr,] - 2 *mat[2:(nr-1),] + mat[1:(nr-2),]) / h^2
    fdiff[1, ] = fdiff[2,]
    fdiff[nr, ] = fdiff[nr-1,]
   }
  fdiff
}

# Compute start and endpoints for runs of 
# 1's in binary sequence.
getRuns = function(index, threshold = 1) {
  if (length(index) == 0) 
    return(NULL)
  rises = which(diff(c(0,index)) >= threshold)
  falls = which(diff(c(index, 0)) <= -threshold)
  cbind(rises, falls)
}


##' Prints basic information about a trend object.
##' 
##' Prints the family, formula and type of trend.
##' @title Print a trend object.
##' @param x A trend object.
##' @param ... Not used.
##' @export
##' @author Jonas Knape
print.trend = function(x, ...) {
  print(x$family)
  cat("Formula: ")
  print(x$formula)
  cat("Trend type: ", x$trendType)
  cat("\n")
  invisible(x)
}

##' Computes a index estimate for each time point in the survey.
##'
##' For a smooth or loglinear trend model the function computes an index from
##' the estimated trend line at each time point in the survey. By default, the reference
##' value is the first time point. Note that if the trend model was fitted with random 
##' effects, the random effects are not included in the index. Thus the index refers
##' to the long-term component of the trend.
##' 
##' For an index trend model the index at each time point is computed.
##' 
##' @title Summary of trend estimates
##' @param object A trend object returned by \code{\link{ptrend}}.
##' @param baseline 
##' @param ... Not used.               
##' @export
##' @author Jonas Knape
summary.trend = function(object, baseline = NULL,  ...) {
  timeVar = object$timeVar
  if (!is.null(baseline)) {
    if (!identical(baseline, object$baseline))
      object = updateBaseline(object, baseline = baseline)
  }
  
  level = attr(object$indexCI, 'level')  
  
  df = cbind(object$indexFrame[, c(timeVar, 'index')], object$indexCI[, c('ci.low.index', 'ci.upp.index')])
  names(df)[names(df) == 'ci.low.index'] =paste0("  ",format((1 - level)/2 * 100, digits = 3), "%")
  names(df)[names(df) == 'ci.upp.index'] =paste0("  ",format((1 + level) /2 * 100, digits = 3), "%")
  
  # if (object$timeRE) {
  #   df = cbind(df, object$indexFrame[, 'index.re', drop = FALSE], object$indexCI[, c('ci.low.re', 'ci.upp.re')])
  #   names(df)[names(df) == 'index.re'] = 'RE-index'
  #   names(df)[names(df) == 'ci.low.re'] = paste0("REI",format((1 - level)/2 * 100, digits = 3), "%")
  #   names(df)[names(df) == 'ci.upp.re'] = paste0("REI",format((1 + level) /2 * 100, digits = 3), "%")
  # }
  
  out = list(formula = object$formula, family = object$family, 
             trendType = object$trendType, estimates = df)
  class(out) = "summary.trend"
  out
}

#' @export
print.summary.trend = function(x, ..., digits = 2) {
  #browser()
  cat("Family: ")
  cat(x$family$family,'\n')
  cat("Formula: ")
  print(x$formula)
  cat("Trend type: ", x$trendType)
  cat("\n\n")
#  if (x$trendType != "index")
#    cat("Trend estimates:\n\n")
#  else
    cat("Index estimates:\n\n")
  print(x$estimates, ..., digits = digits, row.names = FALSE)
  invisible(x)
}
