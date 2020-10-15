##' Simulate population survey data.
##'
##' Simulates count survey data with a non-linear trend, and site and temporal random effects. 
##' The logistic function is used to create a trend the reduces the expected population size to half 
##' its initial value over the time period.
##'
##' @param nyear The number of years in the simulated survey.
##' @param nsite The number of sites in the simulated survey
##' @param mu The expected mean of the counts at the start of the survey.
##' @param timeSD Standard deviation (at log-scale) of annual mean deviation from the trend.
##' @param siteSD Standard deviation (at log-scale) of simulated among site variation.
##' @return A data frame containing simulated data.
##' @export
##' @author Jonas Knape
simTrend = function(nyear = 30, nsite = 40, mu = 3, timeSD = 0.1, siteSD = 0.3){
  if (mu <= 0)
    stop("mu must be positive")
  nyear = nyear
  nsite = nsite
  yeareff = log(.5 + .5*plogis(20 * (nyear/2 - 1:nyear) / nyear)) + timeSD * rnorm(nyear)
  siteeff =  siteSD * rnorm(nsite)
  data = data.frame(count = NA, year = rep(1:nyear, nsite), site = factor(rep(paste0("site", 1:nsite), each = nyear)))
  data$count = rpois(nyear*nsite, lambda  = exp(log(mu) + yeareff[data$year] + siteeff[rep( 1:nsite, each = nyear)]))
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
##' @param base A time point or function used to compute the baseline of the trend. 
##'               If the argument is numeric, the point in the \var{trendGrid} argument of the function \code{\link{ptrend}}
##'               closest to this value will be taken as the baseline (i.e. the estimated trend will be 1 at this point).
##'               If the argument is a function, the function is applied to trends and the resulting value is used as the baseline.
##'               By default, the first time point is taken as the reference.
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
##' @param ... Further arguments passed to \code{\link[graphics]{plot.default}}.
##' @export
##' @author Jonas Knape
plot.trend = function(x, base = NULL, ylab = "abundance index", trendCol = "black", lineCol = adjustcolor("black", alpha.f = 0.05), 
                      shadeCol = adjustcolor("#0072B2", alpha.f = 0.4), incCol = "#009E73", decCol ="#D55E00",
                      plotGrid = TRUE, gridCol = 'gray80', plotLines = FALSE, ranefCI = TRUE, ...) {
  timeVar = x$timeVar
  if (!is.null(base)) {
    if (!identical(base, x$baseline))
      x = updateBaseline(x, baseline = base)
  }
  if (x$trendType != 'index') {
    tGrid = x$trendFrame[[timeVar]]
    trendEst = x$trendFrame[['trend']]
  } else {
    tGrid = x$indexFrame[[timeVar]]
    trendEst = x$indexFrame[['index']]
  }
  
  level = x$ciLevel
  ci = NULL
  if (!is.null(x[['trendCI']])) {
    pGradInd = which(x$trendCI$grad1p > .5 + level/2)
    nGradInd = which(x$trendCI$grad1p < .5 - level/2)
#    if (x$trendType == "smooth" & secDeriv) {
#      pGrad2Ind = getRuns(which(x$trendCI$grad2p > .5 + level/2))
#      nGrad2Ind = getRuns(which(x$trendCI$grad1p < .5 - level/2))
#    } else {
#      pGrad2Ind = nGrad2Ind = NULL
#    }
    ci$low = lowess(tGrid, x$trendCI$ciLow, f = .03)$y
    ci$upp = lowess(tGrid, x$trendCI$ciUpp, f = .03)$y
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
  if (!is.null(x[['bootT']])) {
    bymin = par()$usr[3]
    bymax = bymin + .1 * (par()$usr[4] - bymin)
    if (plotLines) {
      for (i in 1:ncol(x$bootT)) {
        points(tGrid, exp(x$bootT[,i] - x$bootRefValue[i]), type = "l", 
               col = lineCol)
      }
    }
    polygon(c(tGrid,rev(tGrid)), c(ci$low, rev(ci$upp)), col = shadeCol, border = NA)
  }
  if (x$trendType != "index")
    points(tGrid, trendEst , type = "l", lwd = 3, col = trendCol)
  if (!is.null(x[['bootT']])) {
    segments(tGrid[replace(pGradInd-1,pGradInd==1,1)], trendEst[replace(pGradInd-1,pGradInd==1,1)],
             tGrid[pGradInd], trendEst[pGradInd], lwd = 3, col = incCol)
    segments(tGrid[pGradInd], trendEst[pGradInd],
             tGrid[pGradInd+1], trendEst[pGradInd+1], lwd = 3, col = incCol)
    segments(tGrid[replace(nGradInd-1,nGradInd==1,1)], trendEst[replace(nGradInd-1,nGradInd==1,1)],
             tGrid[nGradInd], trendEst[nGradInd], lwd = 3, col = decCol)
    segments(tGrid[nGradInd], trendEst[nGradInd],
             tGrid[nGradInd+1], trendEst[nGradInd+1], lwd = 3, col = decCol)
  }
  if (x$timeRE | x$trendType == "index") {
    if(!is.null(x$bootIRes)) {
      if (plotLines) {
        for (i in 1:ncol(x$bootIRes)) {
          #if (x$timeRE)
            points(x$indexFrame[[timeVar]], exp(x$bootI[,i] + x$bootIRes[, i] - x$bootRefValue[i]), type = "p", pch = 20,cex = .5,
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
    points(x$indexFrame[[timeVar]], x$indexFrame[['index']] , type = "p", pch = 20, col = trendCol)
  }
  invisible(list(x$trendFrame, x$indexFrame))
}


##' Computes the estimated percentual change in the population between two given time points, 
##' and an approximate confidence interval for the change.
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
##' @param alpha alpha-level for approximate confidence interval.
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
change = function(trend, start, end, alpha = .05) {
  if (trend$trendType != "index") {
    tGrid = trend$trendFrame[[trend$timeVar]]
    trendEst = trend$trendFrame$trendRaw
  } else {
    tGrid = trend$indexFrame[[trend$timeVar]]
    trendEst = trend$indexFrame$trendRaw
  }
  sInd = which.min(abs(start - tGrid))
  eInd = which.min(abs(end - tGrid))
  percChange = function(tr) 100 * (exp(tr[eInd] - tr[sInd]) - 1)
  pc = percChange(trendEst)
  if(trend$trendType != "index" & !is.null(trend$bootT)) {
    bootpc = apply(trend$bootT, 2, percChange)
    CI = quantile(bootpc, probs = c(alpha/2, 1 - alpha/2))
  } else if (trend$trendType == "index" & !is.null(trend$bootI)) {
    bootpc = apply(trend$bootI, 2, percChange)
    CI = quantile(bootpc, probs = c(alpha/2, 1 - alpha/2))
  } else
    CI = NULL
  cat("Estimated percent change from ", trend$timeVar, " = ",  tGrid[sInd], " to ", tGrid[eInd], ": ", format(pc, digits = 2), "% ", 
      ifelse(is.null(CI), "", paste0("(", format(CI[1], digits = 2),"%, ", format(CI[2], digits =2), "%)")), "\n",sep = "")
  
  invisible(list(percentChange = pc, CI = CI, start = tGrid[sInd], end = tGrid[eInd]))
}

# Only computes relative magnitude of derivative and assumes regular grid (no 1/h, 1/h^2)
getGradient = function(bootTrend, order = 1) {
  nr = nrow(bootTrend)
  if (order == 1) {
    #grad = apply(bootTrend, 2, function(bt) {diff(bt, lag = 2) / 2})
    #grad = rbind(bootTrend[2,]-bootTrend[1,],grad, bootTrend[nrow(bootTrend),]-bootTrend[nrow(bootTrend)-1,])
    # Below is faster
    grad =  (bootTrend[3:nr,] - bootTrend[1:(nr-2),])/2
    grad = rbind(bootTrend[2,]-bootTrend[1,],grad, bootTrend[nr,]-bootTrend[nr-1,])
    return(grad)
  }
  if (order == 2) {
    secDeriv = (bootTrend[3:nr,] - 2 *bootTrend[2:(nr-1),] + bootTrend[1:(nr-2),])
    secDeriv = rbind(bootTrend[3,]-2*bootTrend[2,] + bootTrend[1,], secDeriv, bootTrend[nr,]-2*bootTrend[nr-1,] + bootTrend[nr-2,])
    #apply(bootTrend, 2, function(bt) {diff(bt, differences = 2)})
    return(secDeriv)
  }
}

getRuns = function(index) {
  if (length(index) == 0) 
    return(NULL)
  breaks = which(diff(index) > 1)
  end = index[breaks]
  start = index[breaks + 1]
  cbind(c(min(index), start), c(end, max(index)))
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

##' Computes a trend or index estimate for each time point in the survey.
##'
##' For a smooth or loglinear trend model the function computes an estimate
##' of the trend value for each time point in the survey. By default, the reference
##' value is the first time point. Note that if the trend model was fitted with random 
##' effects, the random effects are not included in the estimate. Thus the estimate refers
##' to the long-term component.
##' 
##' For an index trend model the index at each time point is computed.
##' 
##' If bootstrap samples are available, bootstrap confidence intervals for the trend 
##' or index values are also computed.
##' @title Summary of trend estimates
##' @param object A trend object returned by \code{\link{ptrend}}.
##' @param baseline A time point or function used to compute the baseline of the trend. 
##'               If the argument is numeric, the point in the \code{trendGrid} argument of the function \code{\link{ptrend}}
##'               closest to this value will be taken as the baseline (i.e. the estimated trend will be 1 at this point).
##'               If the argument is a function, the function is applied to trends and the resulting value is used as the baseline.
##'               By default, the first time point is taken as the reference.
##' @param ... Not used.               
##' @param alpha alpha level for approximate confidence intervals.
##' @export
##' @author Jonas Knape
summary.trend = function(object, baseline = NULL, alpha = 0.05, ...) {
  browser()
  isTP = which(!object$trendFrame$isGridP)
  timeVar = object$timeVar
  tGrid = object$trendFrame[[timeVar]]
  if (object$trendType != "index")
    trendEst = object$trendFrame$trend
  else
    trendEst = object$trendFrame$trendResid
  df = data.frame(object$trendFrame[[timeVar]])
  names(df) = timeVar    
  if (object$trendType != "index")
    df$trend = object$trendFrame$trend[isTP] / tDiv
  if (!is.null(object$bootTrend)) {
    ciEst = t(apply(object$bootTrend[isTP, ], 1, function(row) quantile(row / bDiv, probs = c(alpha/2, 1-alpha/2), type = 1)))
    df[[paste0("  ",alpha/2 * 100, "%")]] = ciEst[,1]
    df[[paste0("  ",(1 - alpha/2) * 100, "%")]] = ciEst[,2]
    
  }
  if (object$trendType != "index") {
    #df$index = object$trendFrame$trendResid[isTP] * df$trend / tDiv
  } else {
    df$index = object$trendFrame$trendResid[isTP] / tDiv
    if(!is.null(object$bootResid)) {
      ciEst = t(apply(object$bootResid[isTP, ], 1, function(row) quantile(row / bDiv, probs = c(alpha/2, 1-alpha/2), type = 1)))
      df[[paste0("  ",alpha/2 * 100, "%")]] = ciEst[,1]
      df[[paste0("  ",(1 - alpha/2) * 100, "%")]] = ciEst[,2]   
    }
  }
  out = list(formula = object$formula, family = object$family, 
             trendType = object$trendType, estimates = df)
  class(out) = "summary.trend"
  out
}

#' @export
print.summary.trend = function(x, ..., digits = 2) {
  #browser()
  cat("Formula: ")
  print(x$formula)
  cat("Trend type: ", x$trendType)
  cat("\n\n")
  if (x$trendType != "index")
    cat("Trend estimates:\n\n")
  else
    cat("Index estimates:\n\n")
  print(x$estimates, ..., digits = digits, row.names = FALSE)
  invisible(x)
}

##' Produces various goodness of fit plots and diagnostic measures.
##'
##' The function simply calls \link[mgcv]{plot.gam} and \link[mgcv]{gam.check} on the 
##' underlying gam model for checking goodness of fit. 
##' @title Check goodness of fit of a trend model.
##' @param trend A fitted object of class trend.
##' @param residuals Should residuals be plotted?
##' @param ... Further arguments passed to \code{\link[mgcv]{plot.gam}}.
##' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link[mgcv]{gam.check}}
##' @export
##' @author Jonas Knape
checkFit = function(trend, residuals = TRUE, ...) {
  if(is.null(trend$gam))
    stop("gam fit not available. Try setting argument keepGam to TRUE in call to ptrend.")
  mgcv::plot.gam(trend$gam, residuals = residuals, ...)
  mgcv::gam.check(trend$gam)
}
