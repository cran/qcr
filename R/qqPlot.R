qqPlot <- function (x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1,
                    start, ...) 
{
  DB = FALSE
  parList = list(...)
  if (is.null(parList[["col"]])) 
    parList$col = 1:2
  if (is.null(parList[["pch"]])) 
    parList$pch = 19
  if (is.null(parList[["lwd"]])) 
    parList$lwd = 1
  if (is.null(parList[["cex"]])) 
    parList$cex = 1
  
  #if (!require(MASS)) 
  #    stop("Package MASS needs to be installed!")
  
  if (class(x) == "distrCollection") {
    distList = x@distr
    for (i in 1:length(distList)) {
      d = distList[[i]]
      do.call(qqPlot, c(list(x = d@x, y = d@name), parList))
    }
    invisible()
  }
  if (missing(y)) 
    y = "normal"
  if(missing(alpha))
    alpha = 0.05
  if (alpha <=0 || alpha >=1) 
    stop(paste("alpha should be between 0 and 1!"))		
  if (missing(main)) 
    main = paste("Q-Q Plot for", deparse(substitute(y)), 
                 "distribution")
  if (missing(xlab)) 
    xlab = paste("Quantiles for", deparse(substitute(x)))
  if (missing(ylab)) 
    ylab = paste("Quantiles from", deparse(substitute(y)), 
                 "distribution")
  if (is.numeric(y)) {
    cat("\ncalling (original) qqplot from namespace stats!\n")
    return(stats::qqplot(x, y, ...))
  }
  qFun = NULL
  theoretical.quantiles = NULL
  xs = sort(x)
  distribution = tolower(y)
  distWhichNeedParameters = c("weibull", "logistic", "gamma", 
                              "exponential", "f", "geometric", "chi-squared", "negative binomial", 
                              "poisson")
  
  
  # new
  threeParameterDistr = c("weibull3", "lognormal3", "gamma3")                
  threeParameter = distribution %in% threeParameterDistr
  if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
  # end new
  
  if (is.character(distribution)) {
    qFun = .charToDistFunc(distribution, type = "q")
    if (is.null(qFun)) 
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  theoretical.probs = ppoints(xs)
  
  xq = NULL
  yq = quantile(xs, prob = c(0.25, 0.75))
  dots <- list(...)
  if (TRUE) {
    if (DB) 
      print("TODO: Pass the estimated parameters correctly")
    fitList = .lfkp(parList, formals(qFun))
    fitList$x = xs
    fitList$densfun = distribution
    if (!missing(start)) 
      fitList$start = start
    if (DB) {
      print(fitList)
      print("Ende")
    }
    # new
    if(!threeParameter){
      fittedDistr = do.call(fitdistr, fitList)
      parameter = fittedDistr$estimate
      
      #save the distribution parameter#
      thethas = fittedDistr$estimate
      # save the cariance-covariance matrix
      varmatrix = fittedDistr$vcov
      # end of my code
      
      # new code for three parameter
    } else {
      parameter = do.call(paste(".",distribution, "3", sep = ""), list(xs) )    ####
      threshold = parameter$threshold
    }
    
    parameter = .lfkp(as.list(parameter), formals(qFun))
    params = .lfkp(parList, formals(qFun))
    parameter = .lfrm(as.list(parameter), params)
    parameter = c(parameter, params)
    theoretical.quantiles = do.call(qFun, c(list(c(theoretical.probs)), 
                                            parameter))
    
    # new
    if(!threeParameter){		
      # array containing names of the distributions, for which conf intervals can be computed
      confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
      getConfIntFun = .charToDistFunc(distribution, type = ".confint")
      # if possible, compute the conf intervals
      if(confbounds == TRUE){
        if(distribution %in% confIntCapable){
          confInt = getConfIntFun(xs, thethas, varmatrix, alpha)
        }
      }# end of my code
    }
    
    xq <- do.call(qFun, c(list(c(0.25, 0.75)), parameter))
    if (DB) {
      print(paste("parameter: ", parameter))
      print(xq)
    }
  }
  else {
    params = .lfkp(parList, formals(qFun))
    params$p = theoretical.probs
    theoretical.quantiles = do.call(qFun, params)
    params$p = c(0.25, 0.75)
    xq = do.call(qFun, params)
  }
  
  params = .lfkp(parList, c(formals(plot.default), par()))	
  
  if(!threeParameter){
    params$y = theoretical.quantiles
  }  else {
    params$y = theoretical.quantiles+threshold
  }
  params$x = xs
  params$xlab = xlab
  params$ylab = ylab
  params$main = main
  if (!(is.null(params$col[1]) || is.na(params$col[1]))) 
    params$col = params$col[1]
  if (!missing(xlim)) 
    params$xlim = xlim
  if (!missing(ylim)) 
    params$ylim = ylim
  params$lwd = 1
  do.call(plot, params)
  pParams = params
  pParams = .lfkp(pParams, list(x = 1, y = 1, col = 1, cex = 1))
  do.call(points, pParams)
  params = .lfkp(parList, c(formals(abline), par()))
  params$a = 0
  params$b = 1
  params$col = border
  do.call(abline, params)
  
  if(!threeParameter){
    # plot the confInt if available
    if(confbounds == TRUE){
      if(distribution %in% confIntCapable){
        params = .lfkp(parList, c(formals(lines), par()))	
        params$x = confInt[[3]]
        params$y = confInt[[1]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)
        
        params$x = confInt[[3]]
        params$y = confInt[[2]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)
      }
    } #end of my function
  }
  
  invisible(list(x = theoretical.quantiles, y = xs, int = params$a, 
                 slope = params$b))
}