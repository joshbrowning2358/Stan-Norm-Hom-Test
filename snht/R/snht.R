snht <- function(data, period, robust=F, time=NULL, ...){
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(2*period>length(data))
    stop("period is too large to compute statistics!")
  if(!is.null(time))
    if(length(time)!=length(data))
      stop("If time is not NULL, it must be the same length as data!")
  
  if(!is.null(time)){
    if(robust)
      out = robustSNHTunequal( data=data, period=period, time=time, estimator=NULL )
    else
      out = robustSNHTunequal( data=data, period=period, time=time
                              ,estimator=function(x){ #Use the mean and sd for non-robust fit
                                  x = x[!is.na(x)]
                                  return( c(mean(x), sd(x)) )
                              } )
    return(out)
  }
  
  if(!robust){
    #Use the "Robust" SNHT function, but supply a non-robust function: mean and sd.
    out = robustSNHT( data, period=period, estimator=function(x){
      x = x[!is.na(x)]
      return( c(mean(x), sd(x)) )
    } )
  } else {
    out = robustSNHT(data, period, ...)
  }
  colnames(out) = c("score", "rightMean", "leftMean")
  return(out)
}
