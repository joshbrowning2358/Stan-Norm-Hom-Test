robustSNHTunequal <-
function(data, period, time, estimator=NULL){
  #Data quality checks
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(2*period>length(data))
    stop("period is too large to compute statistics!")
  if(!is.numeric(time))
    stop("time must be numeric!")
  if(!is.null(estimator) & !is(estimator,"function"))
    stop("estimator must be a function or must be NULL!")
  if(length(data)!=length(time))
    stop("data and time must be of the same length!")
  if(any(is.na(time)))
    stop("time cannot have missing values!")
  if(any(time!=floor(time))){
    warning("Only integer values of time are used!  Rounding down.")
    time = floor(time)
  }

  d = data.frame(data=data, time=time, realObs=1)
  d = d[order(d$time),]
  #Bind on rows to d so we can see times that have no data as well
  dTemp = rbind(d, data.frame(data=NA, time=min(time):max(time), realObs=0) )
  obsPerDay = ddply(dTemp, "time", function(df){sum(df$realObs)})
  #Compute maximum observations for any time unit (called "day", but arbitrary)
  maxObs = max(obsPerDay[,2])
  obsPerDay$toBind = maxObs - obsPerDay[,2]
  #Determine which times are missing, and repeat the appropriate amount.  Then, bind on.
  toBind = rep( obsPerDay$time, times=obsPerDay$toBind )
  if(length(toBind)>0)
    d = rbind(d, data.frame(data=NA, time=toBind, realObs=0) )
  d = d[order(d$time),]
  
  if(is.null(estimator))
    out = robustSNHT(data=d[,1], period=period*maxObs)
  else
    out = robustSNHT(data=d[,1], period=period*maxObs, estimator)
  out$time = d$time
  return(out)
}
