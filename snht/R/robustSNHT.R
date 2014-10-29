robustSNHT <- function(data, period
#      #Tukey estimator:
#      ,estimator=function(x){
#          fit = rlm( x ~ 1, psi=psi.bisquare, c=1.5, maxit=100)
#          return(c(fit$coeff, fit$s))        
#        }
      #Huber estimator:
      ,estimator=function(x, minObs=5){
          #minObs arbitrarilty set to 5, want to ensure a decent number of values
          x = x[!is.na(x)]
          if(length(x)<minObs) #Too many NA values, don't return a result
            return(c(NA,NA))
          if(max(table(x))>length(x)/2) #Too many duplicate values, MAD will be 0
            return(c(NA,NA))
          fit = MASS::huber(x)
          return(c(fit[[1]], fit[[2]]))        
        }
      )  
{
  #Data quality checks
  if(!is.numeric(data))
    stop("data must be numeric!")
  if(2*period>length(data))
    stop("period is too large to compute statistics!")
    
  #Compute rolling means for 1:period, 2:(period+1), etc.
  #by=2 calculates the stat for every day, i.e. every two obs.  Could speed up by increasing this value
  Means = rollapply(data, width=period, by=1, FUN=estimator)
  n = rollapply(data, width=period, by=1, FUN=function(x) sum(!is.na(x)) )
  #Right means start at observation period+2 (first obs to use is at period+1
  # and then means are to the right)
  rMeans = Means[(period+2):nrow(Means),]
  rN = n[(period+2):nrow(Means)]
  #Left means are same length as right means but are at the beginning instead of end
  lMeans = Means[1:nrow(rMeans),]
  lN = n[1:nrow(rMeans)]
  totMean = (lN*lMeans[,1]+rN*rMeans[,1])/(lN+rN)
  scores = data.frame(
    #tukeyR also computes variances.  Variance of differences is sum of variances.
    #Original test stat has (N/2*(mu_L-mu)^2+N/2*(mu_R-mu)^2 )/sigma, but it assumes
    #N/2 observations on each side.  We must adjust for differing counts due to NA's:
     score= (lN*(lMeans[,1]-totMean)^2 + rN*(rMeans[,1]-totMean)^2 ) /
          sqrt((lN*lMeans[,2]^2+rN*rMeans[,2]^2)/(lN+rN))
    ,leftMean=lMeans[,1], rightMean=rMeans[,1])
  
  #Add zeros for rows skipped at beginning/end:
  toBind = data.frame(score=rep(0,period), leftMean=0, rightMean=0)
  scores = rbind(toBind, scores, toBind)
  
  return(scores)
}
