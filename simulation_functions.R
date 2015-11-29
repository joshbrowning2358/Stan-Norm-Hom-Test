library(sn)
library(chron)
library(Rcpp)
library(inline)
library(RcppArmadillo)
library(ggplot2)
library(scales)
library(plyr)
library(reshape)
library(mgcv)
library(MASS)
#Change point libraries:
# library(changepoint)
# library(bcpmeta)


########################################################################################################
##
## Functions for simulating data.
## - fitBaseTS: Fit a GAM model that estimates Temp as a function of s(Hour), s(Day of Year) and Year.
## - simBaseTS: Simulate data and use the model from fitBaseTS to "simulate" the expected value.
## - fitSkewT: Fit a skew-t distribution to the residuals of the TS model.
## - simResiduals: Simulate correlated residuals from the fitted skew-t distribution.
## - simResiduals1: Same as simResiduals, uses vectorized code.
## - simResiduals2: Fastest version.  Implemented in C++.
## - accumulate: Helper function for simResiduals2.
## - contaminateOutlier: Contaminates the simulated data with some specified number of outliers.
## - contaminateBreak: Contaminates the simulated data with jumps/breaks in the data.
##
########################################################################################################

#Fit a GAM to understand how the data changes in time.  The model is essentially:
# data ~ s(Day_of_Year) + s(Hour) + yr + epsilon
#data: the dataset to fit the model to.  Should have Reading, Month, Year, and Hour as columns
fitBaseTS = function(data=test.data)
{
  data$Date = as.Date(paste(data$Year,data$Month, data$Day), format="%Y %m %d")
  data$Day_Of_Year = as.numeric( as.character( data$Date, "%j" ) )
  mod=gam( Reading ~ s(Day_Of_Year) + s(Hour) + Year, data=data)
  return(mod)
}

#Fits a skew-t distribution to the errors from the time-series fit
#fitTS: the object returned from fitBaseTS
#data: the original dataset
getSkewTParams = function(fitTS, data=test.data)
{
  ##--------------------------------
  # Residuals of Original Data
  # w/Skew-T Fitted to it
  ##--------------------------------
  
  res1=data$Reading-predict(fitTS, newdata=data)
  res=res1[is.na(res1)!=TRUE]
  #mst.fit seems to be deprecated:
  #ms=mst.fit(y=res)
  #instead use selm.fit
  ms = try(selm.fit(x=matrix(1,length(res)), y=as.numeric(res), family="ST",
    selm.control = list()))
  if(is(ms, "try-error")){ #Skew-t fails, try skew-normal
    ms = try(selm.fit(x=matrix(1,length(res)), y=as.numeric(res), family="SN",
      selm.control = list()))
    if(is(ms, "try-error")){ #Skew normal fails, estimate mean and sd and set alpha=0, nu=Inf
      out = c(xi=mean(res), omega=sd(res), alpha=0, nu=Inf)
    } else {
      out = c(ms$param$dp, nu=Inf)
    }
  } else {
    out = ms$param$dp
  }
  
  return(out)
}

#Compute acf for unequally spaced data.  This function assumes that data has
#Date and Hour columns which are used to compute time lags.  Currently, ACF is
#computed on a daily scale, and ACF values for 1 up to maxLags days are computed.
#R computes ACF by first computing auto-cov:
# gamma(h) = 1/n*sum_{t=1}^{n-h}(x_{t+h}-\bar{x})(x_t-\bar{x})
#Note that the sum has n-h terms but we divide by n.  Also, \bar{x} is taken
#over the entire time series for both cases, and hence uses all n data points.
#ACF(h)=gamma(h)/gamma(0).  Generalizing this to unequal spacing has complications:
#- What should the denominator be?  If we use n pairs, n+h?
#- Should we compute \bar{x} over the entire series?  Or only over all pairs used?
#I ignore these problems as they presumably become small as n gets large.
#data: the data.frame containing the relevant data.  Should have the TS data
# as well as Date and Hour.
#data.col: Which column of data contains the TS?
#maxLags: How many lags should we compute ACF for?
#lagScale: allows computation of non-daily lags.  ACF will be computed for (1:maxLags)*lagScale
#adjDenom: With R base's acf() function, the autocorrelation at lag h is estimated
# using n-h pairs, but the denominator used in n.  This ensures the estimator is 
# positive definite, although neither denominator leads to an unbiased estimator.
# Generalizing this to the case of unequal spacing is not clear: should the denominator
# be the number of pairs?  The number of pairs + h-1 (to match the acf function)?  The
# total number of observations?  If adjDenom=TRUE, the number of pairs + h-1 is used; if
# FALSE then the number of pairs is used.
acfUnequal = function( data, data.col=5, maxLags=10, lagScale=1, adjDenom=TRUE )
{
  data = data[!is.na(data[,data.col]),] #Can't use for acf anyways, and causes problems in calcs
  reqCols = c("Date", "Hour")
  test = reqCols %in% colnames(data)
  if( any(!test) )
    stop(paste0("Missing the following columns: ", paste(reqCols[!test],collapse=", ")))
  if(ncol(data)<data.col)
    stop("data.col is larger than ncol(data)")

  diffT = as.numeric(diff(data$Date)) + diff(data$Hour)/2400
  diffT = diffT/lagScale
  lags = acfUnequalC(data[,data.col], diffT, maxLags)
  colnames(lags) = c("ACF", "count")
  if(adjDenom)
    lags[,1] = lags[,1]*lags[,2]/(lags[,2]+1:nrow(lags))
  return( data.frame(Lag.Days=1:maxLags*lagScale, lags) )
}

#Helper function for acfUnequal.
#dataR: the TS data, should be a vector.
#diffTR: the time difference between consecutive observations, should be a vector.
#maxLagsR: the number of lags to compute acf for, should be a numeric value.
acfUnequalC = cxxfunction(signature(dataR="numeric", diffTR="numeric", maxLagsR="numeric"),
  body='
  arma::vec data = as<arma::vec>(dataR);
  arma::vec diffT = as<arma::vec>(diffTR);
  int maxLags = as<int>(maxLagsR);
  arma::vec Exy(maxLags);
  arma::vec Exx(maxLags);
  arma::vec Eyy(maxLags);
  arma::vec Ex(maxLags);
  arma::vec Ey(maxLags);
  arma::vec cnt(maxLags);
  Exy.fill(0);
  Exx.fill(0);
  Eyy.fill(0);
  Ex.fill(0);
  Ey.fill(0);
  cnt.fill(0);

  for(unsigned i=0; i<data.size(); i++)
  {
    double diffCurr = 0;
    int j = i;
    for(int lag=1; lag<=maxLags; lag++)
    {
      while( (diffCurr<lag-0.05) & (j>=1) )
      {
        j--;
        diffCurr += diffT[j];
      }
      if( (diffCurr<=lag+.05) & (diffCurr>=lag-0.05)) //time lag within 0.05 days detected, use data to compute acf
      {
        Exy[lag-1] += data[i]*data[j]; //lag-1 because indices go from 0:maxLags-1 but lag loops over 1 to maxLags
        Exx[lag-1] += data[i]*data[i];
        Eyy[lag-1] += data[j]*data[j];
        Ex[lag-1] += data[i];
        Ey[lag-1] += data[j];
        cnt[lag-1]++;
      }
    }
  }
  
  arma::vec sxy = (Exy - Ex%Ey/cnt)/(cnt); //% is element-wise multliplication, / is element-wise division
  arma::vec sxx = (Exx - Ex%Ex/cnt)/(cnt);
  arma::vec syy = (Eyy - Ey%Ey/cnt)/(cnt);
  arma::mat acf(sxy.n_elem, 2);
  acf.col(0) = sxy/sqrt(sxx%syy);
  acf.col(1) = cnt;
/*  arma::mat acf(sxy.n_elem, 7);
  acf.col(0) = sxy/sqrt(sxx%syy);
  acf.col(1) = cnt;
  acf.col(2) = Exx;
  acf.col(3) = Eyy;
  acf.col(4) = Exy;
  acf.col(5) = Ex;
  acf.col(6) = Ey;*/
  return(wrap(acf));
', plugin="RcppArmadillo" )

acfUnequalR = function( data, data.col=5, maxLags=10 )
{
  data = data[!is.na(data[,data.col]),] #Can't use for acf anyways, and causes problems in calcs
  reqCols = c("Date", "Hour")
  test = reqCols %in% colnames(data)
  if( any(!test) )
    stop(paste0("Missing the following columns: ", paste(reqCols[!test],collapse=", ")))
  if(ncol(data)<data.col)
    stop("data.col is larger than ncol(data)")

  diffT = as.numeric(diff(data$Date)) + diff(data$Hour)/2400
  Exy = rep(0,maxLags)
  Exx = rep(0,maxLags)
  Eyy = rep(0,maxLags)
  Ex = rep(0,maxLags)
  Ey = rep(0,maxLags)
  cnt = rep(0,maxLags)

  for(i in 1:nrow(data))
  {
    diffCurr = 0
    j = i
    for(lag in 1:maxLags)
    {
      while( (diffCurr<lag-0.05) & j>=2 )
      {
        j = j-1
        diffCurr = diffCurr + diffT[j];
      }
      if(diffCurr<=lag+.05 & diffCurr>=lag-.05) #time lag within 0.05 days detected, use data to compute acf
      {
        Exy[lag] = Exy[lag] + data[i,data.col]*data[j,data.col]
        Exx[lag] = Exx[lag] + data[i,data.col]*data[i,data.col]
        Eyy[lag] = Eyy[lag] + data[j,data.col]*data[j,data.col]
        Ex[lag] = Ex[lag] + data[i,data.col]
        Ey[lag] = Ey[lag] + data[j,data.col]
        cnt[lag] = cnt[lag] + 1
      }
    }
  }
  
  sxy = (Exy - Ex*Ey/cnt)/(cnt);
  sxx = (Exx - Ex*Ex/cnt)/(cnt);
  syy = (Eyy - Ey*Ey/cnt)/(cnt);
  acf = sxy/sqrt(sxx*syy);
  return(acf);
}

#Using the model from fitBaseTS, simulate data.
#mod: the model returned by fitBaseTS
#data: the data used for fitBaseTS
simBaseTS = function(mod, data=test.data, startYear=1955, endYear=2011)
{
  ##--------------------------------
  # Creating Time Variables
  # to Simulate With
  ##--------------------------------
  #57 years worth of data
  year=array()
  for(i in startYear:endYear){
  	year=c(year,rep(i,365*2))
  }
  year=year[-1]
  sim.n=length(year)
  
  #one.months=seq(1,12,len=365*2)
  one.months=c(rep(1,31*2),rep(2,28*2),rep(3,31*2),rep(4,30*2),rep(5,31*2),rep(6,30*2),rep(7,31*2),rep(8,31*2),rep(9,30*2),rep(10,31*2),rep(11,30*2),rep(12,31*2))
  many.months=rep(one.months,endYear-startYear+1)
  day=rep(rep(1:365,each=2),times=length(unique(year)))
  hourDay=sample(test.data$Hour[test.data$Hour>=500 & test.data$Hour<1700 & !is.na(test.data$Hour)]
          ,sim.n/2,replace=TRUE)
  hourNgt=sample(test.data$Hour[test.data$Hour>=1700 | test.data$Hour < 500 & !is.na(test.data$Hour)]
          ,sim.n/2,replace=TRUE)
  hour = c(hourDay, hourNgt)
  hour = hour[rep(1:(sim.n/2),each=2) + c(0,sim.n/2)] #reorder so we have day, night, day, night, etc.
    
  sim.dat=data.frame(Year=year, Month=many.months, Day_Of_Year=day, Hour=hour)
  ##--------------------------------
  # Fitted Model with New Time
  # Variables
  ##--------------------------------
  
  sim.temp=as.numeric( predict(mod, newdata=sim.dat) )
  sim.dat$sim.temp=sim.temp
  sim.dat$Date = as.POSIXct( paste(sim.dat$Year, sim.dat$Day_Of_Year
                    ,floor(sim.dat$Hour/100), round((sim.dat$Hour %% 100)*60/100,0) )
                ,format="%Y %j %H %M")
  
  #If you hit daylight savings time, it's possible that a particular hour does not exist.  Adding 1 hour fixes the problem
  if(any(is.na(sim.dat$Date))){
    sim.dat$Hour[is.na(sim.dat$Date)] = sim.dat$Hour[is.na(sim.dat$Date)] + 100
    sim.dat$Date = as.POSIXct( paste(sim.dat$Year, sim.dat$Day_Of_Year
                      ,floor(sim.dat$Hour/100), round((sim.dat$Hour %% 100)*60/100,0) )
                  ,format="%Y %j %H %M")
  }
  
  sim.dat = sim.dat[order(sim.dat$Year, sim.dat$Day_Of_Year, sim.dat$Hour),]
  diffT = as.numeric(diff(sim.dat$Date)) + as.numeric(diff(sim.dat$Hour))/2400
  #Some obs have diffT==0, remove those and any that have diffT<0 (should never happen)
  sim.dat = sim.dat[c(T,diffT>0),]

  return(sim.dat)
}

# Simulating Skew-T residuals w/AR(1) structure
#fitST: the skew-t parameters fit to the residuals
#sim.dat: The simulated data without errors (from simBaseTS())
#ar1: the coefficient for the first auto-regressive term
simResiduals2 = function(fitST, sim.dat, ar1=0.15)
{
  if(!"sim.temp" %in% colnames(sim.dat) )
    stop("sim.dat missing required columns!")
  
  diffT = as.numeric(diff(sim.dat$Date)) + as.numeric(diff(sim.dat$Hour))/2400
  resVector = rst(1000+nrow(sim.dat),xi=fitST[1], omega=fitST[2], alpha=fitST[3], nu=fitST[4])
  #The AR model will increase variance of the errors:
  #Err_i <- Err_i + ar1*Err_{i-1} + ar1^2*Err_{i-2} + ...
  #Var(Err_i) <- Var(Err_i) + ar1^2*Var(Err_{i-1}) + ar1^4*Var(Err_{i-2}) + ...
  # <- Var(Err_i)/(1-ar1^2)
  #So, scale down the simulated errors by this proportion
  #resVector=resVector*(1-ar1^2)
  #But, this assumes equal spacing, which is also incorrect...
  update = accumulate(matrix(resVector), ar1, c(rep(1,1000),diffT))
  update = update[-(1:1000)]
  sim.dat$sim.temp=update+sim.dat[,"sim.temp"]
  colnames(sim.dat)[colnames(sim.dat)=="sim.temp"] = "Simulated"
  return(sim.dat)
}

#Created to speed up simResiduals.  Essentially, it implements the main loop of simResiduals in C++
#vectR: The R vector of data to "accumulate"
#ar1R: The ar1 value to use in the accumulation.
#diffTR: Difference of times, used to adjust AR1 for unequally spaced data
accumulate = cxxfunction(signature(vectR="numeric", ar1R="numeric", diffTR="numeric"), body='
  arma::vec vect = as<arma::vec>(vectR);
  arma::vec diffT = as<arma::vec>(diffTR);
  double ar1 = as<double>(ar1R);
  arma::vec output(vect.size(), 1);
  output[0] = vect[0];
  for(int i=1; i<vect.size(); ++i)
  {
    output[i] = output[i-1]*pow(ar1,diffT[i-1]) + vect[i];
  }
  return(wrap(output));
', plugin="RcppArmadillo" )

accumulateR = function(vec, ar1, diffT){
  output = rep(NA, length(vec))
  output[1] = vec[1]
  for(i in 2:length(vec)){
    output[i] = output[i-1]*ar1^diffT[i-1] + vec[i]
  }
  return(output)
}

# 5 different methods for detecting outliers:
# 1.  Take mean and sd for a fixed pressure level but all time.
# 2.  Take mean and sd for an observation using a 45 day window of obs (at a fixed pressure level)
# 3.  Take mean and sd for an observation using a 12 hour and 45 day window (at a fixed pressure level)
# 4.  Step one with values >= 6 sigma removed, then method 2.
# 5.  Step one with values >= 6 sigma removed, then method 3.
#
# Tukey's symmetric estimator vs Huber's Asymmetric estimator

#sim.dat should be the simulated data, output from simResiduals()
#type should be one of "extreme" or "negative", indicating the type of outlier to simulate
#pct should be between 0 and 1, indicating what percent of obs. should be outliers
contaminateOutlier = function(sim.dat, type, pct, muError=c(8,9,10))
{
  if(!type %in% c("extreme", "negative") )
    stop("Invalid outlier type provided.")
  if(pct>1 | pct<0)
    stop("Invalid pct provided.")
  
  rowsToSample = 1:nrow(sim.dat)
  #Don't overwrite current outliers:
  if(!is.null(sim.dat$outlierSize))
    rowsToSample = rowsToSample[sim.dat$outlierSize==0]
  contamRows = sample(1:nrow(sim.dat), size=nrow(sim.dat)*pct)
  sigma = sd(sim.dat$Simulated)
  if(type=="extreme")
  {
    #Randomly select the size of each outlier
    meanError = sample(muError, size=length(contamRows), replace=T)
    sigmaError = rnorm(length(contamRows), mean=meanError, sd=1)
  }
  if(type=="negative")
    sigmaError = rep(-2,length(contamRows))
  sim.dat[contamRows,"Simulated"] = sim.dat[contamRows,"Simulated"] +
    sigmaError*sigma*sample(c(-1,1),size=length(contamRows),replace=T)
  
  sim.dat$outlierSize = 0
  sim.dat$outlierSize[contamRows] = sigmaError
    
  return(sim.dat)
}

contaminateBreak = function(sim.dat, breaks=2, sigmaMult=c(0.2,0.4,0.6))
{
  #Force breaks to not be too close to the edge: 1 year away
  contamRows = sort( sample(365:(nrow(sim.dat)-365), size=breaks) )
  sigmaMult = sample(sigmaMult, size=breaks, replace=T)
  breakSize = rnorm(breaks, sd=sd(sim.dat$Simulated)*sigmaMult)
  
  for(i in 1:breaks)
  {
    sim.dat$Simulated[(contamRows[i]+1):nrow(sim.dat)] = 
      sim.dat$Simulated[(contamRows[i]+1):nrow(sim.dat)] + breakSize[i]
  }
  sim.dat$Breaks = 0
  sim.dat$Breaks[contamRows] = breakSize
  
  return(sim.dat)
}

########################################################################################################
##
## Functions to detect outliers in the data.  There are two flavors implemented: Tukey and Huber
## robust estimators.  There are four functions in total:
## - tukeyScore() and huberScore() are in C++ to efficiently loop over all observations and compute
##   the estimators over a window (hour and day based)
## - tukey() and huber() are base R functions that call tukeyScore() and huberScore().  The base
##   C++ functions should generally not be called directly.  If no windowing is required, these
##   base R functions check for outliers without call the *Score() functions.
##
########################################################################################################


#Assumes the data is sorted by hour group and then by date (never by year!)
#dataR: dataset to be checked for outliers
#dayR: days corresponding to data.  Used for windowing dataR.
#hourR: hours corresponding to data.  Should be a bucket value, as only values with the same hour are ever grouped together.
#dayWindowR: single value specifying the width of the window to calculate the tukey weights.
#Note: tukey scores (at least currently) are being computed across years.  Rather
# than trying to code multiple periods in C++, I solved this problem by requiring
# the dataset be sorted by hour and then day (not year!).  Thus, the first few 
# values will be all the values with hour=(first bucket), day=1, and any year.  The
# wrapper function ensures this sorting takes place, but be careful with direct
# calls to this function.
#Returns a matrix of three columns.  The first is the Tukey score for that observation,
# the second has the Tukey mean estimate, and the final has the Tukey SD estimate.
tukeyScore = cxxfunction( signature( dataR="numeric", dayR="numeric", hourR="numeric"
    ,dayWindowR="numeric"), includes="#include <cmath>", body='
  //Use vec from Armadillo because they have a median function implemented
  //Include <cmath> for abs() overloaded to accept non-int arguments
  arma::vec data = as<arma::vec>(dataR);
  arma::vec day = as<arma::vec>(dayR);
  arma::vec hour = as<arma::vec>(hourR);
  double dayWindow = as<double>(dayWindowR);

  //Use start/endWindowRow to track the rows used in computing median
  //std::nth_element will do this for you, but you need first and last element.
  //Array should exist in elements >=start and <end.
  int startWindowRow(0), endWindowRow(0); 
  while(day[endWindowRow]<day[0]+dayWindow)
    endWindowRow++;

  //Declare objects to be used in the loop
  double M;
  double MAD;
  double biweightMean;
  double biweightSD;
  double const a = 1.4826; //constant to adjust MAD for consistency
  double const c = 6.;
  //score has Tukey score, biweightMean, biweightSd
  arma::mat score(data.size(),3);

  for(int i=0; i<data.size(); i++) //Loop over the observations
  {
    //day[i] is small, EX: day[i]=3 but dayWindow=7.  Then, you need to
    //force day[startWindowRow]>=365+3-7 or day[startWindowRow]<day[i].
    if(day[i]-dayWindow<=0)
    {
      while( (day[startWindowRow]<day[i]-dayWindow+365) & (day[startWindowRow]>day[i]) )
        startWindowRow++;
    }
    else
    {
      while(day[startWindowRow]>day[i]) //Roll 365 over to 1 if day[i]>dayWindow
        startWindowRow++;
      while(day[startWindowRow]<day[i]-dayWindow)
        startWindowRow++;
    }
    //Iterate up startWindowRow if hours do not match
    while(hour[startWindowRow]!=hour[i])
      startWindowRow++;

    //day[i] is close to 365, EX: day[i]=360 and dayWindow=7.  Then, you need to
    //force day[endWindowRow]>=360+7-365 and day[endWindowRow]<day[i].
    if(day[i]+dayWindow>=365)
    {
      while((day[endWindowRow]>day[i]) & (endWindowRow<day.size()))
        endWindowRow++;
      while((day[endWindowRow]<=day[i]-365+dayWindow) & endWindowRow<day.size())
        endWindowRow++;
    }
    else
    {
      while((day[endWindowRow]<=day[i]+dayWindow) & (hour[endWindowRow]==hour[i]) & (endWindowRow<day.size()))
        endWindowRow++;
    }

    //u stores the weights for the biweight estimator
    arma::vec u(endWindowRow-startWindowRow);
    //d stores the subsetted data
    arma::vec d(endWindowRow-startWindowRow);
    
    d = data(arma::span(startWindowRow,endWindowRow-1));
    M = arma::median(d);
    //Median absolute deviation
    MAD = a*arma::median(abs(d-M));
    u = (d-M)/(c*MAD);
    for(int j=0; j<u.size(); j++)
    {
      if(u[j]>1 | u[j]<-1)
        u[j] = 1;
    } //clip u between -1 and 1
    biweightMean = M + dot((d-M),pow((1-pow(u,2.)),2.)) / dot((1-pow(u,2.)),(1-pow(u,2.)));
    biweightSD = sqrt((endWindowRow-startWindowRow)*dot( pow(d-M,2.), pow(1-pow(u,2.),4.))) /
      std::abs( dot( 1-pow(u,2.), 1-5*pow(u,2.) ) ); //Use std::abs for a non-int abs function
    score(i,0) = (data[i]-biweightMean) / biweightSD;
    score(i,1) = biweightMean;
    score(i,2) = biweightSD;
    //Below code is useful for testing purposes:
/*    if(i==1234-1)
    {
      std::cout << "Data: " << data[i] << std::endl;
      std::cout << "start: " << startWindowRow << std::endl;
      std::cout << "end: " << endWindowRow << std::endl;
      std::cout << "Median: " << M << std::endl;
      std::cout << "MAD: " << MAD << std::endl;
      std::cout << "biweightMean " << biweightMean << std::endl;
      std::cout << "biweightSD " << biweightSD << std::endl;
      std::cout << "u: ";
      for(int j=0; j<u.size(); j++)
        std::cout << u[j] << " ";
    }*/
  }

  return(wrap(score));
', plugin="RcppArmadillo" )

#Basic R implementation of Tukey's robust estimators.
#data: the data to compute the estimator on.
#c: value after which weights become 0.  c=6 is equivalent to 4 sigma in the normal case.
#Returns the robust estimator for the mean and sd
tukeyR = function(data, c=6, na.rm=T){
  if(na.rm)
    data = data[!is.na(data)]
  
  M = median(data)
  #Median absolute deviation
  MAD = 1.4826*mad(data) #1.4826 adjusts for normal consistency
  u = (data-M)/(c*MAD)
  u[abs(u)>1] = 1
  biweightMean = M + sum( (data-M)*(1-u^2)^2) / sum((1-u^2)^2)
  biweightSD = sqrt(length(data)*sum( (data-M)^2 * (1-u^2)^4)) /
    abs( sum( (1-u^2)^2*(1-5*u^2) ) )
  return(c(biweightMean=biweightMean, biweightSD=biweightSD))
}

huberR = function(data, tolerance = 1E-6, k = 1.5, na.rm=T){
  if(na.rm)
    data = data[!is.na(data)]

  h = 2*pnorm(k)-1
  beta = h+k^2*(1-h)-2*k*dnorm(k)

  yBar = median(data)
  sLower = median( abs(data-yBar) )*1.4826
  sLower0 = sLower
  sUpper = sLower
  converged = FALSE
  iterations = 0
  while(!converged)
  {
    y = data
    y[y<yBar-k*sLower] = yBar-k*sLower
    y[y>yBar+k*sUpper] = yBar+k*sUpper
  
    yBarNew = mean(y);
    #Update sLower and sUpper:
    nU = sum(y>yBarNew)
    nL = sum(y<yBarNew)
    sUpperNew = sqrt( sum( ((y-yBarNew)^2)[y>yBarNew] ) / ((nU-1)*beta) )
    sLowerNew = sqrt( sum( ((y-yBarNew)^2)[y<yBarNew] ) / ((nL-1)*beta) )
  
    #Check if values have converged:
    if( abs(yBarNew-yBar) < tolerance*sLower0
      & abs(sLowerNew-sLower) < tolerance*sLower0
      & abs(sUpperNew-sUpper) < tolerance*sLower0 )
      converged = TRUE
    iterations = iterations + 1
#      if(iterations>=99)
#        converged = TRUE

    #Set old values to new
    yBar = yBarNew;
    sLower = sLowerNew;
    sUpper = sUpperNew;
  }
  return(list(mu=yBar, sdLower = sLower, sdUpper = sUpper) )
}

#Wrapper to tukeyScore, the C++ function.
#This function checks for outliers in sim.dat using Tukey's robust estimator.
#sim.dat: the simulated dataset
#dayWindow: How many total days should be used for computing median and MAD?
#hourBuckets: How many buckets should be used to partition the hours?
# Set to 1 to skip hourBucketing.
# Buckets will be centered on 0:00 + k*24/hourBuckets
tukeyWrap = function(sim.dat, dayWindow=45, hourBuckets=2, returnEst=F)
{
  #Data checks
  if(is(sim.dat)[1]!="data.frame")
    stop("sim.dat is not a data.frame!")
  reqCols = c("Year", "Hour", "Day_Of_Year", "Simulated")
  test = reqCols %in% colnames(sim.dat)
  if( any(!test) )
    stop(paste0("Missing the following columns: ", paste(reqCols[!test],collapse=", ")))
  if( max(sim.dat$Day_Of_Year) <= 31)
    stop("Maximum value for day is <= 31.  Need day of year in Day.")
  if(dayWindow<=0)
    stop("Invalid dayWindow value.")
  
  sim.dat$HourBucket = ceiling( (sim.dat$Hour-2400/(2*hourBuckets)) / (2400/hourBuckets) )
  sim.dat$HourBucket[sim.dat$HourBucket==hourBuckets] = 0
  ord = order( sim.dat[,"HourBucket"], sim.dat[,"Day_Of_Year"] )
  sim.dat = sim.dat[ord,]
  if(dayWindow>=365)
  {
    if(dayWindow!=365)
      warning("dayWindow>365 invalid.  dayWindow=365 used instead.")
    if(hourBuckets!=1)
      warning("HourBucket ignored if no day bucketing is present")
    est = tukeyR(sim.dat$Simulated)
    biweightMean = est["biweightMean"]
    biweightSD = est["biweightSD"]
    score = (sim.dat$Simulated-biweightMean) / biweightSD
    #Unsort the score (so it corresponds with the original data)
    out = rep(0,length(score))
    out[ord] = score
    if(returnEst)
      out = cbind(out, biweightMean, biweightSD)
    return(out)
  }
  
  score = tukeyScore( sim.dat$Simulated, sim.dat$Day_Of_Year, sim.dat$HourBucket, floor( (dayWindow-1)/2 ) )
  #Unsort the score (so it corresponds with the original data)
  if(returnEst){
    out = matrix(0,nr=nrow(score), nc=ncol(score))
    out[ord,] = score
  } else {
    out = rep(0,nrow(score))
    out[ord] = score[,1]    
  }
  return(out)
}

#Assumes the data is sorted by hour group and then by date (never by year!)
#dataR: dataset to be checked for outliers
#dayR: days corresponding to data.  Used for windowing dataR.
#hourR: hours corresponding to data.  Should be a bucket value, as only values with the same hour are ever grouped together.
#dayWindowR: single value specifying the width of the window to calculate the tukey weights.
#toleranceR: single value specifying the tolerance (i.e. when the algorithm is deemed to have converged).
# Note that after 100 iterations the algorithm terminates anyways, as convergence did not occur in some cases in testing.
#betaR: Coefficient for computation of scores, see Ashley's temperature paper.
#Note: huber scores (at least currently) are being computed across years.  Rather
# than trying to code multiple periods in C++, I solved this problem by requiring
# the dataset be sorted by hour and then day (not year!).  Thus, the first few 
# values will be all the values with hour=(first bucket), day=1, and any year.  The
# wrapper function ensures this sorting takes place, but be careful with direct
# calls to this function.
huberScore = cxxfunction( signature( dataR="numeric", dayR="numeric", hourR="numeric"
    ,dayWindowR="numeric", toleranceR="numeric", betaR="numeric"), includes="#include <cmath>", body='
  //Use vec from Armadillo because they have a median function implemented
  //Include <cmath> for abs() overloaded to accept non-int arguments
  arma::vec data = as<arma::vec>(dataR);
  arma::vec day = as<arma::vec>(dayR);
  arma::vec hour = as<arma::vec>(hourR);
  double dayWindow = as<double>(dayWindowR);
  double tolerance = as<double>(toleranceR);
  double beta = as<double>(betaR);

  //Use start/endWindowRow to track the rows used in computing median
  //std::nth_element will do this for you, but you need first and last element.
  //Array should exist in elements >=start and <end.
  int startWindowRow(0), endWindowRow(0); 
  while(day[endWindowRow]<day[0]+dayWindow)
    endWindowRow++;

  //Declare objects to be used in the loop
  double yBar, yBarNew;
  double sUpper, sUpperNew;
  double sLower, sLowerNew;
  double sLowerCnt, sLowerSum;
  double sUpperCnt, sUpperSum;
  double const pi = 3.141593; //constant to adjust MAD for consistency
  double const k = 1.5;
  bool converged(false);
  //score will store the huber score, mean, left SD, and rightSD
  arma::mat score(data.size(),4);

  for(int i=0; i<data.size(); i++) //Loop over the observations
  {
    //day[i] is small, EX: day[i]=3 but dayWindow=7.  Then, you need to
    //force day[startWindowRow]>=365+3-7 or day[startWindowRow]<day[i].
    if(day[i]-dayWindow<=0)
    {
      while( (day[startWindowRow]<day[i]-dayWindow+365) & (day[startWindowRow]>day[i]) )
        startWindowRow++;
    }
    else
    {
      while(day[startWindowRow]>day[i]) //Roll 365 over to 1 if day[i]>dayWindow
        startWindowRow++;
      while(day[startWindowRow]<day[i]-dayWindow)
        startWindowRow++;
    }
    //Iterate up startWindowRow if hours do not match
    while(hour[startWindowRow]!=hour[i])
      startWindowRow++;

    //day[i] is close to 365, EX: day[i]=360 and dayWindow=7.  Then, you need to
    //force day[endWindowRow]>=360+7-365 and day[endWindowRow]<day[i].
    if(day[i]+dayWindow>=365)
    {
      while((day[endWindowRow]>day[i]) & (endWindowRow<day.size()))
        endWindowRow++;
      while((day[endWindowRow]<=day[i]-365+dayWindow) & endWindowRow<day.size())
        endWindowRow++;
    }
    else
    {
      while((day[endWindowRow]<=day[i]+dayWindow) & (hour[endWindowRow]==hour[i]) & (endWindowRow<day.size()))
        endWindowRow++;
    }

    //d stores the subsetted data
    arma::vec d(endWindowRow-startWindowRow);
    
    d = data(arma::span(startWindowRow,endWindowRow-1));
    arma::vec y(d.size()); //y stores the winsorized values
    yBar = arma::median(d);
    //Mean absolute deviation
    sLower = arma::sum(abs(d-yBar))/d.size()*2/sqrt(2*pi);
    double sLower0 = sLower; //Used for tolerance computation
    sUpper = sLower;
    converged = false;
    int iterations = 0;
    while(!converged)
    {
      iterations++;
      for(int j=0; j<y.size(); j++)
      {
        if( (d[j]>yBar-k*sLower) & (d[j]<yBar+k*sUpper) )
          y[j] = d[j];
        if( d[j] > yBar + k*sUpper )
          y[j] = yBar + k*sUpper;
        if( d[j] < yBar - k*sLower )
          y[j] = yBar - k*sLower;
      }
      yBarNew = mean(y);
      //Update sLower and sUpper:
      sLowerCnt=0; sUpperCnt=0;
      sLowerSum=0; sUpperSum=0;
      for(int j=0; j<y.size(); j++)
      {
        if(y[j]<yBarNew)
        {
          sLowerCnt++;
          sLowerSum += pow(y[j]-yBarNew, 2.);
        }
        else
        {
          sUpperCnt++;
          sUpperSum += pow(y[j]-yBarNew, 2.);
        }
      }
      sUpperNew = sqrt( sUpperSum/(beta*(sUpperCnt-1)) );
      sLowerNew = sqrt( sLowerSum/(beta*(sLowerCnt-1)) );

      //Check if values have converged:
      if( std::abs(yBarNew-yBar) < tolerance*sLower0
        & std::abs(sLowerNew-sLower) < tolerance*sLower0
        & std::abs(sUpperNew-sUpper) < tolerance*sLower0 )
        converged = true;

      //Go to next if taking too long
      if(iterations>=99)
        converged=true;

      //Set old values to new
      yBar = yBarNew;
      sLower = sLowerNew;
      sUpper = sUpperNew;

    }

    if(data[i] < yBar)
      score(i,0) = (data[i]-yBar) / sLower;
    else
      score(i,0) = (data[i]-yBar) / sUpper;
    score(i,1) = yBar;
    score(i,2) = sLower;
    score(i,3) = sUpper;
  }

  return(wrap(score));
', plugin="RcppArmadillo" )


#Wrapper to huberScore, the C++ function.
#This function checks for outliers in sim.dat using Huber's robust asymmetric estimator.
#sim.dat: the simulated dataset
#dayWindow: How many total days should be used for computing the estimators?
#hourBuckets: How many buckets should be used to partition the hours?
# Set to 1 to skip hour bucketing.
# Buckets will be centered on 0:00 + k*24/hourBuckets
huber = function(sim.dat, dayWindow=45, hourBuckets=2, returnEst=F)
{
  #Data checks
  if(is(sim.dat)[1]!="data.frame")
    stop("sim.dat is not a data.frame!")
  reqCols = c("Year", "Hour", "Day_Of_Year", "Simulated")
  test = reqCols %in% colnames(sim.dat)
  if( any(!test) )
    stop(paste0("Missing the following columns: ", paste(reqCols[!test],collapse=", ")))
  if( max(sim.dat$Day_Of_Year) <= 31)
    stop("Maximum value for day is <= 31.  Need day of year in Day.")
  if(dayWindow<=0)
    stop("Invalid dayWindow value.")
  
  sim.dat$HourBucket = ceiling( (sim.dat$Hour-2400/(2*hourBuckets)) / (2400/hourBuckets) )
  sim.dat$HourBucket[sim.dat$HourBucket==hourBuckets] = 0
  #Save ordering
  ord = order( sim.dat[,"HourBucket"], sim.dat[,"Day_Of_Year"] )
  sim.dat = sim.dat[ord,]

  tolerance = 1E-6
  k = 1.5
  h = 2*pnorm(k)-1
  beta = h+k^2*(1-h)-2*k*dnorm(k)

  if(dayWindow>=365)
  {
    if(dayWindow!=365)
      warning("dayWindow>365 invalid.  dayWindow=365 used instead.")
    if(hourBuckets!=1)
      warning("HourBucket ignored if no day bucketing is present")
    yBar = median(sim.dat$Simulated)
    sLower = median( abs(sim.dat$Simulated-yBar) )*1.4826
    sLower0 = sLower
    sUpper = sLower
    converged = FALSE
    iterations = 0
    while(!converged)
    {
      y = sim.dat$Simulated
      y[y<yBar-k*sLower] = yBar-k*sLower
      y[y>yBar+k*sUpper] = yBar+k*sUpper
      
      yBarNew = mean(y);
      #Update sLower and sUpper:
      nU = sum(y>yBarNew)
      nL = sum(y<yBarNew)
      sUpperNew = sqrt( sum( ((y-yBarNew)^2)[y>yBarNew] ) / ((nU-1)*beta) )
      sLowerNew = sqrt( sum( ((y-yBarNew)^2)[y<yBarNew] ) / ((nL-1)*beta) )
    
      #Check if values have converged:
      if( abs(yBarNew-yBar) < tolerance*sLower0
        & abs(sLowerNew-sLower) < tolerance*sLower0
        & abs(sUpperNew-sUpper) < tolerance*sLower0 )
        converged = TRUE
      iterations = iterations + 1
#      if(iterations>=99)
#        converged = TRUE
    
      #Set old values to new
      yBar = yBarNew;
      sLower = sLowerNew;
      sUpper = sUpperNew;
    }
    
    score = ifelse( sim.dat$Simulated<yBar
                  ,(sim.dat$Simulated-yBar)/sLower
                  ,(sim.dat$Simulated-yBar)/sUpper )
    #Unsort the scores (so they're ordered in the same way as the original data)
    out = rep(0,length(score))
    out[ord] = score
    if(returnEst)
      out = cbind(out, yBar, sLower, sUpper)
    return(out)
  }

  score = huberScore(sim.dat$Simulated, sim.dat$Day_Of_Year, sim.dat$HourBucket, dayWindow, tolerance, beta)
  #Unsort the scores (so they're ordered in the same way as the original data)
  if(returnEst){
    out = matrix(0,nr=nrow(score), nc=ncol(score))
    out[ord,] = score
  } else {
    out = rep(0,nrow(score))
    out[ord] = score[,1]    
  }
  return(out)
}

rolling = function(data, day, hour, dayWindow, estFunc){
  if(!is(estFunc,"function")){
    stopifnot(estFunc %in% c("huber", "tukey"))
    if(estFunc=="huber")
      estFunc=huberR
    else #if(estFunc=="tukey")
      estFunc=function(data){
        est = tukeyR(data)
        return(list(mean=est$biweightMean, sdUpper=est$biweightSD, sdLower=est$biweightSD) )
      }
  }
    
  out = data.frame(data=data)
  out$mean = NA
  out$sdUpper = NA
  out$sdLower = NA
  date = as.Date( paste0("2000-",day), format="%Y-%j")
  dateUnique = unique(date)
  for(i in 1:length(dateUnique)){
    i = dateUnique[i]
    for(j in unique(hour)){
      filter = rep(FALSE, length(day))
  	  if(i-as.Date("1998-01-01")<dayWindow)
  	    filter[date>i-dayWindow+366] = TRUE
  	  if(i+dayWindow>as.Date("1999-01-01"))
  	    filter[date<i+dayWindow-366] = TRUE
  	  filter[abs(i-date)<dayWindow] = TRUE
  	  filter[hour!=j] = FALSE
  	  est = estFunc(data[filter])
  	  out$mean[date==i & hour==j] = est$mu
  	  out$sdUpper[date==i & hour==j] = est$sdUpper
  	  out$sdLower[date==i & hour==j] = est$sdLower
	  }
  }
  out$statistic = ifelse(out$data<out$mean
                        ,(out$mean-out$data)/out$sdLower
                        ,(out$data-out$mean)/out$sdUpper )
  return(out)
}

########################################################################################################
##
## Functions to check for breaks in the data.  A "break" is a systematic error.  This code looks for
## changes in the mean value of a time series and removes those jumps by adding a constant.  Most of the
## work is done by homogenizationCpp, and homogenization is a R function that's basically a wrapper to
## the C++ version.
##
########################################################################################################

#dataR: vector of the data to be homogenized
#dayR: vector of the days corresponding to the times
#periodR: the time period to use as a windowing period for the homogenization model.
homogenizationCpp = cxxfunction(signature(dataR="numeric", dayR="numeric", periodR="numeric")
  ,body="
  //Initializations
  arma::vec d = as<arma::vec>(dataR);
  arma::vec day = as<arma::vec>(dayR);
  int period = as<int>(periodR);
  double leftSum(0), rightSum(0), leftMean(0), rightMean(0), sumX2(0);
  double mean, sigma2, testStat;
  arma::mat output(d.n_rows,3, arma::fill::zeros);

  //center is the row of the current obs we're computing the SNHT for.
  //leftBound is the first row in the window, i.e. the first row s.t. day(leftBound)>day(center)-period
  //rightBound is the last row in the window, i.e. the last row s.t. day(rightBound)<day(center)+period
  int rightBound(0), leftBound(0), center(0);
  while(day(center+1)-day(0)<period)
  {
    leftSum += d(center);
    sumX2 += pow(d(center),2.);
    center++;
  }
  rightBound = center+1;
  while(day(rightBound)<day(center)+period)
  {
    rightSum += d(rightBound);
    sumX2 += pow(d(rightBound),2.);
    rightBound++;
  }
  
  //Now, increment center until we no longer have a full window (i.e. rightBound=d.n_rows)
  while(rightBound<(d.n_rows-1))
  {
    //Adjust leftBound as necessary, and update leftMean and sumX2 along the way
    while(day(leftBound) < (day(center) - period))
    {
      leftSum -= d(leftBound);
      sumX2 -= pow(d(leftBound),2.);
      leftBound++;
    }
    //Adjust rightBound as necessary, and update rightMean and sumX2 along the way
    while( day(rightBound+1) < (day(center) + period) )
    {
      rightSum += d(rightBound+1);
      sumX2 += pow(d(rightBound+1),2.);
      rightBound++;
      if(rightBound == (d.n_rows-1))
        break; //quit so that we don't attempt day(d.nrows-1+1), which would fail
    }

    //Recalculate test statistic
    //rightBound-leftBound is number of elements between rightBound and leftBound,
    //inclusive, minus one.  We don't include center, so this is perfect.
    mean = (leftSum + rightSum)/(rightBound-leftBound);
    leftMean = leftSum/(center-leftBound-1);
    rightMean = rightSum/(rightBound-center-1);
    sigma2 = sumX2/(rightBound-leftBound)-pow(mean,2.); //E(X^2)-E(X)^2;
    testStat = (double)period/(sqrt(sigma2)) * ( pow(leftMean-mean,2.) + pow(rightMean-mean,2.) );
    output(center,0) = testStat;
    output(center,1) = leftMean;
    output(center,2) = rightMean;

    leftSum += d(center);
    center++;
    rightSum -= d(center);
  }

  return(wrap(output));
", plugin="RcppArmadillo")

#data: the time-series data to homogenize
#type: What kind of homogeneity test should be performed?  BinSeg and PELT both come from
#     the changepoint package.
#period: the length of the window period.  A period of one year means that for each point, we'll 
#     compare the previous year and next year to see if a jump occured.
#crit: The critical value for the SNHT.  Monte Carlo suggests 15, Haimberger suggests 100 in
#     "On the homogeneity of radiosonde wind time series."  Only used for SNHT, robust.
#minHomWidth: After a break has been detected, how close can other breaks be?  Only used for
#     SHNT, robust.
#scaled: If TRUE, divide the statistic in Haimberger by s^2 instead of s.
#Test data:
#  sim.dat = data.frame(Simulated=rnorm(365*4)+rep(0:1,each=365*2), Day_Of_Year=rep(1:365,each=4)
#     ,Year=2011+rep(0:3,each=365), Hour=2300 )
#  sim.dat$Date = as.POSIXct( paste(sim.dat$Year, sim.dat$Day_Of_Year, floor(sim.dat$Hour/100), round((sim.dat$Hour %% 100)*60/100,0) ),format="%Y %j %H %M")
#  sim.dat$Outlier.Fl = rbinom(nrow(sim.dat), size=1, prob=.01)
#benjamini: Should the test statistics be adjusted via benjamini?  Only for SNHT, robust
#adj: What value should the test statistics be divided by to correct for autocorrelation? Only for SNHT, robust
#alpha: If using the benjamini adjustment, this dictates the alpha level at which to
#  reject the null hypthesis.
homogenization = function(sim.dat, period, type=c("SNHT", "robust", "BinSeg", "PELT")
        ,crit=100, rmPeriod=F, rmTrend=F, rmAC=F, scaled=F, benjamini=T, adj=1, alpha=.05
        ,rmSeasonalPeriod=365, ...)
{
  if(is(sim.dat)[1]!="data.frame")
    stop("sim.dat is not a data.frame!")
  reqCols = c("Year", "Day_Of_Year", "Hour", "Date", "Simulated")
  test = reqCols %in% colnames(sim.dat)
  if( any(!test) )
    stop(paste0("Missing the following columns: ", paste(reqCols[!test],collapse=", ")))
  if( max(sim.dat$Hour)<25 )
    stop("Largest sim.dat$Hour is <25, expect sim.dat$Hour to be of the form HHMM")
  if( max(sim.dat$Day_Of_Year)<32 )
    stop("Largest sim.dat$Day_Of_Year is <32, expect sim.dat$Day_Of_Year to be in 1 to 365")
  if(length(type)>1)
    type = type[1]
  if(!type %in% c("SNHT", "robust", "BinSeg", "PELT") )
    stop("Invalid type specified!")
  
  #Require homogenizations to be at least period observations apart:
  minHomWidth=period
  
  #Remove period, if desired
  if(rmPeriod | rmTrend){
    if(rmPeriod & rmTrend)
      fit = gam( Simulated ~ s(Day_Of_Year) + Year, data=sim.dat )
    if(rmPeriod & !rmTrend)
      fit = gam( Simulated ~ s(Day_Of_Year), data=sim.dat )
    if(!rmPeriod & rmTrend)
      fit = lm( Simulated ~ Year, data=sim.dat )
    modelAdj = predict(fit, newdata=sim.dat)
    sim.dat$Simulated = as.numeric( sim.dat$Simulated - modelAdj )
  }

  #Remove autocorrelation, if desired
  if(rmAC){
    warning("No method currently implemented to remove autocorrelation.")
  }

  #Order sim.dat by time:
  ord = order(sim.dat$Year*365 + sim.dat$Day_Of_Year + sim.dat$Hour/2400)
  sim.dat = sim.dat[ord,]

  #Only use non-outliers for the homogenization (if outliers have been labeled)
  if("outlierFl" %in% colnames(sim.dat)){
    sim.dat.orig = sim.dat
    sim.dat = sim.dat[!sim.dat$outlierFl,]
  }
  
  #SNHT and robust work similarly, BinSeg and PELT need an entirely different implementation
  #This code block computes the breaks and adjusts sim.dat and sim.dat.orig appropriately
  if(type %in% c("SNHT", "robust")){
    if(type=="SNHT")
#       scores = homogenizationCpp(sim.dat$Simulated
#         ,sim.dat$Year*365 + sim.dat$Day_Of_Year + sim.dat$Hour/2400 - 365*1900, period)
      scores = snht(data=sim.dat$Simulated, period, scaled=scaled, robust=FALSE
                ,rmSeasonalPeriod = rmSeasonalPeriod, time=as.numeric(sim.dat$Date)/(24*60*60))
    if(type=="robust"){
      scores = snht(data=sim.dat$Simulated, period, scaled=scaled, robust=TRUE
                ,rmSeasonalPeriod = rmSeasonalPeriod, time=as.numeric(sim.dat$Date)/(24*60*60))
      scores$time = NULL
    }
    colnames(scores) = c("score", "leftMean", "rightMean")
    scores[, 1] = scores[, 1]/adj
    if(benjamini){
      #The default value of crit=100 is fairly arbitrary.  We instead use the Benjamini
      #approach here, and so we set crit=.0001 and manually set insignificant points to
      #values of 0.
      crit = .0001
      pvals = 1-pchisq(scores[,1], df=1)
      reject = sort(pvals) < alpha/(length(sort(pvals)):1)
      minp = max(sort(pvals)[reject])
      scores[pvals>=minp & !is.na(pvals), 1] = 0
    }
    while(max(scores[, 1], na.rm = TRUE)>crit)
    {
      #Which row is most likely to be the break?
      maxScoreRow = which.max(scores[,1])
      diff = scores[maxScoreRow,"leftMean"] - scores[maxScoreRow,"rightMean"]
      
      #Adjust data using the break
      sim.dat$Simulated[maxScoreRow:nrow(sim.dat)] =
        sim.dat$Simulated[maxScoreRow:nrow(sim.dat)] + diff
      
      #Adjust dataset with outliers as well, if applicable
      if("outlierFl" %in% colnames(sim.dat))
      {
        filt = sim.dat.orig$Date > sim.dat$Date[maxScoreRow] |
              (sim.dat.orig$Hour >= sim.dat$Hour[maxScoreRow] &
               sim.dat.orig$Date == sim.dat$Date[maxScoreRow])
        sim.dat.orig$Simulated[filt] = sim.dat.orig$Simulated[filt] + diff
      }
      toBind = cbind(sim.dat[maxScoreRow,c("Year","Day_Of_Year","Hour","Date")],difference=diff)
      if(exists("breaksFound")){
        breaksFound = rbind( breaksFound, toBind )
      } else {
        breaksFound = data.frame(toBind)
        colnames(breaksFound) = c("Year", "Day_Of_Year", "Hour", "Date", "difference")
      }
      
      #Set scores to 0 if they are within minHomWidth of a current break 
      filt = lapply(breaksFound$Date, function(x){
          return(as.Date(x, origin="MST") < as.Date(sim.dat$Date, origin="MST") + minHomWidth & 
                 as.Date(x, origin="MST") > as.Date(sim.dat$Date, origin="MST") - minHomWidth)
        })
      filt = apply(do.call("cbind", filt),1,max)==1
      scores[filt,] = 0
    }
    if(!exists("breaksFound"))
      breaksFound = NULL
  }
  
  if(type %in% c("BinSeg","PELT")){
    brkRows = cpt.mean( sim.dat$Simulated, method=type, ... )@cpts
    brkRows = brkRows[brkRows>1]
    #Create a variable to split sim.dat by and compute means:
    grps = sapply(1:nrow(sim.dat), function(i){sum(brkRows<=i)+1} )
    means = as.numeric(by( sim.dat$Simulated, grps, FUN=mean ))
    breaksFound = cbind(sim.dat[brkRows,c("Year","Day_Of_Year","Hour","Date")],difference=diff(means))

    #Adjust all means to first mean
    sim.dat$Simulated = sim.dat$Simulated - means[grps] + means[1]

    #Also adjust sim.dat.orig
    if("outlierFl" %in% colnames(sim.dat)){
      grps.orig = sapply(sim.dat.orig$Date, function(i){sum(breaksFound$Date<=i)+1})
      sim.dat.orig$Simulated = sim.dat.orig$Simulated - means[grps.orig] + means[1]
    }
  }
  
  #Add back in seasonality and trend, if it was removed
  if(rmPeriod | rmTrend){
    if("outlierFl" %in% colnames(sim.dat)){
      sim.dat$Simulated = sim.dat$Simulated + modelAdj[!sim.dat.orig$outlierFl]
      sim.dat.orig$Simulated = sim.dat.orig$Simulated
    } else {
      sim.dat$Simulated = sim.dat$Simulated + modelAdj
    }
  }
  
  #Add back in autocorrelation, if it was removed
  if(rmAC){warning("No autocorrelation model currently implemented")}

  #Return sim.dat ordered by original ordering:
  if("outlierFl" %in% colnames(sim.dat))
  {
    sim.dat.orig[ord,] = sim.dat.orig
    return(list(sim.dat=sim.dat.orig, breaksFound=breaksFound))
  } else {
    sim.dat[ord,] = sim.dat
    return(list(sim.dat=sim.dat, breaksFound=breaksFound))
  }
}