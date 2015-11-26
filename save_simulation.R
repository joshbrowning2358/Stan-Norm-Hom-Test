library(ggplot2)
library(data.table)

## This file if very similar to run_simulation, but it is instead designed to
## save 10 simulations for further examination. 

# cArgs = commandArgs(trailingOnly=TRUE)
cArgs = c("35121", "850", "40")

if(length(cArgs)!=3)
  stop("Exactly three args are required: station, pressure and years!")
if(!as.numeric(cArgs[2]) %in% c(100,300,850))
  stop("Pressure level must be in 100, 300, 850!")
if(!as.numeric(cArgs[3]) %in% c(10,20,40))
  stop("Years must be in 10,20,40!")

station = as.numeric(cArgs[1])
dataset = "trh" #allowable values are "trh", "trhc", "wind", "windc"
pressure = as.numeric(cArgs[2])
startYear = 2000
endYear = startYear + as.numeric(cArgs[3]) - 1
cat("Using station",station,"\n")
cat("Using pressure",pressure,"\n")
cat("Using # of years", as.numeric(cArgs[3]))
runId = paste(cArgs, collapse="_")

if(Sys.info()[1]=="Windows" & Sys.info()[4]=="JOSH_LAPTOP"){
  setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Wind QC")
  source("Code/simulation_functions.R")
  library(snht)
  removeSeasonalPeriod = snht:::removeSeasonalPeriod
}
if(Sys.info()[1]=="Linux" & grepl("(ch120|bb136)",Sys.info()[4]) ){
  setwd("~/Research/Wind_QC")
  source("Code/simulation_functions.R")
  source("~/Github/Stan-Norm-Hom-Test/snht/R/robustSNHT.R")
  source("~/Github/Stan-Norm-Hom-Test/snht/R/robustSNHTunequal.R")
  source("~/Github/Stan-Norm-Hom-Test/snht/R/snht.R")
}
if(Sys.info()[1]=="Linux" & Sys.info()[4]=="jb" ){
  setwd("/media/storage/Professional Files/Mines/Research/Wind QC")
  source("Code/simulation_functions.R")
  source("/media/storage/Github/Stan-Norm-Hom-Test/snht/R/robustSNHT.R")
  source("/media/storage/Github/Stan-Norm-Hom-Test/snht/R/robustSNHTunequal.R")
  source("/media/storage/Github/Stan-Norm-Hom-Test/snht/R/snht.R")
}

test.data = read.csv(file=paste("Data/uadb", dataset, station, "parsed_temp_cleaned.csv", sep="_") )
test.data = test.data[test.data$Pressure==pressure,]
test.data$Pressure = NULL

#Load Temperature Data
colnames(test.data)=c("Year","Month","Day","Hour","Reading")
test.data$Date = as.Date( paste(test.data$Year, test.data$Month, test.data$Day), format="%Y %m %d")
test.data$Day_Of_Year = as.numeric( as.character( test.data$Date, format="%j" ) )
#Replace -999 with NA for temperature:
test.data$Reading[test.data$Reading==-999] = NA

#Fit time series and residuals
fitTS = fitBaseTS(data=test.data)
fitST = getSkewTParams(fitTS, test.data)
test.data$Preds = predict(fitTS, newdata=test.data)
test.data$Err = test.data$Reading - test.data$Preds
ar0.5 = acfUnequal(test.data, data.col=which(colnames(test.data)=="Err"), lagScale=.5)[1,2]

#Bounds for seed are +/-2147483647 (at least on my machine)
seeds = round(runif(3, min=-21474836, max=21474836))

print("Beginning model building process...")

start = Sys.time()
for(i in 1:length(seeds) )
{
  seed = seeds[i]
  set.seed(seed)
  
  #Simulate new dataset
  sim.dat = simBaseTS(fitTS, data=test.data, startYear=startYear, endYear=endYear)
  sim.dat = simResiduals2(fitST, sim.dat, ar1=ar0.5^2)
  orig = sim.dat
  
  #Contaminate with outliers and breaks
  sim.dat = contaminateOutlier(sim.dat, "extreme", sample(c(0, 1,2,5,10)/100,size=1), muError=c(8,9,10))
  sim.dat = contaminateBreak(sim.dat, breaks=sample(1:3*as.numeric(cArgs[3])/10,size=1), sigmaMult=c(.2,.4,.6))
  contam = sim.dat
  #Set the outlierFl now, and adjust it once outlier detection algo is run
  #We do this because we need this var for the homogenization, but we don't
  #always run the outlier detection before running the homogenization.
  contam$outlierFl = F
  
  outlierType = 5
  #Type 1: Global robust mean/sd ("Tier 1")
  #Type 2: Robust mean/sd using day windows ("Tier 2")
  #Type 3: Robust mean/sd using day windows and hourly bins ("Hourly Tier 2")
  #Type 4: Type 1 and then Type 2 ("Combined")
  #Type 5: Type 1 and then Type 3 ("Hourly Combined")
  homogType = 2
  #Type 1: SNHT
  #Type 2: Robust SNHT
  #Type 3: PELT (implemented in package changepoint)
  #Type 4: BinSeg (implemented in package changepoint)
  #Type 5: SNHT (Trend and seasonality removed)
  #Type 6: Robust SNHT (Trend and seasonality removed)
  #Type 7: PELT (Trend and seasonality removed)
  #Type 8: BinSeg (Trend and seasonality removed)
  #Future Types:
  #Felix's implementation, from the conference
  for(order in c("hFirst", "oFirst"))
  {
    for(qcIterations in 2:3 )
    {
      #for(estimator in c("tukey", "huber"))
      estimator = "huber"
      #Reset sim.dat to the contaminated data, otherwise labels and
      #corrections from the homogenization and outlier detection will be
      #in place.
      sim.dat = contam
      #Also remove currently computed breaks:
      rm(breaks)
      #Use testNo to loop over the applications of the tests.  First, we'll apply
      #either the homogenization test or the outlier test.  Then, in the second
      #iteration, we'll apply the next test.  So, if testNo==1 and order==hFirst
      #then homogenization, if testNo==1 and order==oFirst then outlier, if
      #testNo==2 and order==hFirst the outlier, etc.
      for(testNo in 1:qcIterations)
      {
        
        #Homogenization Test
        if((order=="hFirst" & testNo%%2==1) | (order=="oFirst" & testNo%%2==0))
        {
          if(homogType==1) out = homogenization(sim.dat, period=2*365, minHomWidth=90)
          if(homogType==2) out = homogenization(sim.dat = contam,
                                                      period = 365,
                                                      type = "robust",
                                                      rmPeriod = FALSE,
                                                      rmTrend = FALSE,
                                                      scaled = T,
                                                      benjamini = T,
                                                      alpha = 0.01, adj = 1)
          if(homogType==3) out = homogenization( sim.dat, period=2*365, type="BinSeg" )
          if(homogType==4) out = homogenization( sim.dat, period=2*365, type="PELT" )
          if(homogType==5) out = homogenization( sim.dat, period=2*365, minHomWidth=90, rmPeriod=T, rmTrend=T )
          if(homogType==6) out = homogenization( sim.dat, period=2*365, minHomWidth=90, type="robust", rmPeriod=T, rmTrend=T )
          if(homogType==7) out = homogenization( sim.dat, period=2*365, type="BinSeg", rmPeriod=T, rmTrend=T )
          if(homogType==8) out = homogenization( sim.dat, period=2*365, type="PELT", rmPeriod=T, rmTrend=T )
          ## Save outlierFl so it isn't overwritten
          outlierFl = sim.dat$outlierFl
          sim.dat = out[[1]]
          sim.dat$outlierFl = outlierFl
          if(is.null(out[[2]])){
            out[[2]] = matrix(0, nrow=0, ncol=5)
            colnames(out[[2]]) = c("Year", "Day_Of_Year", "Hour", "Date", "difference")
          }
          #breaks will already exist if this is the second homogenization.
          #Use code below to store existing breaks as well as new ones.
          if(exists("breaks")){
            breaks = rbind(breaks, data.frame(out[[2]]))
          } else {
            breaks = data.frame(out[[2]])
          }
        }
          
        #Outlier Test
        if((order=="oFirst" & testNo%%2==1) | (order=="hFirst" & testNo%%2==0))
        {
          #(Re)compute the outlierFl.  Set to F to clear any previous detections
          sim.dat$outlierFl = F
          #types 1, 4, and 5 all have "Tier 1"
          if(outlierType %in% c(1, 4, 5))
          {
            if(estimator=="tukey"){
              score = tukey(sim.dat, dayWindow=365, hourBuckets=1)
            } else {
              score = huber(sim.dat, dayWindow=365, hourBuckets=1)
            }
            sim.dat$outlierFl[abs(score)>6] = T
          }
          #type 2 and 4 have "Tier 2"
          if(outlierType %in% c(2, 4))
          {
            if(estimator=="tukey"){
              score = tukey(sim.dat, dayWindow=45, hourBuckets=1)
            } else {
              score = huber(sim.dat, dayWindow=45, hourBuckets=1)
            }
            sim.dat$outlierFl[abs(score)>5] = T
          }
          #type 3 and 5 have "Hourly Tier 2"
          if(outlierType %in% c(3, 5))
          {
            if(estimator=="tukey"){
              score = tukey(sim.dat, dayWindow=45, hourBuckets=2)
            } else {
              score = huber(sim.dat, dayWindow=45, hourBuckets=2)
            }
            sim.dat$outlierFl[abs(score)>5] = T
          }
        }
      } #Close testNo for loop
      save(sim.dat, orig, contam, file = paste0("simulation_sample_iter.",
           qcIterations, "_order.", order, "_ID.", seed, "_n.",
           endYear - startYear, ".RData"))
    } #Close qcIterations loop
  } #Close order for loop
} #Close seed for loop
Sys.time()-start





files = dir(getwd(), pattern = ".RData")

vals = gsub(".*ID.|.RData", "", files)
nVals = stringr::str_extract(string = vals, pattern = "[-0-9]*$")
IDVals = stringr::str_extract(string = vals, pattern = "^[-0-9]*")
id = unique(cbind(nVals, IDVals))
for(i in 1:nrow(id)){
    n = id[i, "nVals"]
    ID = id[i, "IDVals"]
    
    load(paste0(paste0("simulation_sample_iter.",
                 2, "_order.", "hFirst", "_ID.", ID, "_n.", n, ".RData")))
    sr = as.data.table(sim.dat)
    sr[, time := as.POSIXct(paste(Year,
                formatC(Day_Of_Year, width = 3, flag = "0"),
                formatC(Hour, width = 4, flag = "0")),
                                format = "%Y %j %H%M")]
    sr[, type := "SR"]
    load(paste0(paste0("simulation_sample_iter.",
                 2, "_order.", "oFirst", "_ID.", ID, "_n.", n, ".RData")))
    rs = data.table(sim.dat)
    rs[, time := as.POSIXct(paste(Year,
                formatC(Day_Of_Year, width = 3, flag = "0"),
                formatC(Hour, width = 4, flag = "0")),
                                format = "%Y %j %H%M")]
    rs[, type := "RS"]
    load(paste0(paste0("simulation_sample_iter.",
                 3, "_order.", "hFirst", "_ID.", ID, "_n.", n, ".RData")))
    srs = data.table(sim.dat)
    srs[, time := as.POSIXct(paste(Year,
                formatC(Day_Of_Year, width = 3, flag = "0"),
                formatC(Hour, width = 4, flag = "0")),
                                format = "%Y %j %H%M")]
    srs[, type := "SRS"]
    load(paste0(paste0("simulation_sample_iter.",
                 3, "_order.", "oFirst", "_ID.", ID, "_n.", n, ".RData")))
    rsr = data.table(sim.dat)
    rsr[, time := as.POSIXct(paste(Year,
                formatC(Day_Of_Year, width = 3, flag = "0"),
                formatC(Hour, width = 4, flag = "0")),
                                format = "%Y %j %H%M")]
    rsr[, type := "RSR"]
    dat = rbind(sr, rs, srs, rsr)
    
    ggplot(dat, aes(x = time, y = Simulated)) + geom_point(aes(color = outlierFl)) +
        facet_wrap( ~ type)
    
    dat[, trueOutlier := outlierSize > 0]
    dat[, .N, by = c("outlierFl", "trueOutlier", "type")]
}