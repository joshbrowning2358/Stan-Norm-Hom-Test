if(Sys.info()[1]=="Windows" & Sys.info()[4]=="JOSH_LAPTOP"){
  # Only run this code if it hasn't been run already:
  if(getwd() != "C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Wind QC"){
    setwd("C:/Users/rockc_000/Documents/Professional Files/Mines/Research/Wind QC")
    source("Code/simulation_functions.R")
    library(snht)
    removeSeasonalPeriod = snht:::removeSeasonalPeriod
  }
}
if(Sys.info()[1]=="Linux" & grepl("(ch120|bb136)",Sys.info()[4]) ){
  # Only run this code if it hasn't been run already:
  if(getwd() != "~/Research/Wind_QC"){
    setwd("~/Research/Wind_QC")
    source("Code/simulation_functions.R")
    source("~/Github/Stan-Norm-Hom-Test/snht/R/robustSNHT.R")
    source("~/Github/Stan-Norm-Hom-Test/snht/R/robustSNHTunequal.R")
    source("~/Github/Stan-Norm-Hom-Test/snht/R/snht.R")
  }
}
if(Sys.info()[1]=="Linux" & Sys.info()[4]=="jb" ){
  if(getwd() != "/media/storage/Professional Files/Mines/Research/Wind QC"){
    setwd("/media/storage/Professional Files/Mines/Research/Wind QC")
    source("Code/simulation_functions.R")
    source("/media/storage/Github/Stan-Norm-Hom-Test/snht/R/robustSNHT.R")
    source("/media/storage/Github/Stan-Norm-Hom-Test/snht/R/robustSNHTunequal.R")
    source("/media/storage/Github/Stan-Norm-Hom-Test/snht/R/snht.R")
  }
}
if(Sys.info()[1]=="Linux" & Sys.info()[4]=="joshua-Ubuntu-Linux"){
  # Only run this if code hasn't been run already
  if(getwd() != "~/Documents/Github/Stan-Norm-Hom-Test/"){
    setwd("~/Documents/Github/Stan-Norm-Hom-Test/")
    source("simulation_functions.R")
    library(snht)
    removeSeasonalPeriod = snht:::removeSeasonalPeriod
  }
  suppressPackageStartupMessages(library(doParallel))
  library(foreach)
  registerDoParallel(cores=detectCores(all.tests=TRUE))
}

# cArgs = commandArgs(trailingOnly=TRUE)
completedRuns = dir(getwd(), pattern = "Simulations.*RData")
params = gsub(".*(JOSH_LAPTOP|jb|ch120|bb136|joshua-Ubuntu-Linux)_", "", completedRuns)
params = gsub(".RData", "", params)
params = data.frame(do.call("rbind", strsplit(params, split = "_")))
colnames(params) = c("station", "pressure", "n")
n = gsub("^[A-Za-z_]*", "", completedRuns)
n = as.numeric(gsub("_.*", "", n))
completedCounts = aggregate(n, by = params, FUN = sum)
allRows = expand.grid(c("35121", "74794", "82332", "85543", "94150", "51777",
                        "70133", "70261", "72387", "72456"),
                      c("100", "300", "850"), c("10", "20", "40"))
colnames(allRows) = colnames(params)
counts = merge(allRows, completedCounts, by = colnames(allRows), all.x = TRUE)
counts$x[is.na(counts$x)] = 0
counts$station = as.character(counts$station)
counts$pressure = as.character(counts$pressure)
counts$n = as.character(counts$n)
stillToRun = (1:90)[counts$x < 500]
cArgs = counts[sample(stillToRun, size = 1), 1:3]

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
cat("Using # of years", as.numeric(cArgs[3]), "\n")
runId = paste(cArgs, collapse="_")

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

## Determine completed files and start from there
# files = dir(".", pattern = paste0("Simulations_.*_", runId, "\\.RData"))
# simNumber = as.numeric(gsub("(Simulations_|_ch.*|_bb.*)", "", files))
# completedSims = 0
# if(length(simNumber) >= 1){
#     bestFile = files[grep(paste0(max(simNumber), "_(ch|bb)"), files)]
#     load(bestFile)
#     completedSims = nrow(dataStat)
# }

#Bounds for seed are +/-2147483647 (at least on my machine)
# seeds = round(runif(500 - completedSims, min=-21474836, max=21474836))
# Run 50 simulations:
seeds = sample(-21474836:21474836, size = 100, replace = FALSE)

cat("Beginning model building process...\n")

start = Sys.time()
results = foreach(i = 1:length(seeds)) %dopar%
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
  
  #Return statistics about the simulated dataset:
  # toBind = data.frame( ID=seed
  # If in parallel, need to rm detectionStat
  rm(detectionStat)
  dataStat = data.frame(ID = seed
    ,outlierPct  = mean( as.numeric( sim.dat$outlierSize!=0 ) )
    ,outlierPct4 = mean( as.numeric( sim.dat$outlierSize==4 ) ) #% of outliers at 4 sigma
    ,outlierPct5 = mean( as.numeric( sim.dat$outlierSize==5 ) ) #% of outliers at 5 sigma
    ,outlierPct6 = mean( as.numeric( sim.dat$outlierSize==6 ) ) #% of outliers at 6 sigma
    ,avgTimeBtwnBreaks = mean( diff( (1:nrow(sim.dat))[sim.dat$Breaks!=0] ) )
    ,simBr1Size = sim.dat$Breaks[sim.dat$Breaks!=0][1]
    ,simBr2Size = sim.dat$Breaks[sim.dat$Breaks!=0][2]
    ,simBr3Size = sim.dat$Breaks[sim.dat$Breaks!=0][3]
    ,simBr1Time = as.Date(paste(sim.dat$Year,sim.dat$Day_Of_Year),"%Y %j")[sim.dat$Breaks!=0][1]
    ,simBr2Time = as.Date(paste(sim.dat$Year,sim.dat$Day_Of_Year),"%Y %j")[sim.dat$Breaks!=0][2]
    ,simBr3Time = as.Date(paste(sim.dat$Year,sim.dat$Day_Of_Year),"%Y %j")[sim.dat$Breaks!=0][3]
    ,n = nrow(sim.dat)
  )
#   if(exists("dataStat"))
#     dataStat = rbind(dataStat, toBind)
#   else
#     dataStat = toBind
  
  for(outlierType in 5)
  #Type 1: Global robust mean/sd ("Tier 1")
  #Type 2: Robust mean/sd using day windows ("Tier 2")
  #Type 3: Robust mean/sd using day windows and hourly bins ("Hourly Tier 2")
  #Type 4: Type 1 and then Type 2 ("Combined")
  #Type 5: Type 1 and then Type 3 ("Hourly Combined")
  {
    for(homogType in 2)
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
    {
      for(order in c("hFirst", "oFirst"))
      {
        for( qcIterations in 2:3 )
        {
          #for(estimator in c("tukey", "huber"))
          for(estimator in c("huber"))
          {
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
              if( (order=="hFirst" & testNo%%2==1) | (order=="oFirst" & testNo%%2==0) )
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
                if(exists("breaks"))
                  breaks = rbind( breaks, data.frame(out[[2]]) )
                else
                  breaks = data.frame(out[[2]])
              }
  
              #Outlier Test
              if( (order=="oFirst" & testNo%%2==1) | (order=="hFirst" & testNo%%2==0) )
              {
                #(Re)compute the outlierFl.  Set to F to clear any previous detections
                sim.dat$outlierFl = F
                #types 1, 4, and 5 all have "Tier 1"
                if(outlierType %in% c(1, 4, 5))
                {
                  if(estimator=="tukey")
                    score = tukey(sim.dat, dayWindow=365, hourBuckets=1)
                  else
                    score = huber(sim.dat, dayWindow=365, hourBuckets=1)
                  sim.dat$outlierFl[abs(score) > 6] = T
                }
                #type 2 and 4 have "Tier 2"
                if(outlierType %in% c(2, 4))
                {
                  if(estimator=="tukey")
                    score = tukey(sim.dat, dayWindow=45, hourBuckets=1)
                  else
                    score = huber(sim.dat, dayWindow=45, hourBuckets=1)
                  sim.dat$outlierFl[abs(score) > 5] = T
                }
                #type 3 and 5 have "Hourly Tier 2"
                if(outlierType %in% c(3, 5))
                {
                  if(estimator=="tukey")
                    score = tukey(sim.dat, dayWindow=45, hourBuckets=2)
                  else
                    score = huber(sim.dat, dayWindow=45, hourBuckets=2)
                  sim.dat$outlierFl[abs(score) > 5] = T
                }
              }
            } #Close testNo for loop
  
            #Compute some statistics to save:
            breakDates = merge( sim.dat, breaks, by=c("Year", "Day_Of_Year", "Hour") )$Date.x[1:5]
            rawRMSE = sqrt( mean( (contam$Simulated - orig$Simulated)[sim.dat$outlierSize==0]^2 ) )
            homoRMSE = sqrt( mean( (sim.dat$Simulated - orig$Simulated)[sim.dat$outlierSize==0]^2 ) )
            #RMSE with outliers included in calc.  Problematic since orig is before outliers were added.
            #rawRMSE = sqrt( mean( (contam$Simulated - orig$Simulated)[sim.dat$outlierSize==0]^2 ) )
            #homoRMSE = sqrt( mean( (sim.dat$Simulated - orig$Simulated)^2 ) )
            
            toBind = data.frame(
              #seed as ID
               ID=seed
              #outlier type, converted to human-readable form
              ,outlierType=ifelse(outlierType==1, "Tier 1", ifelse(outlierType==2, "Tier 2", ifelse(outlierType==3, "Tier 2 Hourly", ifelse(outlierType==4, "Combined", "Combined Hourly") ) ) )
              #homogeneity type
              ,homogType=homogType
              #Which test was done first?
              ,order=ifelse(order=="oFirst", "Outlier First", "Homogenization First")
              #Huber vs. Tukey
              ,estimator=estimator
              #Window size for days
              ,windowSize=ifelse(outlierType==1,365,45)
              #False positive count
              ,FP_Cnt=sum(sim.dat$outlierFl & sim.dat$outlierSize==0)
              #True positive count
              ,TP_Cnt=sum(sim.dat$outlierFl & sim.dat$outlierSize!=0)
              #False negative count
              ,FN_Cnt=sum(!sim.dat$outlierFl & sim.dat$outlierSize!=0)
              #True negative count
              ,TN_Cnt=sum(!sim.dat$outlierFl & sim.dat$outlierSize==0)            
              #True positive rates for specific levels:
              ,TPR_4=sum(sim.dat$outlierFl & sim.dat$outlierSize==4)/sum(sim.dat$outlierSize==4)
              ,TPR_5=sum(sim.dat$outlierFl & sim.dat$outlierSize==5)/sum(sim.dat$outlierSize==5)
              ,TPR_6=sum(sim.dat$outlierFl & sim.dat$outlierSize==6)/sum(sim.dat$outlierSize==6)
              ,detBr1Time =breakDates[1]
              ,detBr2Time =breakDates[2]
              ,detBr3Time =breakDates[3]
              ,detBr4Time =breakDates[4]
              ,detBr5Time =breakDates[5]
              ,detBr1Size =ifelse(is.null(breaks[1,5][[1]]),NA,breaks[1,5][[1]])
              ,detBr2Size =ifelse(is.null(breaks[2,5][[1]]),NA,breaks[2,5][[1]])
              ,detBr3Size =ifelse(is.null(breaks[3,5][[1]]),NA,breaks[3,5][[1]])
              ,detBr4Size =ifelse(is.null(breaks[4,5][[1]]),NA,breaks[4,5][[1]])
              ,detBr5Size =ifelse(is.null(breaks[5,5][[1]]),NA,breaks[5,5][[1]])
              ,efficiencyRMSE=(rawRMSE - homoRMSE)/rawRMSE
              ,homoRMSE=homoRMSE
              ,rawRMSE=rawRMSE
              ,qcIterations
            )
            if(exists("detectionStat"))
              detectionStat = rbind(detectionStat,toBind)
            else
              detectionStat = toBind

          } #Close estimator for loop
        } #Close qcIterations loop
      } #Close order for loop
    } #Close homogType for loop
  } #Close outlierType for loop


#   if(i %% 10==0)
#     cat("Run",i,"completed.  nrow(dataStat):",nrow(dataStat),"  nrow(det...):", nrow(detectionStat), "\n")
#   if(i %% 50==0)
#     save(detectionStat, dataStat, file=paste0("Simulations_nonRobust_",nrow(dataStat),"_",Sys.info()[4],"_",runId,".RData"))
  cat(".")
  return(merge(dataStat, detectionStat))
}

cat("Saving results...\n")
results2 = do.call("rbind", results)
save(results, results2, file=paste0("Simulations_nonRobust_",length(results),"_",Sys.info()[4],"_",runId,".RData"))
cat("Elapsed time:", difftime(Sys.time(), start, units = "hours"), "hours")

# Do a new simulation by re-sourcing itself
source("run_simulation_server.R")