library(data.table)
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
if(Sys.info()[1]=="Linux" & Sys.info()[4]=="jb"){
  if(getwd() != "/media/storage/Professional Files/Mines/Research/Wind QC"){
    setwd("/media/storage/Professional Files/Mines/Research/Wind QC/")
    source("/media/storage/Github/Stan-Norm-Hom-Test/simulation_functions.R")
    library(snht)
    removeSeasonalPeriod = snht:::removeSeasonalPeriod
  }
  suppressPackageStartupMessages(library(doParallel))
  library(foreach)
  registerDoParallel(4)
}
if(Sys.info()[1]=="Linux" & Sys.info()[4]=="jmds"){
  # Only run this if code hasn't been run already
  if(getwd() != "~/GitHub/Stan-Norm-Hom-Test/"){
    setwd("~/GitHub/Stan-Norm-Hom-Test/")
    source("simulation_functions.R")
    library(snht)
    removeSeasonalPeriod = snht:::removeSeasonalPeriod
  }
  suppressPackageStartupMessages(library(doParallel))
  library(foreach)
  registerDoParallel(cores=detectCores(all.tests=TRUE))
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

# load("Results/Simulation_20151222 (robust)/All_Results.RData")
# robustResults = finalResults
# load("Results/Simulation_20151222 (nonrobust)/All_Results.RData")
# nonrobustResults = finalResults
# robustResults$robust = TRUE
# nonrobustResults$robust = FALSE
# detection = rbind(robustResults, nonrobustResults)
# rm(nonrobustResults, robustResults, finalResults, params)
# detection$simBrCnt = NA
load("Results/Simulation_20151222 (robust)/comparison_dataframe.RData")
detection = data.table(detection)

uniqueLevels = detection[, .N, c("station", "pressure", "n")][, N := NULL]

simBrkCnt = foreach(i = 1:nrow(uniqueLevels)) %dopar% {
    station = uniqueLevels[i, station]
    dataset = "trh" #allowable values are "trh", "trhc", "wind", "windc"
    pressure = uniqueLevels[i, pressure]
    startYear = 2000
    endYear = startYear + as.numeric(uniqueLevels[i, n]/365) - 1

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

    seeds = detection[uniqueLevels[i, ], unique(ID),
                      on = c("station", "pressure", "n")]

    results = lapply(seeds, function(seed){
        set.seed(seed)
        
        #Simulate new dataset
        sim.dat = simBaseTS(fitTS, data=test.data, startYear=startYear, endYear=endYear)
        sim.dat = simResiduals2(fitST, sim.dat, ar1=ar0.5^2)
        orig = sim.dat
        
        #Contaminate with outliers and breaks
        sim.dat = contaminateOutlier(sim.dat, "extreme", sample(c(0, 1,2,5,10)/100,size=1), muError=c(8,9,10))
        sim.dat = contaminateBreak(sim.dat, breaks = sample(1:3*(endYear - startYear + 1)/10,size=1),
                                   sigmaMult=c(.2,.4,.6))
        data.frame(seed, simBrCnt = length(unique(sim.dat$Breaks)) - 1)
    })
    out = do.call("rbind", results)
    out = merge(out, uniqueLevels[i, ])
    save(out, file = paste(station, pressure, (endYear-startYear+1), ".RData", sep = "_"))
    out
}
save(simBrkCnt, file = "simBrkCnt.RData")