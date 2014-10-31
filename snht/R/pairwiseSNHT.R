pairwiseSNHT <- function(data, dist, k, period, crit=100, returnStat=FALSE, ...){
  #data quality checks
  stopifnot(is(data,"data.frame"))
  if(ncol(data)==2)
    stopifnot(colnames(data)==c("data","location"))
  if(ncol(data)==3)
    stopifnot(colnames(data)==c("data","location","time"))
  stopifnot(ncol(data) %in% c(2,3))
  locs = as.character(unique(data$location))
  stopifnot(all(rownames(dist) %in% locs))
  stopifnot(all(colnames(dist) %in% locs))
  stopifnot(all(locs %in% rownames(dist)))
  stopifnot(all(locs %in% colnames(dist)))
  stopifnot(k>=1) #Must have at least one neighbor
  stopifnot(k<=length(locs)-1) #Can have at most length(locs)-1 neighbor, since self can't be used  
  stopifnot(diag(dist)==0)
  if(any(dist[row(dist)!=col(dist)]<=0))
    stop("Off diagonal elements of dist must be >0")
  
  #A column of the pairs object defines the neighbors for that location
  pairs = lapply(1:nrow(dist), function(i){
    x = dist[i,]
    filt = rank(x)>1 & rank(x)<=k+1.5
    colnames(dist)[filt]
  })
  names(pairs) = colnames(dist)
  #Change pairsDf to uniquePairs, postprocess pairs so it's always the same object type
  uniquePairs = data.frame(loc1 = rep(colnames(dist), lapply(pairs,length))
                          ,loc2 = do.call("c", pairs), stringsAsFactors=F )
  #Remove non-unique pairs.  Start at nrow(pairsDf) and work down.  Otherwise, the
  #iterator would break when you delete a row, as the nrow(pairsDf) would reduce but
  #the iterator would still run to the original nrow(pairsDf)
  for(i in nrow(uniquePairs):1)
    if(any(uniquePairs[,1]==uniquePairs[i,2] & uniquePairs[,2]==uniquePairs[i,1]))
      uniquePairs = uniquePairs[-i,]
  
  #Add times if they don't already exist (just 1:nrow()).
  if(!"time" %in% colnames(data)){
    if(length( unique( table( data$location ) ) ) != 1){
      stop("All locations must have the same number of obs if time is not provided!
           May need to remove unused levels in data.")
    }
    data$order = 1:nrow(data) #ensure original ordering is preserved
    data = ddply(data, "location", function(df){
      df = df[order(df$order),]
      df$time = 1:nrow(df)
      return(df)
    } )
    data$order = NULL
  }
  
  #Restructure data
  data = cast(data, formula = time ~ location, value="data")
  diffs = data.frame(time=data$time)
  for(i in 1:nrow(uniquePairs)){
    diffs = cbind(diffs, data[,uniquePairs[i,1]] - data[,uniquePairs[i,2]])
    colnames(diffs)[ncol(diffs)] = paste0(uniquePairs[i,1],"-",uniquePairs[i,2])
  }
  
  #Compute snht statistics
  statistics = apply(diffs[,-1], 2, snht, period=period, time=diffs[,1], ...)
  avgDiff = do.call("cbind", lapply(statistics, function(x) x$rightMean-x$leftMean ) )
  statistics = do.call("cbind", lapply(statistics, function(x) x$score))
  if(returnStat)
    return(statistics)

  #Create candidate, a matrix where the (i,j)th entry corresponds to the number
  #of changepoints in difference series for location j that occured at time i.
  #For example, suppose location i_1 was paired with i_2, i_3, and i_4.  If the
  #statistic for i_1-i_2 and i_1-i_4 exceeded the threshold at time j, then
  #candidate_{i,j} = 2.
  candidate = matrix(0, nrow=nrow(data), ncol=length(locs))
  colnames(candidate) = locs
  for( j in 1:ncol(statistics) ){
    name = colnames(statistics)[j]
    name = strsplit(name, "-")[[1]]
    delta = as.numeric(statistics[,j]>crit)
    delta[is.na(delta)] = 0
    if(name[2] %in% pairs[name[1]][[1]])
      candidate[,name[2]] = candidate[,name[2]] + delta
    if(name[1] %in% pairs[name[2]][[1]])
      candidate[,name[1]] = candidate[,name[1]] + delta
  }
  
  #"Unconfound" the candidate matrix by assigning a changepoint to the location of
  #largest count.
  #Algorithm:
  #maxCols = locations which all attain the max count
  #If (|maxCols|>1)
  #  Examine all difference series for all elements of maxCols but restricted to rows
  #  where the max count occurs.  Find the statistic with the largest value.  Call the
  #  two locations forming this series col_A, col_B.  Set the break time to the
  #  time of this statistic, call it brkT.
  #    brkCol = col_B or col_A, respectively
  #  else
  #    At time brkT, examine the other difference statistics for col_A and col_B.  Assign
  #    a break to whichever pair has the second largest statistic at time brkT
  #else
  #  brkCol = column with largest count
  #  Amongst all points with max count, find the time where the largest statistic
  #  occurs.  Call this time brkT
  #
  #Update candidate: Set t_window=brkT + -period:period.
  #candidate[t_window, brkCol] = 0
  #brkCol_pairs = set of all time series that use maxCol in their difference series.
  #candidate[t_window, brkCol_pairs] = pmax( candidate[t_window, brkCol_pairs]-1, 0)
  #Note: statistics must also be updated, to ensure that times from other stations are
  #not choosen when the difference is due to the homogenized station.
  #The assumption is that the changepoint in the maxCol series was the issue in the 
  #other difference series.  This seems pretty reasonable.
  
  breaks = matrix(nrow=0, ncol=3)
  while(any(candidate>0)){
    colMax = apply(candidate, 2, max)
    maxVal = max(colMax)
    maxCols = colMax==maxVal
    
    #Determine which difference series we'll need to examine
    pairsMax = pairs[maxCols]
    pairsMax = data.frame(loc1 = rep(names(pairsMax), lapply(pairsMax,length))
                         ,loc2 = do.call("c", pairsMax) )
    #Could occur in either order, create both possibilities
    pairsMax$diff = paste0(pairsMax[,1], "-", pairsMax[,2])
    pairsMax$diff = ifelse(!pairsMax$diff %in% colnames(statistics)
                          ,paste0(pairsMax[,2], "-", pairsMax[,1])
                          ,pairsMax$diff )
    
    #Determine brkT, the time at which the current break is detected, and brkCol.
    #Note: if columns A and B have the same value in candidate and the maximum is
    #between A and B, then the first column will be used.  Tie-breaking in this case
    #is non-trivial (what if each only have candidate=1?  how can we pull other values?)
    brkTimes = lapply(names(maxCols), function(location){
      if(!maxCols[location])
        return(NULL)
      rows = which(candidate[,location]==maxVal)
      cols = pairsMax$diff[pairsMax$loc1==location]
      currCol = which.max(apply(statistics[rows,cols,drop=F], 2, max))
      currColTime = rows[which.max(statistics[rows,currCol])]
      return(data.frame(time=currColTime, stat=max(statistics[currColTime,]) ) )
    } )
    brkTimes = do.call("rbind", brkTimes)
    brkT = brkTimes$time[which.max(brkTimes$stat)]
    brkCol = names(maxCols)[maxCols][which.max(brkTimes$stat)]
    
    #Update candidate matrix
    tWindow = brkT + -period:period #No changepoints within period observations
    candidate[tWindow, brkCol] = 0
    adjCandCols = lapply(pairs, function(x){brkCol %in% x})
    adjCandCols = names(pairs)[do.call("c",adjCandCols)]
    candidate[tWindow, adjCandCols] = pmax( candidate[tWindow, adjCandCols]-1, 0)
    
    #Update statistics matrix
    adjStatCols = c(paste0(brkCol,"-",adjCandCols), paste0(adjCandCols,"-",brkCol))
    adjStatCols = adjStatCols[adjStatCols %in% colnames(statistics)]
    #Set to zero as we're assuming large values are caused by brkCol
    statistics[tWindow, adjStatCols] = 0
    
    #Update data
    adjMeanCols = c(paste0(brkCol,"-",pairs[[brkCol]]), paste0(pairs[[brkCol]],"-",brkCol))
    adjMeanCols = adjMeanCols[adjMeanCols %in% colnames(statistics)]
    shift = mean(avgDiff[brkT,adjMeanCols], na.rm=T)
    data[brkT:nrow(data),brkCol] = data[brkT:nrow(data),brkCol] - shift
    
    #Append detected break to breaks
    breaks = rbind(breaks, c(brkT, brkCol, shift) )
  }
  breaks = data.frame(breaks)
  colnames(breaks) = c("time", "location", "shift")
  breaks$time = as.numeric(breaks$time)
  breaks$shift = as.numeric(breaks$shift)
  
  data = melt(data, id.vars=c("time"))
  rownames(data) = NULL
  colnames(data)[colnames(data)=="value"] = "data"
  data = data[,c("data", "location", "time")]
  return(list(data=data, breaks=breaks))
}
