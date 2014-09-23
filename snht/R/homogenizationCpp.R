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
    leftMean = leftSum/(center-leftBound);
    rightMean = rightSum/(rightBound-center);
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