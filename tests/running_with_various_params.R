install.packages("~/GitHub/Stan-Norm-Hom-Test/snht_1.0.1.tar.gz", type = "src",
                 repo = NULL)
library(snht)

## Functions to test:
##      snht
##      robustSNHT
##      robustSNHTunequal
##      pairwiseSNHT

i = 0
while(TRUE){
    n = round(exp(runif(1,min=1.61, max=10)))
    data = rnorm(n)
    period = runif(1, min = 1, max = n/2)
    while(2 * round(period) >= n - 1)
        period = period - 1 # don't let period get too big
    robust = sample(c(T, F), size = 1)
    if(n <= 10) # Huber estimator assumes it has >= 5 obs
        robust = FALSE
    if(sample(c(T, F), size = 1)){
        time = 1:n + rnorm(n, sd = .3)
    } else {
        time = NULL
    }
    scaled = sample(c(T, F), size = 1)
    if(sample(c(T, F), size = 1)){
        rmSeasonalPeriod = Inf
    } else {
        rmSeasonalPeriod = runif(1, min = 1, max = n/2)
    }
    stat = snht(data, period, robust, time, scaled, rmSeasonalPeriod)
    i = i + 1
    cat("Run number", i, "completed.\n")
}