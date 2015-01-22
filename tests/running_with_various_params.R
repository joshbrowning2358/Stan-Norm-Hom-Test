# install.packages("~/GitHub/Stan-Norm-Hom-Test/snht_1.0.1.tar.gz", type = "src",
#                  repo = NULL)
library(snht)

## Functions to test:
##      snht
##      robustSNHT
##      robustSNHTunequal
##      pairwiseSNHT

i = 0
while(TRUE){
    n = round(exp(runif(1,min=0, max=7)))
    data = rnorm(n)
    period = runif(1, min = 1, max = n/2)
    robust = sample(c(T, F), size = 1)
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
    snht(data, period, robust, time, scaled, rmSeasonalPeriod)
    i = i + 1
    cat("Run number", i, "completed.\n")
}