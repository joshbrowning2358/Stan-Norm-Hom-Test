context("Basic SNHT Functionality")

test_that("Functions return appropriate values", {
    stat = snht( rnorm(500), period = 100)
    expect_that( stat, is_a("data.frame") )
    expect_that( colnames(stat), equals(c("score", "leftMean", "rightMean") ) )
    expect_that( dim(stat), equals(c(500, 3) ) )
    stat = robustSNHT( rnorm(500), period = 100)
    expect_that( stat, is_a("data.frame") )
    expect_that( colnames(stat), equals(c("score", "leftMean", "rightMean") ) )
    expect_that( dim(stat), equals(c(500, 3) ) )
    stat = robustSNHTunequal( rnorm(500), period = 100, time = 1:500)
    expect_that( stat, is_a("data.frame") )
    expect_that( colnames(stat),
        equals(c("score", "leftMean", "rightMean", "time") ) )
    expect_that( dim(stat), equals(c(500, 4) ) )
})

library(mgcv)
test_that("Arguments work", {
    x = rnorm(100)
    stat = snht(data = x, period = 10, robust = T, time = 1:100 + rnorm(100),
        scaled = FALSE, rmSeasonalPeriod = 30 )
    expect_that( colnames(stat),
        equals(c("score", "leftMean", "rightMean", "time") ) )
    stat = snht(data = x, period = 30, robust = T, scaled = FALSE,
        rmSeasonalPeriod = 10 )
    expect_that( colnames(stat),
        equals(c("score", "leftMean", "rightMean") ) )
})

test_that("Scale gives same answer on scaled data", {
    x = rnorm(100)
    stat1 = snht( data = x, period = 10, scaled = TRUE )
    stat2 = snht( data = 100*x, period = 10, scaled = TRUE )
    expect_that( stat1[,1], equals(stat2[,1]) )
    expect_that( stat1[,2:3], equals(stat2[,2:3]/100) )
})

