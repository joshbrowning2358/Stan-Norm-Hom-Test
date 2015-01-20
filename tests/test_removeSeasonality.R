context("Remove Seasonality")

test_that("Pure sin is reduced to 0", {
    x = sin(1:100*pi/10)
    out = removeSeasonalPeriod(x, period = 20)
    expect_true( max( abs(out) ) <= .01 )
})