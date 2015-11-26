library(mvtnorm)
library(ggplot2)
n = 10000
rho = 0.0

d = rmvnorm(n = n, sigma = matrix(c(1, rho, rho, 1), nrow = 2))
cor(d)
d = data.frame(d)
fit = lm(I(X1 - X2) ~ X2 + 0, data = d)
fit

ggplot(d, aes(x = X2, y = X1 - X2)) + geom_point(alpha = 1)
