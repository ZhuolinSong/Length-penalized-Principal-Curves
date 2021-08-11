library(princurve)

devtools::load_all() # our package
# para setup
s_n <- 100L
s_K <- 25L
s_q <- 2L


# Example 1. bivariate spherical gaussial data
m_y <- matrix(rnorm(s_n * s_q, 0, 1), s_n, s_q, byrow = FALSE)

result <- f_ppslp(m_y, s_K, plot_update = T)
fit <- principal_curve(m_y, cv = F)
lines(fit, type = "l", col = "green", lwd = 2, lty = 1)

# Example 2.Circle data withou noise
m_y <- y_generator(1)

result <- f_ppslp(m_y, s_K, plot_update = T, start_end = xy_fix(1))
fit <- principal_curve(m_y)
lines(fit, type = "l", col = "green", lwd = 2, lty = 1)

spar <- cv_smooth(y_generator(3), s_K, maxit = 1e2, thresh = 1e-3,
                s_range = 10, l_range = c(0, exp(-5), exp(-4)))
len <- cv_length(y_generator(3), s_K, maxit = 1e2, thresh = 1e-3,
                s_range = unique(spar), l_range = c(0, exp(-5), exp(-4)))

len <- cv_length(y_generator(3), s_K, maxit = 1e2, thresh = 1e-3,
                s_range = 0.1 * 3:10, l_range = 10)
spar <- cv_smooth(y_generator(3), s_K, maxit = 1e2, thresh = 1e-3,
                s_range = 0.1 * 4:10, l_range = unique(len))
result <- f_ppslp(y_generator(3), s_K,
                    unique(spar), unique(len), plot_update = T)

# Example 3.Circle data with noise
init <- xy_truth(3, s_K)
m_y <- y_generator(2)

# 3.A init as principal component
result <- f_ppslp(m_y, s_K, 0.5, 0.01, plot_update = T)
fit <- principal_curve(m_y)
lines(fit, type = "l", col = "green", lwd = 2, lty = 1)
lines(init, col = "black", lty = 3)

# 3.B init as the circle
result <- f_ppslp(m_y, s_K, plot_update = T, init = init, start_end = xy_fix(3))
fit <- principal_curve(m_y, start = init, thresh = 1e-10, cv = F)
lines(fit, type = "l", col = "green", lwd = 2, lty = 1)
lines(init, col = "black", lty = 3)



# Example 4.Hastie bias examples
m_y <- y_generator(4)

result <- f_ppslp(m_y, s_K, plot_update = T, start_end = xy_fix(4))
fit <- principal_curve(m_y)
lines(fit, type = "l", col = "green", lwd = 2, lty = 1)
lines(xy_truth(4), col = "black", lty = 3)



# Example 5.HuanCHen Three MNIST dataset
m_y <- y_generator(7)
result <- f_ppslp(m_y, s_k = 64, start_end = xy_fix(7))

library(ppclp) # HuanChen's method
library(tidyverse) # required by ppclp
tmpCurve <- ppclp2D(threeExample$x, threeExample$y, threeExample$xFix, threeExample$yFix)

# Hastie's method
fit <- principal_curve(m_y, start = cbind(tmpCurve$xFit, tmpCurve$yFit))

# plotting for comparison
pp_plot(result)
points(threeExample$xFix, threeExample$yFix, pch = 16, cex = 1.5, col = "red")
lines(tmpCurve$xFit, tmpCurve$yFit, type = "l", col = "red", lwd = 2, lty = 1)
lines(fit, type = "l", col = "green", lwd = 2, lty = 1)

legend("bottomleft",
    c("Our method", "HuanChen", "Hastie"),
    lwd = c(1, 2, 1),
    col = c("blue", "red", "green"),
    lty = 1)


# CV Tunning for Princurve
m_y <- y_generator(7)
result <- f_ppslp(m_y, s_K, plot_update = T)
fit <- cv_princurve(m_y, s_K, maxit = 1e3, thresh = 1e-3)
lines(fit, type = "l", col = "green", lwd = 2, lty = 1)