library(devtools)
devtools::load_all()
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(999983)

# para setup
s_q <- 2L
ncores <- 1
s_n <- 128L
s_k <- 32L


cv_test <- grid_simulation(ncores, s_n, s_k, s_q,
                maxit = 5, thresh = 1e-3,
                s_range = 10, l_range = 8)


save(cv_test, file = "cv_test.RData")

analysis_plot(cv_test, 1:7)
