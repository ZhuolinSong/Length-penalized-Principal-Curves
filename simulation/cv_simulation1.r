library(devtools)
devtools::load_all()

library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(999983)


# para setup
s_q <- 2L
ncores <- 16
s_n <- 128L
s_k <- 32L


cv1 <- cv_simulation(ncores, s_n, s_k, s_q,
                maxit = 1e4, thresh = 1e-4,
                s_range = 5, l_range = 10, lambda_l = 5)


save(cv1, file = "cv1.RData")
