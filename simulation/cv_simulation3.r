library(devtools)
devtools::load_all()

library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(999983)


# para setup
s_q <- 2L
ncores <- 16
s_n <- 128L
s_k <- 64L


cv3 <- cv_simulation(ncores, s_n, s_k, s_q,
                maxit = 1e4, thresh = 1e-3,
                s_range = 4, l_range = 8)

save(cv3, file = "cv3.RData")
