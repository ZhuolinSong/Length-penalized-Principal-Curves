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


cv_prin <- princurve_simulation(ncores, s_n, s_k, s_q,
                maxit = 1e2, thresh = 1e-5,
                s_range = 10)


save(cv_prin, file = "cv_prin.RData")

