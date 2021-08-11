#' Grid Cross Validation simulation
#'
#' parallel simulation
#'
#' @param ncores number of cores used in simulation
#' @param s_n number of observations
#' @param s_k number of knots
#' @param s_q dimension
#' @param ... other args for cv_grid
#'              (maxit, thresh, l_range, s_range, folds, lambda_l)
#' @import parallel
#' @return a list consisting fitted ppslp object for each simulation
#' @keywords Principal curve
#' @export
grid_simulation <- function(ncores = 16, s_n = 100L,
                    s_k = 25L, s_q = 2L, ...) {

  grid <- mclapply(1:7, inner_loop <- function(case) {
        init <- NULL

        if (case == 3 || case == 6) {
          init <- xy_truth(case, s_k)
        }

        m_y <- y_generator(case, s_n, s_q)
        output <- ppslp::cv_grid(m_y, s_k, init = init, ...)
        return(output)

  }, mc.cores = ncores)

  return(grid)
}

#' Princurve Cross Validation simulation
#'
#' parallel simulation
#'
#' @param ncores number of cores used in simulation
#' @param s_n number of observations
#' @param s_k number of knots
#' @param s_q dimension
#' @param ... other args for cv_grid
#'              (maxit, thresh, s_range, folds, lambda_l)
#' @import parallel
#' @return a list consisting fitted ppslp object for each simulation
#' @keywords Principal curve
#' @export
princurve_simulation <- function(ncores = 16, s_n = 100L,
                    s_k = 25L, s_q = 2L, ...) {

  fit <- mclapply(1:7, inner_loop <- function(case) {
        init <- NULL

        if (case == 3 || case == 6) {
          init <- xy_truth(case, s_k)
        }

        m_y <- y_generator(case, s_n, s_q)
        output <- ppslp::cv_princurve(m_y, s_k, start = init, ...)
        return(output)

  }, mc.cores = ncores)

  return(fit)
}


#' Smoothing parameter Cross Validation simulation
#'
#' parallel simulation
#'
#' @param ncores number of cores used in simulation
#' @param s_n number of observations
#' @param s_k number of knots
#' @param s_q dimension
#' @param ... other args for cv_smooth
#'              (maxit, thresh, l_range, s_range, folds, lambda_l)
#' @import parallel
#' @return a list consisting fitted ppslp object for each simulation
#' @keywords Principal curve
#' @export
smooth_simulation <- function(ncores = 16, s_n = 100L,
                    s_k = 25L, s_q = 2L, ...) {

  out <- mclapply(1:7, inner_loop <- function(case) {
        init <- NULL

        if (case == 3 || case == 6) {
          init <- xy_truth(case, s_k)
        }

        m_y <- y_generator(case, s_n, s_q)
        output <- ppslp::cv_smooth(m_y, s_k, init = init, ...)
        return(output)

  }, mc.cores = ncores)

  return(out)
}



#' Smoothing parameter Cross Validation simulation
#'
#' parallel simulation
#'
#' @param ncores number of cores used in simulation
#' @param s_n number of observations
#' @param s_k number of knots
#' @param s_q dimension
#' @param ... other args for cv_smooth
#'              (maxit, thresh, l_range, s_range, folds, lambda_l)
#' @import parallel
#' @return a list consisting fitted ppslp object for each simulation
#' @keywords Principal curve
#' @export
length_simulation <- function(ncores = 16, s_n = 100L,
                    s_k = 25L, s_q = 2L, ...) {

  out <- mclapply(1:7, inner_loop <- function(case) {
        init <- NULL

        if (case == 3 || case == 6) {
          init <- xy_truth(case, s_k)
        }

        m_y <- y_generator(case, s_n, s_q)
        output <- ppslp::cv_length(m_y, s_k, init = init, ...)
                  return(output)

  }, mc.cores = ncores)

  return(out)
}


#' Simulation Examples
#'
#' @param case number indicates the examples(1-7):
#'              1: 3/4 circle
#'              2: 3/4 circle with noise
#'              3: 3/4 circle with noise, initialize with the truth
#'              4: 1/2 circle
#'              5: 1/2 circle with noise
#'              6: 1/2 circle with noise, initialize with the truth
#'              7: the shape of number 3 from mnist dataset
#' @param n number of obervation in the simulation
#' @param q dimension of the simulation
#' @return the simulation examples
#' @keywords Principal curve
#' @export
y_generator <- function(case, n=100L, q=2L, seeds=999983) {
  if (case == 1) {
    noise <- runif(n, 0, 1.5 * pi)
    return(matrix(5 * c(sin(noise), cos(noise)), n, q, byrow = FALSE))
  } else if (case == 2 || case == 3) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seeds)
    noise <- runif(n, 0, 1.5 * pi)
    truth <- matrix(5 * c(sin(noise), cos(noise)), n, q, byrow = FALSE)
    return(truth + matrix(rnorm(q * n), n, q))
  } else if (case == 4) {
    mean <- runif(n, 0, pi)
    return(matrix(5 * c(cos(mean), sin(mean)), n, q, byrow = FALSE))
  } else if (case == 5 || case == 6) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seeds + 1)
    mean <- runif(n, 0, pi)
    mean <- 5 * c(cos(mean), sin(mean))
    return(matrix(rnorm(n = q * n, mean), n, q, byrow = FALSE))
  } else if (case == 7) {
    return(cbind(threeExample$x, threeExample$y))
  }
}


#' Start and end point of the simulation Examples
#'
#' @param case number indicates the examples(1-7):
#'              1: 3/4 circle
#'              2: 3/4 circle with noise
#'              3: 3/4 circle with noise, initialize with the truth
#'              4: 1/2 circle
#'              5: 1/2 circle with noise
#'              6: 1/2 circle with noise, initialize with the truth
#'              7: the shape of number 3 from mnist dataset
#' @return the start_end point for each case
#' @keywords Principal curve
#' @export
#' @return start_end points (matrix):
#'         [start, end]^T
xy_fix <- function(case) {
  if (case == 1 || case == 2 || case == 3) {
    return(rbind(c(0, 5), c(-5, 0)))
  } else if (case == 4 || case == 5 || case == 6) {
    return(rbind(c(5, 0), c(-5, 0)))
  } else if (case == 7) {
    return(cbind(threeExample$xFix, threeExample$yFix))
  }
}

#' The truth curve of the simulation Examples
#'
#' @param case number indicates the examples(1-7):
#'              1: 3/4 circle
#'              2: 3/4 circle with noise
#'              3: 3/4 circle with noise, initialize with the truth
#'              4: 1/2 circle
#'              5: 1/2 circle with noise
#'              6: 1/2 circle with noise, initialize with the truth
#'              7: the shape of number 3 from mnist dataset
#' @param n number of the truth
#' @param q dimension of the simulation
#' @return the start_end point for each case
#' @keywords Principal curve
#' @export
xy_truth <- function(case, n=10000, q=2) {
  if (case == 1 || case == 2 || case == 3) {
    noise <- 1.5 * pi / n * seq_len(n)
    return(matrix(5 * c(sin(noise), cos(noise)), n, q, byrow = FALSE))
  } else if (case == 4 || case == 5 || case == 6) {
    noise <- pi / n * seq_len(n)
    return(matrix(5 * c(cos(noise), sin(noise)), n, q, byrow = FALSE))
  } else if (case == 7) {
    return(matrix(0, n, q, byrow = FALSE))
  }
}


mse_calculate <- function(case, theta, n=1e5) {
    truth <- xy_truth(case, n)
    proj <- ppslp::project(truth, theta, diff(theta))
    return(proj$dist / n)
}


#' Analysis plot
#'
#' plot simulation results
#'
#' @param cv list of object ppslp of the simulation
#' @param cases vector of the example simulations
#' @param benchmark indicate whether to plot hastie
#' @return mse
#' @import princurve
#' @keywords Principal curve
#' @export
analysis_plot <- function(cv, cases, benchmark = T) {
    mse <- sapply(cases, outer_loop <- function(case) {
                output <- cv[[case]]
                if (is.atomic(output)) { # check bugs !
                    print(case)
                }

                pp_plot(output)

                theta <- output$theta
                dist <- mse_calculate(case, theta, 1e5)

                if (benchmark) {
                    y <- output$y
                    fit <- principal_curve(y)
                    lines(fit, type = "l", col = "green", lwd = 2, lty = 1)
                    legend("topleft",
                        c("Our method", "Hastie"),
                        lwd = 2,
                        col = c("blue",  "green"),
                        lty = 1)

                    theta <- fit$s
                    dist[2] <- mse_calculate(case, theta)
                }

                lines(xy_truth(case), lty = 3)
                return(dist)
    })

    rownames(mse) <- c("ppslp", "Hastie's")
    return(mse)
}