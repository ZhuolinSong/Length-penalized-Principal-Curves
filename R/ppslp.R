#' Principal Curve with Length Penalty
#'
#' This function allows you to fit a principal curve with length penalty.
#' @param m_y your data matrix y
#' @param s_k number of knots
#' @param maxit maximum number of iterations, default to be 60
#' @param spar smoother parameter, usually in (0, 1]
#' @param scaled_lam scaled length parameter,
#'                     $=\lambda_l/\lambda_{\max}$, default to be 0.01
#' @param init either a previously fit principal curve,
#'               or starting value of matrix theta.
#'               If missing or NULL, then the first principal component is used.
#'              If  then a circle is used as the start.
#' @param thresh convergence treshhold
#' @param plot_update plot the update
#' @param start_end start and end constraint points(matrix):
#'                  [start, end]^T, default to NULL no constraints
#' @param lambda_s smoother parameter, default to be NULL
#' @param lambda_l length parameter, default to be NULL
#' @import splines
#' @return \code{out} An prinicipal curve opject includes
#'        \item{y}{matrix y}
#'        \item{fitted_y}{fitted y on principal curve}
#'        \item{B}{matrix B}
#'        \item{theta}{spline matrix theta}
#'        \item{iter}{number of iteration}
#'        \item{lambda_l}{length parameter}
#'        \item{lambda_s}{smoother parameter}
#'        \item{dist}{difference of theta}
#'        \item{converged}{A logical indicating whether the algorithm converged
#'              or not.}
#'        \item{num_iterations}{Number of iterations
#'                              completed before returning.}
#'        \item{call}{the call that created this object; allows it to be
#'              \code{updated()}.}
#' @keywords Principal curve
#' @export
#' @examples
#' \dontrun{
#' devtools::load_all() # our package
#' data("threeExample")
#' y <- cbind(threeExample$x, threeExample$y)
#' fit <- f_ppslp(y, k = 50, lambda_s = 10, lambda_l = 0.1)
#' plot(fit)
#' lines(fit)
#' points(fit)
#' whiskers(y, fit$fitted_y)
#' points(y, pch = 16, cex = 0.8, col = "red")
#' }
f_ppslp <- function(m_y,
                    s_k,
                    spar = 0.6,
                    scaled_lam = 0.01,
                    maxit = 1e3,
                    init = NULL,
                    thresh = 1e-3,
                    plot_update = F,
                    lambda_l = NULL,
                    lambda_s = NULL,
                    start_end=NULL) {

  #-------------------------------------------------------
  ##### check up #####
  #-------------------------------------------------------
  if (is.null(lambda_l)) {
    if (scaled_lam >= 1 || scaled_lam < 0) {
      stop("Invalid input for scaled length parameter")
    }
    if (scaled_lam >= 0.5) {
        print("Warning: high scaled length parameter!")
    }
  }

  #-------------------------------------------------------
  ##### set up #####
  #-------------------------------------------------------
  function_call <- match.call()
  dist_old <- sum(diag(stats::var(m_y))) * (nrow(m_y) - 1)

  m_x <- f_mx(s_k) # conversion matrix m_x
  m_d <- diff(diag(s_k - 1)) %*% diff(diag(s_k))
  m_dtd <- crossprod(m_d)  # matrix D^TD
  m_ybar  <- center_colmeans(m_y)

  m_w <- matrix(0, nrow = 2, ncol = s_k) #constraint matrix
  m_w[1, 1] <- 1
  m_w[2, s_k] <- 1

  iter <- 0

  ## Initialization ##
  if (is.null(init) || missing(init)) {# case 1(principle component)
    init <- f_initialization(m_y, s_k, start_end)
  } else if (is.matrix(init)) {        # case 2(user's theta matrix)
    theta <- init # Set the initialized theta
    init <- project(m_y, theta, diff(theta))
    init$theta <- theta
  } else if (!inherits(init, "ppslp")) {# case 3(supply ppslp object)
    stop("Invalid initial value: should be a matrix or principal_curve")
  }

  # projection step
  proj <- init
  theta <- init$theta

  out <- list(y = m_y,
              theta = theta,
              lambda_l = lambda_l,
              scaled_lam = scaled_lam,
              lambda_s = lambda_s,
              spar = spar,
              dist = proj$dist,
              num_iterations = as.integer(iter),
              call = function_call)

  if (plot_update) {
    pp_plot(out)
  }

  converge <- abs(dist_old - proj$dist) <= thresh * dist_old

  ### outer loop ###
  while (!converge && iter < maxit) {
    iter <- iter + 1

#ptm <- proc.time()# Start the clock!
    # Optimization
    if (is.null(lambda_s)) {
      s_r <- sum(diag(crossprod(proj$m_b))) / sum(diag(m_dtd))
      m_ldtd <- s_r * 256 ^ (3 * spar - 1) * m_dtd
    } else {
      m_ldtd <- lambda_s * m_dtd
    }
    theta <- optimization(m_y, m_ybar, proj$m_b, m_x, m_ldtd,
                          lambda_l, scaled_lam)

    if (!is.null(start_end)) {
      g_matrix <- solve(crossprod(proj$m_b) + m_ldtd)
      inverse_m <- solve(m_w %*% g_matrix %*% t(m_w))
      constraints <- m_w %*% theta - start_end
      theta <- theta - g_matrix %*% crossprod(m_w, inverse_m) %*% constraints
    }
#print(proc.time() - ptm)# Stop the clock

    dist_old  <- proj$dist
    # projection step
    proj <- project(m_y, theta, diff(theta))
    # Evaluation
    converge <- abs(dist_old - proj$dist) <= thresh * dist_old

    # plot the updates
    if (plot_update) {
      out$theta <- theta
      out$dist <- proj$dist
      out$num_iterations <- as.integer(iter)

      pp_plot(out)
    }
  }

  out <- list(y = m_y,
              fitted_y = proj$m_b %*% theta,
              m_b = proj$m_b,
              theta = theta,
              lambda_l = lambda_l,
              scaled_lam = scaled_lam,
              lambda_s = lambda_s,
              spar = spar,
              dist = proj$dist,
              converge = converge,
              num_iterations = as.integer(iter),
              start_end = start_end,
              m_x = m_x,
              m_ldtd = m_ldtd,
              call = function_call)

  class(out) <- "ppslp"
  return(out)
}




#' @rdname f_ppslp
#' @export
#' @importFrom graphics lines
lines.ppslp <- function(pp, ...) {
  graphics::lines(pp$theta[, ], ...)
}

#' @rdname f_ppslp
#' @export
#' @importFrom graphics plot
plot.ppslp <- function(pp, ...) {
  graphics::plot(pp$theta[, ], ..., type = "l")
}

#' @rdname f_ppslp
#' @export
#' @importFrom graphics points
points.ppslp <- function(pp, ...) {
  graphics::points(pp$fitted_y, ...)
}

#' @rdname f_ppslp
#' @param s a parametrized curve, represented by a polygon.
#' @importFrom graphics segments
#' @export
whiskers <- function(x, s, ...) {
  graphics::segments(x[, 1], x[, 2], s[, 1], s[, 2], ...)
}
