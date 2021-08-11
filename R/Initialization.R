#' Default Initialization
#'
#' initialize the weight matrix m_Theta
#'
#' @param y data matrix y
#' @param k number of knots
#' @param start_end start and end constraint points(matrix):
#'                  [start, end]^T, default to NULL no constraints
#' @return \code{m_Theta}, output m_Theta
#'
#' @keywords Default initialization, internal
#' @examples
#' \dontrun{
#' y <- matrix(c(13, -4, 2, -4, 11, -2, 2, -2, 8, 10), 5, 2, byrow = TRUE)
#' k <- 6L
#' init <- f_initialization(y, k)
#' plot(y[, 1], y[, 2], type = "o")
#' lines(init$theta[, 1], init$theta[, 2], col = "blue")
#' }
f_initialization <- function(y, k, start_end) {

  if (!is.null(start_end)) {
    slope <- (start_end[1, ] - start_end[2, ]) / (k - 1)
    theta <- outer(seq(k - 1, 0), slope) + start_end[2, ]
    proje <- project(y, theta, diff(theta))
    proje$theta <- theta
    return(proje)
  }
  
  # svd
  y_mean <- colMeans(y)
  ystar <- scale(y, center = y_mean, scale = FALSE)
  svd_y <- svd(ystar)
  dd <- svd_y$d
  v_pi <- svd_y$u[, 1] * dd[1]

  # m_Theta
  slope <- (max(v_pi) - min(v_pi)) / (k - 1)
  intercept <- min(v_pi)
  s <- seq(0, k - 1) * slope + intercept
  theta <- rep(y_mean, rep.int(k, ncol(y))) + outer(s, svd_y$v[, 1])
  # m_b
  v_pi <- (v_pi - min(v_pi)) / (max(v_pi) - min(v_pi))
  v_knots <- 1 / (k - 1) * c(-1:k) # vector of knots

  m_b <- splines::splineDesign(knots = v_knots, x = v_pi, ord = 2, outer.ok = T)

  dist <- sum((dd^2)[-1])

  return(list(theta = theta,
              m_b = m_b,
              dist = dist)
        )
}
