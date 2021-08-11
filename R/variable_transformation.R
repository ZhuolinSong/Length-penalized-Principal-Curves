
#' Conversion matrix
#'
#' formulate m_x,
#'
#' @param s_k number of knots
#' @return \code{m_x}
#'
#' @keywords Variable transformation internal
f_mx <- function(s_k) {
  # form m_x
  m_x <- matrix(0L, nrow = s_k, ncol = (s_k - 1))
  m_x[lower.tri(m_x)] <- 1L

  return(m_x)
}


#' Inverse variable transformation
#'
#' Change of variables from Beta to Theta, using optimal <U+03B3>
#'
#' @param m_beta matrix m_beta
#' @param m_y matrix m_y
#' @param m_b matrix m_b
#' @param m_x matrix conversion X
#' @return \code{m_theta} output m_theta
#' @keywords Variable transformation
#' @noRd
#' @examples
#' m_y <- matrix(c(13, -4, 2, -4, 11, -2, 2, -2, 8, 10), 5, 2, byrow = TRUE)
#' m_b <- matrix(c(0, 0, 0, 0.27, 0.72, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.28, 0.72, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow = 5, ncol = 6, byrow = TRUE)
#' m_theta <- matrix(c(4.92, -2.67, 5.83, -1.76, 6.74, -0.85, 7.65, 0.054, 8.56, 0.96, 9.47, 1.87), 6, 2, byrow = TRUE)
#' m_beta <- diff(m_theta)
#' f_Beta2Theta(m_y, m_b, m_beta)
f_Beta2Theta <- function(m_y, m_b, m_beta, m_x) {
  # para setup
  s_k <- ncol(m_b)

  # calculate optimal gamma
  opt_gamma <- colMeans(m_y - m_b %*% m_x %*% m_beta)
  # Transform back to m_theta
  m_theta <- tcrossprod(rep(1, s_k), opt_gamma) + m_x %*% m_beta

  return(m_theta)
}