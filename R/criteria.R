#' Evaluation
#'
#' Evaluate the difference b/w coefficient m_Theta
#'
#' @param m_Y data matrix m_Y
#' @param m_Theta matrix m_Theta
#' @return \code{diff}, difference b/w coefficient m_Theta
#'
#' @keywords Evaluation
theta_diff <- function(new_theta, old_theta) {
  return(f_2nrm(new_theta - old_theta))
}


#' Evaluation
#'
#'  Length Penalty
#'
#' @param theta matrix m_Theta
#' @param s_lambda_L length parameter
#' @return Length penalty
#'
#' @keywords Evaluation internal
f_lenpen <- function(theta, lambda_l) {
  beta <- diff(theta)
  return(lambda_l * sum(sqrt(diag(tcrossprod(beta)))))
}

#' Evaluation
#'
#' Evaluate Distance + Length Penalty (use for tunning the smoothing parameter)
#'
#' @param pp an prinicipal curve opject includes
#'        \item{Y}: matrix Y
#'        \item{Theta}: spline matrix m_Theta
#'        \item{iter} number of iteration
#'        \item{lambda_l} length parameter
#'        \item{lambda_s} smoother parameter
#'        \item{dist} difference of m_Theta
#' @param new_data External data matrix
#'
#' @return \code{diff}, difference b/w coefficient m_Theta
#'
#' @keywords Evaluation internal
len_criteria <- function(pp, new_data = NULL) {
  theta <- pp$theta; lambda_l <- pp$lambda_l
  if (is.null(lambda_l)) {
    m_ybar <- center_colmeans(pp$y)
    m_bbar  <- center_colmeans(pp$m_b)
    lambda_l <- pp$scaled_lam * bcgd(m_ybar, m_bbar, pp$m_x, pp$m_ldtd,
                                        lambda_l, pp$scaled_lam, max_lam = T)
  }

  dist <- orig_criteria(pp, new_data)
  return(dist + f_lenpen(theta, lambda_l))
}


#' Evaluation
#'
#' Evaluate Distance
#'
#' @param pp an prinicipal curve opject includes
#'        \item{Y}: matrix Y
#'        \item{Theta}: spline matrix m_Theta
#'        \item{iter} number of iteration
#'        \item{lambda_l} length parameter
#'        \item{lambda_s} smoother parameter
#'        \item{dist} difference of m_Theta
#' @param new_data External data matrix
#'
#' @return \code{diff}, difference b/w coefficient m_Theta
#'
#' @keywords Evaluation internal
orig_criteria <- function(pp, new_data = NULL) {

  if (missing(new_data) || is.null(new_data)) {
    dist <- pp$dist;
  } else {
    theta <- pp$theta;
    out <- project(new_data, theta, diff(theta))
    dist <- out$dist
  }

  return(dist)
}