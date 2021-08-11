#' Optimization
#' @param m_y matrix y
#' @param m_ybar matrix mean Y
#' @param m_b matrix B
#' @param m_x matrix X
#' @param m_ldtd matrix lam_p * D^TD
#' @param lambda_l length parameter
#' @param scaled_lambda scaled length parameter
#' @param max_lam indicator whether return the lenght parameter
#'
#' @return list \code{m_Theta}, status indicator for if it collaspe
optimization <- function(m_y, m_ybar, m_b, m_x, m_ldtd,
                        lambda_l, scaled_lambda) {

    # Block coordinate descent
    m_bbar  <- center_colmeans(m_b)
    m_beta  <- bcgd(m_ybar, m_bbar, m_x, m_ldtd,
                    lambda_l, scaled_lambda, maxit = 1e8)
    theta <- f_Beta2Theta(m_y, m_b, m_beta, m_x)
    return(theta);
}
