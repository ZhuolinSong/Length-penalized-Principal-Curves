#' Block Coordinate descence
#'
#'
#' An group lasso solver
#'
#'
#' @param m_ybar matrix Y
#' @param m_bbar matrix B
#' @param m_x matrix X
#' @param m_ldtd matrix D
#' @param lambda_l length parameter
#' @param scaled_lambda scaled length parameter
#' @param maxit maximum number of iterations, default to be 60
#' @param max_lam indicator whether return the lenght parameter
#'
#' @return \code{m_beta} matrix Beta
#' @keywords group lasso
#'
#'
#' @examples
#' library(ppslp)
#' m_y <- matrix(c(13, -4, 2, -4, 11, -2, 2, -2, 8, 10), 5, 2, byrow = TRUE)
#' s_k     <- 5L;
#' m_ybar  <- center_colmeans(m_y)
#' # initialize m_Theta
#' m_Theta <- f_initialization(m_y, s_k, w=6)
#' m_beta  <- diff(m_Theta) # Transform m_Theta to m_beta
#' # projection step
#' m_b     <- f_projection(m_y, m_Theta, m_beta)
#' m_bbar  <- center_colmeans(m_b);
#' lambda_s <- 1;
#' lambda_l <- 1;
#' bcgd(m_ybar, m_bbar, lambda_l, lambda_s)

bcgd <- function(m_ybar, m_bbar, m_x, m_ldtd, lambda_l, scaled_lambda,
                maxit=1e8, max_lam=F) {

  # Get column# and row# from the inputs
  s_k <- ncol(m_bbar)
  s_q <- ncol(m_ybar)


  # Initialize
  v_active <- rep(0L, s_k - 1) # 0 indicates inactive
  m_beta <- matrix(0L, nrow = (s_k - 1), ncol = s_q)
  # Pre-calculate
  m_c <- crossprod(m_bbar %*% m_x, m_ybar)
  m_h <- crossprod(m_bbar) + m_ldtd
  m_o <- crossprod(m_x, m_h %*% m_x)
  v_a <- diag(m_o)



  if (is.null(lambda_l)) {
    lambda_l <- scaled_lambda * sqrt(max(rowSums(m_c^2)))
    if (max_lam) {
      return(lambda_l)
    }
  }

  if (lambda_l == 0) {
    ls_beta <- solve(m_o, m_c)
    return(ls_beta)
  }

  # Outer loop
  for (i in 1:maxit) {

    # Inner loop
    for (k in 1:(s_k - 1)) {
      if (v_active[k] == 1L) {
        m_bnk <- m_beta[-k, ]
        v_s <- m_c[k, ] - crossprod(m_o[k, -k], m_bnk)
        s_a <- v_a[k]
        s_coef <-  1 - lambda_l / f_2nrm(v_s)

        if (s_coef > 0) {
          m_beta[k, ] <- 1 / s_a * s_coef * v_s
        }
        else { # set the beta_k as zero and remove from active set
          m_beta[k, ] <- 0L
          v_active[k] <- 0L
        }
      }
    }


    #Check KKT condition:
    v_kkt <- rep(0L, s_k - 1)
    m_s <- m_c - m_o %*% m_beta
    for (k in 1:(s_k - 1)) {
      if (v_active[k] == 0) {
        v_kkt[k] <- f_2nrm(m_s[k, ])
      }
    }

    # Add a new beta_k into the active set
    if (max(v_kkt) > lambda_l) {
      v_active[which.max(v_kkt)] <- 1L;
    } else {
      break;
    }

  }

  return(m_beta)

}


# m_y <- matrix(c(13, -4, 2, -4, 11, -2, 2, -2, 8, 10), 5, 2, byrow = TRUE)
# s_k     <- 5L;
# m_ybar  <- center_colmeans(m_y)
# # initialize m_Theta
# m_Theta <- f_initialization(m_y, s_k, w=6)
# m_beta  <- diff(m_Theta) # Transform m_Theta to m_beta
# # projection step
# m_b     <- f_projection(m_y, m_Theta, m_beta)
# m_bbar  <- center_colmeans(m_b);
# lambda_s <- 1;
# lambda_l <- 1;


# bcgd(m_ybar, m_bbar, lambda_l, lambda_s)
