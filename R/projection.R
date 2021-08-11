#' Projection
#'
#' Compute scalar p_i for each y_i.
#' @param y data matrix y
#' @param theta weight matrix theta
#' @param m_Beta transform weight matrix m_Beta
#' @return A structure represent the projected points, including
#'         \code{m_B} the projected spline matrix
#'         \item{v_pi}{The projection index for each point}
#'         \item{dist}{The total squared distance from the curve}
#'         \item{dist_ind}{The squared distances from the curve to each of the respective points}

#'
#' @keywords Projection, Spline
#' @noRd
#' @examples
#' y <- matrix(c(13, -4, 2, -4, 11, -2, 2, -2, 8, 10), 5, 2, byrow = TRUE)
#' theta <- matrix(c(4.92, -2.67, 5.83, -1.76, 6.74, -0.85,7.65, 0.054, 8.56, 0.96, 9.47, 1.87), 6, 2, byrow = TRUE)
#' m_Beta <- diff(theta)
#' out <- f_projection(y, theta, m_Beta)
#' proj_p <- out$m_B %*% theta
#' plot(y[, 1], y[, 2], type = "o")
#' lines(theta[, 1], theta[, 2], col = "blue")
#' points(proj_p[1, 1], proj_p[1, 2], col = "red")
#' points(proj_p[2, 1], proj_p[2, 2], col = "red")
#' points(proj_p[3, 1], proj_p[3, 2], col = "red")
#' points(proj_p[4, 1], proj_p[4, 2], col = "red")
#' points(proj_p[5, 1], proj_p[5, 2], col = "red")
f_projection <- function(y, theta, m_Beta) {
  # Set Paras
  s_n <- nrow(y)
  s_k <- nrow(theta)
  s_q <- ncol(y)
  s_h <- 1 / (s_k - 1)



  #Precompute values
  v_knots <- s_h * c(-1:s_k) # vector of knots
  v_tk <- v_knots[2:(s_k)]
  denorm <- rowSums(m_Beta^2)
  denorm[which(denorm == 0)] <- 1

# ptm <- proc.time()# Start the clock!


  v_proj <- sapply(c(1:s_n), outer_loop <- function(i) {
    y_i <- y[i, ]
    numera <- s_h * diag(m_Beta %*% (y_i - t(theta)))
    thresh <- numera / denorm

    v_pi <- thresh + v_tk # a_k
    v_pi[which(thresh <= 0)] <- v_tk[which(thresh <= 0)] # lower value
    v_pi[which(thresh >= s_h)] <- v_tk[which(thresh >= s_h)] + s_h # upper value

    # #Inner Loop
    # for (k in 1:(k-1)) {#loop over each piece of linear line
    #   t_k = (k-1)*s_h;
    #
    #   if(f_2nrm(m_Beta[k,])!=0) {# Check if the piece of line collaspe
    #       #Calculate a_k
    #       temp = crossprod(m_Beta[k,]);
    #       a_k = crossprod(m_Beta[k,], y_i-theta[k,]);
    #       a_k = a_k * s_h/temp + t_k
    #
    #       #Project a_k back into [t_k, t_(k+1))
    #       if (t_k <= a_k && a_k < (t_k+s_h)) {
    #         v_pi[k] = a_k;
    #       } else if ( abs(a_k-t_k) > abs(a_k - (t_k + s_h))) {
    #           v_pi[k] = t_k+s_h;
    #       } else {
    #           v_pi[k] = t_k;
    #       }
    #     }else{v_pi[k]= t_k;}
    #
    # }
    # check
    if (sum(v_pi < 0 | v_pi > 1)) print(0)

    # Calculate the distance
    b_pi <- t(splines::splineDesign(v_knots, v_pi, 2))
    dis <- colSums((y_i - crossprod(theta, b_pi))^2) # previously, diag(crossprod(y_i-crossprod(theta, b_pi)))
    return(c(v_pi[which.min(dis)], min(dis))) # find the optimal projected point
  })



# print(proc.time() - ptm)# Stop the clock



  ## prepare output
  v_pi <- v_proj[1, ]
  dist_ind <- v_proj[2, ]
  dist <- sum(dist_ind)
  m_B <- splines::splineDesign(knots = v_knots, x = v_pi, ord = 2)

  out <- list(
    m_B = m_B,
    #v_pi = v_pi, dist_ind = dist_ind,
    dist = dist
  )

  return(out)
}