#' Plot Principal Curve with Length Penalty
#'
#' This function allows you to fit a principal curve with length penalty.
#' @param pp an prinicipal curve opject includes
#'        \item{Y}: matrix Y
#'        \item{Theta}: spline matrix m_Theta
#'        \item{iter} number of iteration
#'        \item{lambda_l} length parameter
#'        \item{lambda_s} smoother parameter
#'        \item{dist} difference of m_Theta
#'
#' @keywords Principal curve
#' @export
pp_plot <- function(pp) {
  y <- pp$y; theta <- pp$theta; lam_s <- pp$lambda_s; lam_l <- pp$lambda_l
  iter <- pp$num_iterations; dist <- pp$dist; scaled_l <- pp$scaled_lam

  if (is.null(lam_l)) {
    main <- paste("scaled_L=", round(scaled_l, 5))
  } else {
    main <- paste("lam_L=", round(lam_l, 3))
  }

  if (is.null(lam_s)) {
    main <- paste(main, " Spar=", round(pp$spar, 3))
  } else {
    main <- paste(main, " lam_S=", round(lam_s, 3))
  }

  plot(
    y[, 1:2],
    type = "p", pch = 16, cex = 0.8, col = "grey",
    xlim = grDevices::extendrange(y[, 1], f = .2),
    ylim = grDevices::extendrange(y[, 2], f = .2),
    xlab = "", ylab = "",
    sub = paste("Iter:", iter, ", dist=", round(dist, 3)),
    main = main
  )
  lines(theta[, 1:2], col = "blue", lwd = 2)
  points(theta[, 1:2], col = "black", cex = 1, pch = 16)
}





#' Principal Curve with Length Penalty with animation
#'
#' Wrapper function allows you to visualize the fitting procedure of a principal curve with length penalty.
#' @param m_Y your data matrix Y
#' @param s_K number of knots
#' @param s_q dimension of the data, default to be \code{ncol(m_Y)}
#' @param pen penalty term
#' @param increment increment step of length penalty parameter, default to be 0.1
#' @param maxit maximum number of iterations, default to be 60
#' @param s_lambdap smoother parameter, default to be 1
#' @param init_w initialization weights of the first principal component, default to be 6
#' @param init_Theta manual start point for initialization, default to be \code{NULL}, so initializae with first principal component curve
#' @param tol tolerance or stopping criteria, if change in coefficient is less than this value the algorithm will terminate
#' @import splines
#' @return \code{m_Theta} output from the last equation, plot the principal curve
#' @keywords Principal curve
#' @export
f_Animeppslp <- function(m_Y, s_K, s_q, pen = 1,
                         increment = 0.1, maxit = 60, s_lambdap = 2,
                         init_w = 6, init_Theta = NULL, tol = 1e-3) {
  # dir.create("examples")
  setwd("examples")
  png(file = "foo%02d.png")

  result <- f_ppslp(m_Y, s_K, s_q, s_lambdap,
    pen, increment, maxit, init_w,
    init_Theta, tol,
    draw = TRUE
  )
  dev.off()

  system("magick convert -delay 1 foo*.png L0_P100.gif")
  file.remove(list.files(pattern = ".png"))
  setwd("..")
}

# make.mov <- function(){
#   unlink('plot.gif')
#   system('magick convert -delay 1 foo*.png plot.gif')
# }