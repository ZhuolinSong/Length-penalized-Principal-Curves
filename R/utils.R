#' The 2norm or Frobenius norm function
#'
#'
#' @param x vector or matrix to calculate its 2norm or Frobenius norm
#' @return 2norm or Frobenius norm
#'
#' @keywords norm internal
f_2nrm <- function(x) sqrt(sum(x^2))
#' function from GLODE
#'
#' prepare some functions(from GLODE code)
#' @param m matrix
#'
#' @keywords norm internal
as.doublematrix <- function(m) {
  matrix(as.double(m), nrow(m), ncol(m))
}


#' function
#'
#' center the col of matrix x
#' @param x a matrix
#'
#' @keywords norm internal
center_colmeans <- function(x) {
  xcenter <- colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}