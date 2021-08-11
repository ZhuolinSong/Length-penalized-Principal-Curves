#' Cross Validation Grid search
#' on the loss
#'
#' @param l_range scaled length range(from 0 to 1):
#'              either a vector of length parameters
#'              or an integer indicating the range of the selection range
#'              in exponential:
#'                         exp(-range:-1)
#' @param s_range spar range(from 0 to 1):
#'              either a vector of smoothing parameters
#'              or an integer indicating the range of the selection range
#'              : (1 - 0.3) / range * (0:range) + 0.3
#' @param folds the number of fold to split Y, default to be 5
#' @param m_y your data matrix Y
#' @param s_k number of knots
#' @param lambda_l length parameter(from 0 to inf), default to be NULL
#'                either a vector of length parameters
#'                or an integer indicating the range of the selection range
#'                in exponential:
#'                         exp(-range:range)
#' @param ... other args for f_ppslp(maxit, thresh, init)
#' @import splines
#' @return a fitted class ppslp using the selected parameters
#' @keywords Principal curve
#' @export
cv_grid <- function(m_y, s_k,
                    l_range = 0,
                    s_range = 10,
                    folds = 10,
                    lambda_l = NULL, ...) {
  s_n <- nrow(m_y) # Define the number of observation
  if (s_n < folds) {
    stop(paste("Can't do K-folds cv, your folds:", folds, "is greater than n:", s_n))
  }

  if (length(s_range) == 1) {# define the range of smoothing parameter
    v_smoothing <- (1 - 0.3) / s_range * (0:s_range) + 0.3
  } else {
    v_smoothing <- s_range
  }

  if (is.null(lambda_l)) {
    if (length(l_range) == 1) {# define the range of scaled length parameter
      v_length <- c(0, exp(-l_range:-1))
    } else {
      v_length <- l_range
    }
  } else {
    if (length(lambda_l) == 1) {# define the range of length parameter
      v_length <- c(0, exp(-lambda_l:lambda_l))
    } else {
      v_length <- lambda_l
    }
  }

  loss <- matrix(0, nrow = length(v_length), ncol = length(v_smoothing))


  m_y <- m_y[sample(s_n), ] # Randomly shuffle the data
  v_folds <- cut(seq(s_n), breaks = folds, labels = FALSE) # Create k folds

  for (k in 1:folds) { # Perform k fold cross validation
    test_idx <- which(v_folds == k, arr.ind = TRUE) # Segement data by fold
    test <- m_y[test_idx, ]
    train <- m_y[-test_idx, ]

    # Loop the range
    for (i in seq_along(v_length)) {
      len <- v_length[i]
      for (j in seq_along(v_smoothing)) {
        spar <- v_smoothing[j]
        # Build the model
        if (is.null(lambda_l)) {
          pp <- f_ppslp(train, s_k, spar, len, ...)
        } else {
          pp <- f_ppslp(train, s_k, spar, lambda_l = len, ...)
        }
        # Store the loss
        loss[i, j] <- loss[i, j] + (orig_criteria(pp, test) / folds)
      }
    }
  }

  cv_idx <- arrayInd(which.min(loss), dim(loss))
  len <- v_length[cv_idx[1]]
  lambda_s <- v_smoothing[cv_idx[2]]

  if (is.null(lambda_l)) {
    pp <- f_ppslp(m_y, s_k, lambda_s, len, ...)
  } else {
    pp <- f_ppslp(m_y, s_k, lambda_s, lambda_l = len, ...)
  }

  pp$loss <- loss
  return(pp)
}


#' Cross Validation Tune Smoothing parameter
#' on loss
#'
#' @param l_range scaled length range(from 0 to 1):
#'               a vector of length parameters
#' @param s_range smoothing range(from 0 to inf):
#'              either a vector of smoothing parameters
#'              or an integer indicating the range of the selection range
#'              in exponential:
#'                         exp(-range:range)
#' @param folds the number of fold to split Y, default to be 5
#' @param m_y your data matrix Y
#' @param s_k number of knots
#' @param lambda_l length parameter(from 0 to inf), default to be NULL
#'                either a vector of length parameters
#'                or an integer indicating the range of the selection range
#'                in exponential:
#'                         exp(-range:range)
#' @param ... other args for f_ppslp(maxit, thresh, init)
#' @import splines
#' @return the selected smoothing parameters
#' @keywords Principal curve
#' @export
cv_smooth <- function(m_y, s_k,
                    l_range = 0,
                    s_range = 10,
                    folds = 10,
                    lambda_l = NULL, ...) {
  s_n <- nrow(m_y) # Define the number of observation
  if (s_n < folds) {
    stop(paste("Can't do K-folds cv, your folds:", folds, "is greater than n:", s_n))
  }

  if (length(s_range) == 1) {# define the range of smoothing parameter
    v_smoothing <- (1 - 0.3) / s_range * (0:s_range) + 0.3
  } else {
    v_smoothing <- s_range
  }

  if (is.null(lambda_l)) {
    v_length <- l_range
  } else {
    v_length <- lambda_l
  }

  loss <- matrix(0, nrow = length(v_length), ncol = length(v_smoothing))


  m_y <- m_y[sample(s_n), ] # Randomly shuffle the data
  v_folds <- cut(seq(s_n), breaks = folds, labels = FALSE) # Create k folds

  for (k in 1:folds) { # Perform k fold cross validation
    test_idx <- which(v_folds == k, arr.ind = TRUE) # Segement data by fold
    test <- m_y[test_idx, ]
    train <- m_y[-test_idx, ]

    # Loop the range
    for (i in seq_along(v_length)) {
      len <- v_length[i]
      for (j in seq_along(v_smoothing)) {
        lambda_s <- v_smoothing[j]
        # Build the model
        if (is.null(lambda_l)) {
          pp <- f_ppslp(train, s_k, lambda_s, len, ...)
        } else {
          pp <- f_ppslp(train, s_k, lambda_s, lambda_l = len, ...)
        }
        # Store the loss
        loss[i, j] <- loss[i, j] + orig_criteria(pp, test) / folds
      }
    }
  }

  cv_idx <- apply(loss, 1, which.min)
  lambda_s <- v_smoothing[cv_idx]
  return(lambda_s)
}


#' Cross Validation Tune length parameter
#' on the loss (need to use cv_smooth)
#'
#' @param l_range scaled length range(from 0 to 1):
#'              either a vector of length parameters
#'              or an integer indicating the range of the selection range
#'              in exponential:
#'                         exp(-range:-4)
#' @param s_range smoothing range(from 0 to inf):
#'              either a vector of smoothing parameters
#' @param folds the number of fold to split Y, default to be 5
#' @param m_y your data matrix Y
#' @param s_k number of knots
#' @param lambda_l length parameter(from 0 to inf), default to be NULL
#'                either a vector of length parameters
#'                or an integer indicating the range of the selection range
#'                in exponential:
#'                         exp(-range:range)
#' @param ... other args for f_ppslp(maxit, thresh, init)
#' @import splines
#' @return a fitted class ppslp using the selected parameters
#' @keywords Principal curve
#' @export
cv_length <- function(m_y, s_k,
                    l_range = 10,
                    s_range = 0.6,
                    folds = 10,
                    lambda_l = NULL, ...) {
  s_n <- nrow(m_y) # Define the number of observation
  if (s_n < folds) {
    stop(paste("Can't do K-folds cv, your folds:", folds, "is greater than n:", s_n))
  }

  if (is.null(lambda_l)) {
    if (length(l_range) == 1) {# define the range of scaled length parameter
      v_length <- c(0, exp(-l_range:-1))
    } else {
      v_length <- l_range
    }
  } else {
    if (length(lambda_l) == 1) {# define the range of length parameter
      v_length <- c(0, exp(-lambda_l:lambda_l))
    } else {
      v_length <- lambda_l
    }
  }

  v_smoothing <- s_range

  loss <- matrix(0, nrow = length(v_length), ncol = length(v_smoothing))

  m_y <- m_y[sample(s_n), ] # Randomly shuffle the data
  v_folds <- cut(seq(s_n), breaks = folds, labels = FALSE) # Create k folds

  for (k in 1:folds) { # Perform k fold cross validation
    test_idx <- which(v_folds == k, arr.ind = TRUE) # Segement data by fold
    test <- m_y[test_idx, ]
    train <- m_y[-test_idx, ]

    # Loop the range
    for (i in seq_along(v_length)) {
      len <- v_length[i]
      for (j in seq_along(v_smoothing)) {
        lambda_s <- v_smoothing[j]
        # Build the model
        if (is.null(lambda_l)) {
          pp <- f_ppslp(train, s_k, lambda_s, len, ...)
        } else {
          pp <- f_ppslp(train, s_k, lambda_s, lambda_l = len, ...)
        }
        # Store the loss
        loss[i, j] <- loss[i, j] + len_criteria(pp, test) / folds
      }
    }
  }

  cv_idx <- apply(loss, 2, which.min)
  len <- v_length[cv_idx]
  return(len)
}


#' Cross Validation Princurve
#' on loss
#'
#' @param s_range smoothing range(from 0 to inf):
#'              either a vector of smoothing parameters
#'              or an integer indicating the range of the selection range
#'              in exponential:
#'                         1 / s_range * (1:s_range)
#' @param folds the number of fold to split Y, default to be 5
#' @param m_y your data matrix Y
#' @param s_k number of knots
#' @param ... other args for f_ppslp(maxit, thresh, init)
#' @import splines
#' @import princurve
#' @return the selected smoothing parameters
#' @keywords Principal curve
#' @export
cv_princurve <- function(m_y, s_k,
                    s_range = 10,
                    folds = 10,
                    ...) {
  s_n <- nrow(m_y) # Define the number of observation
  if (s_n < folds) {
    stop(paste("Can't do K-folds cv, your folds:", folds, "is greater than n:", s_n))
  }


  if (length(s_range) == 1) {# define the range of smoothing parameter
    v_smoothing <- 1 / s_range * (1:s_range)
  } else {
    v_smoothing <- s_range
  }

  loss <- matrix(0, nrow = length(v_smoothing))


  m_y <- m_y[sample(s_n), ] # Randomly shuffle the data
  v_folds <- cut(seq(s_n), breaks = folds, labels = FALSE) # Create k folds

  for (k in 1:folds) { # Perform k fold cross validation
    test_idx <- which(v_folds == k, arr.ind = TRUE) # Segement data by fold
    test <- m_y[test_idx, ]
    train <- m_y[-test_idx, ]

    # Loop the range
    for (i in seq_along(v_smoothing)) {
        spar <- v_smoothing[i]
        # Build the model
        fit <- princurve::principal_curve(m_y, spar = spar, ...)
        # Store the loss
        fit$theta <- fit$s
        loss[i] <- loss[i] + orig_criteria(fit, test) / folds
    }

  }

  cv_idx <- which.min(loss)
  spar <- v_smoothing[cv_idx]
  fit <- princurve::principal_curve(m_y, spar = spar, ...)
  fit$spar <- spar
  fit$y <- m_y
  return(fit)
}