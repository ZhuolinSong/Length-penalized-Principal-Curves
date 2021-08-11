devtools::load_all()
library(ppclp) # Huanchen's package
library(princurve) # Hastie's package
library(splines)
library(splines2)
library(tidyverse) # required by ppclp

#Load spect data into R from the package ppclp
data("spectExample")
y <- cbind(spectExample$x, spectExample$y, spectExample$z)
y <- y[1:100, ]

# Create time storage
time.table <- matrix(ncol = 5, nrow = 3)

# Compute the principal curve
ptm <- proc.time() # Start the clock!
fitlength <- f_ppslp(y, s_k = 200, maxit = 1e4, thresh = 1e-3,
    lambda_s = 1e6, lambda_l = 5, plot_update = F
    )
time.table[1, ] <- proc.time() - ptm # Stop the clock


# HuanChen's method
ptm <- proc.time() # Start the clock!
tmpCurve <- ppclp3D(y[, 1], y[, 2], y[, 3],
            spectExample$xFix, spectExample$yFix, spectExample$zFix)

time.table[2, ] <- proc.time() - ptm # Stop the clock

# Hastie's method
ptm <- proc.time() # Start the clock!
fit <- principal_curve(y, thresh = 1e-3, maxit = 100)
time.table[3, ] <- proc.time() - ptm # Stop the clock


save(time.table, fitlength, tmpCurve, fit, file = "colon.RData")