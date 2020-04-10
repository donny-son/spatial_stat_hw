library(sp)
library(gstat)
library(fields)
library(classInt)
library(maps)

load("hw1/CAtemps.RData")

ploteqc <- function(spobj, z, breaks, ...){
  pal <- tim.colors(length(breaks)-1)
  fb <- classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}

## Plotting
par(mfrow = c(1,2))
range(CAtemp$avgtemp)
breaks <- 40:75
# x11()
ploteqc(CAtemp, CAtemp$avgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")

range(CAgrid$elevation)
breaks <- seq(-100, 3600, by = 100)
# x11()
ploteqc(CAgrid, CAgrid$elevation, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Elevations at prediction locations, m")

##########################
#### Solution for Q2 #####
##########################

# a
CAtemp_with_coordinates <- cbind(CAtemp, coordinates(CAtemp))
ordinary_least_squares <- lm(avgtemp ~ lon + lat + elevation, data = CAtemp_with_coordinates)
predictions <- predict(ordinary_least_squares, newdata = CAtemp_with_coordinates)
residuals <- CAtemp_with_coordinates$avgtemp - predictions
CAtemp_with_residuals <- cbind(CAtemp_with_coordinates, residuals); names(CAtemp_with_residuals)[5] <- 'residuals'
CAtemp_with_predictions <- cbind(CAtemp_with_coordinates, predictions); names(CAtemp_with_predictions)[5] <- 'predictions'

range(CAtemp_with_residuals$residuals)
breaks <- -7:7
ploteqc(CAtemp_with_residuals, CAtemp_with_residuals$residuals, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Residual plot")

par(mfrow = c(1,2))
range(CAtemp$avgtemp)
breaks <- 40:75
ploteqc(CAtemp, CAtemp$avgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")
ploteqc(CAtemp_with_predictions, CAtemp_with_predictions$predictions, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Prediction plot")

# b
variogram <- variogram(residuals ~ 1, data = CAtemp_with_coordinates)
variogram.fit <- fit.variogram(variogram, vgm(1, "Exp", 100, 2))
plot(variogram, variogram.fit)

sigma_sqrd_hat <- variogram.fit$psill[2]
tau_sqrd_hat <- variogram.fit$psill[1]
rho_hat <- variogram.fit$range[2]

# c
# Get distance matrix.
distance_matrix <- rdist(coordinates(CAtemp_with_coordinates))
dim(distance_matrix) # (200, 200)

# Create cov matrix.
exponential_covariance <- function(distance_matrix, tau_sqrd_hat, sigma_sqrd_hat, rho_hat){
  n = dim(distance_matrix)[1]
  matrix_with_no_nugget <- matrix(rep(0, n*n), ncol=n)
  print(dim(matrix_with_no_nugget))
  for (i in 1:n){
    for (j in 1:n){
      h = distance_matrix[i,j]
      matrix_with_no_nugget[i,j] <- sigma_sqrd_hat * exp(-h/rho_hat)
    }
  }
  matrix_with_nugget <- (sigma_sqrd_hat + tau_sqrd_hat) * diag(n)
  print(dim(matrix_with_nugget))
  return(matrix_with_no_nugget + matrix_with_nugget)
}
covariance_matrix_hat <- exponential_covariance(distance_matrix, tau_sqrd_hat, sigma_sqrd_hat, rho_hat)

# Invert covariance matrix and store.
inverse_covariance_matrix_hat <- solve(covariance_matrix_hat)

# Create X
X <- cbind(CAtemp$elevation, coordinates(CAtemp)); colnames(X) <- c('elevation', 'lon', 'lat')
y <- CAtemp$avgtemp

# Form beta_gls
(beta_gls <- solve(t(X) %*% inverse_covariance_matrix_hat %*% X) %*% t(X)%*% inverse_covariance_matrix_hat %*%y)
