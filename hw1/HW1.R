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

range(CAtemp_with_residuals$residuals)
breaks <- -7:7
ploteqc(CAtemp_with_residuals, CAtemp_with_residuals$residuals, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Residual plot")