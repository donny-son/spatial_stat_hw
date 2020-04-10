# Objective : TODO::Write multivariate normal creating function
# Created by: ericson
# Created on: 2020/04/10


# Q1 - (b)
library(mvtnorm)
generate_multivariate_normal <- function (mu, Sigma, seed= 1) {
  set.seed(seed)
  lower_triangle_matrix <- t(chol(Sigma))
  dimension <- nrow(Sigma)
  Z <- rmvnorm(n = 1, mean = rep(0, dimension), sigma = diag(dimension))
  Y = mu + lower_triangle_matrix %*% t(Z)
  return(Y)
}

# Q1 - (c)
library(mvtnorm)
library(clusterGeneration)
library(yarrr)
generate_z <- function(dimension, seed= 1) {
  set.seed(seed)
  Z <- rmvnorm(n = 1, mean = rep(0, dimension), sigma = diag(dimension))
  return(Z)
}

generate_multivariate_normal_with_Z <- function (mu, Sigma, Z, seed=1) {
  set.seed(seed)
  lower_triangle_matrix <- t(chol(Sigma))
  Y = mu + (lower_triangle_matrix %*% t(Z))
  return(Y)
}


dimension <- 100
Z_sample <- generate_z(dimension)

plot(x=1,
     type="n",
     xlim = c(1,100),
     ylim = c(-10,10),
     pch = 16,
     xlab="Sample Index",
     ylab="Y_sample",
     )
grid()

set.seed(1)
Sigma <- genPositiveDefMat(dimension, covMethod="eigen")$Sigma
Y_sample <- generate_multivariate_normal_with_Z(mu = rep(0,100), Sigma = Sigma, Z = Z_sample)
points(1:100, Y_sample,
       pch = 16,
       col = transparent("coral2", trans.val = .8))

set.seed(1)
Sigma <- genPositiveDefMat(dimension, covMethod="onion")$Sigma
Y_sample <- generate_multivariate_normal_with_Z(mu = rep(0,100), Sigma = Sigma, Z = Z_sample)
points(1:100, Y_sample,
       pch = 16,
       col = transparent("coral", trans.val = .5))

set.seed(1)
Sigma <- genPositiveDefMat(dimension, covMethod="unifcorrmat")$Sigma
Y_sample <- generate_multivariate_normal_with_Z(mu = rep(0,100), Sigma = Sigma, Z = Z_sample)
points(1:100, Y_sample,
       pch = 16,
       col = transparent("coral3", trans.val = .3))

Sigma <- diag(dimension)
Y_sample <- generate_multivariate_normal_with_Z(mu = rep(0,100), Sigma = Sigma, Z = Z_sample)
points(1:100, Y_sample,
       pch = 12,
       col = transparent("blue", trans.val = .1))