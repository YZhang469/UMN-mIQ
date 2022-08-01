library(radiant.data)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(gtable)
library(lemon)
library(scales)

source("functions.R")

###############################
## data generative mechanism ##
###############################

# notation in the manuscript: C := c2, gamma := c1

sigma <- 1

# stage 1 heterogeneous effects
calcOut.homoA2 <- function(X1, A1, X2, A2, C){
  H2 <- cbind.data.frame("intercept" = rep(1, length(X1)), "X1" = X1, "A1" = A1, "X2" = X2)
  beta20 <- c(3, -1, 0.1, -0.1)
  beta21 <- C * c(-6, 0, 5, 0)
  Y.exp <- as.matrix(H2)%*%beta20 + A2*(as.matrix(H2)%*%beta21)
  return(Y.exp)
}
calcOut.heteroA2 <- function(X1, A1, X2, A2, C){
  H2 <- cbind.data.frame("intercept" = rep(1, length(X1)), "X1" = X1, "A1" = A1, "X2" = X2)
  beta20 <- c(3, -1, 0.1, -0.1)
  beta21 <- C * c(-6, -4, 5, -0.2)
  Y.exp <- as.matrix(H2)%*%beta20 + A2*(as.matrix(H2)%*%beta21)
  return(Y.exp)
}

# function to generate a training dataset
generateTrainData.uncorUV <- function(n = 250, calcOut = calcOut.hetero, C = 1.5, gamma = 0){
  X1 <- rnorm(n, -2, 1)
  A1 <- 2 * rbinom(n, 1, 0.5) - 1
  A2 <- 2 * rbinom(n, 1, 0.5) - 1
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma) # sd(phi) is too small, 7-12 is appropriate
  # unmeasured variable uncorrelated with stage 2 main effects predictors
  V <- rnorm(n, -1, 1)
  Y <- calcOut(X1, A1, X2, A2, C) + gamma*V + phi
  train <- data.frame("ID" = 1:n, "X1" = X1, "A1" = A1, "X2" = X2, "A2" = A2, "Y" = Y)
  return(train)
}
generateTrainData.corUV <- function(n = 250, calcOut = calcOut.hetero, C = 1.5, gamma = 0){
  X1 <- rnorm(n, -2, 1)
  A1 <- 2 * rbinom(n, 1, 0.5) - 1
  A2 <- 2 * rbinom(n, 1, 0.5) - 1
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma)
  # unmeasured variable correlated with stage 2 main effects predictors
  V <- rnorm(n, 2*X1*X2, 1)
  Y <- calcOut(X1, A1, X2, A2, C) + gamma*V + phi
  train <- data.frame("ID" = 1:n, "X1" = X1, "A1" = A1, "X2" = X2, "A2" = A2, "Y" = Y)
  return(train)
}
generateTrainData.heteroA1 <- function(n = 250, calcOut = calcOut.hetero, C = 1.5, gamma = 0){
  X1 <- rnorm(n, -2, 1)
  A1 <- 2 * rbinom(n, 1, 0.5) - 1
  A2 <- 2 * rbinom(n, 1, 0.5) - 1
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma)
  Y <- calcOut(X1, A1, X2, A2, C) + gamma*X1*A1 + phi
  train <- data.frame("ID" = 1:n, "X1" = X1, "A1" = A1, "X2" = X2, "A2" = A2, "Y" = Y)
  return(train)
}

# function to generate a test dataset with all potential outcomes
generateTestData.uncorUV <- function(n = 10000, calcOut = calcOut.hetero, C = 1.5, gamma = 0){
  X1 <- rnorm(n, -2, 1)
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma)
  V <- rnorm(n, -1, 1)
  Y.1.1 <- calcOut(X1, A1 = rep(1, n), X2, A2 = rep(1, n), C) + gamma*V + phi
  Y.1.n1 <- calcOut(X1, A1 = rep(1, n), X2, A2 = rep(-1, n), C) + gamma*V + phi
  Y.n1.1 <- calcOut(X1, A1 = rep(-1, n), X2, A2 = rep(1, n), C) + gamma*V + phi
  Y.n1.n1 <- calcOut(X1, A1 = rep(-1, n), X2, A2 = rep(-1, n), C) + gamma*V + phi
  test <- data.frame("ID" = 1:n, "X1" = X1, "X2" = X2, 
                     "Y(1,1)" = Y.1.1, "Y(1,-1)" = Y.1.n1, 
                     "Y(-1,1)" = Y.n1.1, "Y(-1,-1)" = Y.n1.n1, 
                     "Yopt" = pmin(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1, na.rm = TRUE), 
                     check.names = FALSE)
  trt.ind <- which.pmin(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1, na.rm = TRUE)
  # optimal treatments
  test$A1opt <- ifelse(trt.ind == 1|trt.ind == 2, 1, -1)
  test$A2opt <- ifelse(trt.ind == 1|trt.ind == 3, 1, -1)
  return(test)
}
generateTestData.corUV <- function(n = 10000, calcOut = calcOut.hetero, C = 1.5, gamma = 0){
  X1 <- rnorm(n, -2, 1)
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma)
  V <- rnorm(n, 2*X1*X2, 1)
  Y.1.1 <- calcOut(X1, A1 = rep(1, n), X2, A2 = rep(1, n), C) + gamma*V + phi
  Y.1.n1 <- calcOut(X1, A1 = rep(1, n), X2, A2 = rep(-1, n), C) + gamma*V + phi
  Y.n1.1 <- calcOut(X1, A1 = rep(-1, n), X2, A2 = rep(1, n), C) + gamma*V + phi
  Y.n1.n1 <- calcOut(X1, A1 = rep(-1, n), X2, A2 = rep(-1, n), C) + gamma*V + phi
  test <- data.frame("ID" = 1:n, "X1" = X1, "X2" = X2, 
                     "Y(1,1)" = Y.1.1, "Y(1,-1)" = Y.1.n1, 
                     "Y(-1,1)" = Y.n1.1, "Y(-1,-1)" = Y.n1.n1, 
                     "Yopt" = pmax(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1, na.rm = TRUE), 
                     check.names = FALSE)
  trt.ind <- which.pmin(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1, na.rm = TRUE)
  # optimal treatments
  test$A1opt <- ifelse(trt.ind == 1|trt.ind == 2, 1, -1)
  test$A2opt <- ifelse(trt.ind == 1|trt.ind == 3, 1, -1)
  return(test)
}
generateTestData.heteroA1 <- function(n = 10000, calcOut = calcOut.hetero, C = 1.5, gamma = 0){
  X1 <- rnorm(n, -2, 1)
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma)
  Y.1.1 <- calcOut(X1, A1 = rep(1, n), X2, A2 = rep(1, n), C) + gamma*X1 + phi
  Y.1.n1 <- calcOut(X1, A1 = rep(1, n), X2, A2 = rep(-1, n), C) + gamma*X1 + phi
  Y.n1.1 <- calcOut(X1, A1 = rep(-1, n), X2, A2 = rep(1, n), C) - gamma*X1 + phi
  Y.n1.n1 <- calcOut(X1, A1 = rep(-1, n), X2, A2 = rep(-1, n), C) - gamma*X1 + phi
  test <- data.frame("ID" = 1:n, "X1" = X1, "X2" = X2, 
                     "Y(1,1)" = Y.1.1, "Y(1,-1)" = Y.1.n1, 
                     "Y(-1,1)" = Y.n1.1, "Y(-1,-1)" = Y.n1.n1, 
                     "Yopt" = pmax(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1, na.rm = TRUE), 
                     check.names = FALSE)
  trt.ind <- which.pmin(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1, na.rm = TRUE)
  # optimal treatments
  test$A1opt <- ifelse(trt.ind == 1|trt.ind == 2, 1, -1)
  test$A2opt <- ifelse(trt.ind == 1|trt.ind == 3, 1, -1)
  return(test)
}

####################
#### simulation ####
####################
## Aim 1 (preliminary study, results in Web Appendix): impact of an unmeasured variable on stage 1 rule identification
set.seed(2021)
sigma <- 1

sim.uv <- function(n.train, n.test, generateTestData, generateTrainData, calcOut, gamma, n.sim){
  test <- generateTestData(n = n.test, calcOut = calcOut, C = 1.5, gamma = gamma)
  A1.perc <- data.frame(matrix(ncol = 4, nrow = n.sim))
  colnames(A1.perc) = c("sQ", "mQ", "IQ", "mIQ")
  for (i in 1:n.sim){
    train <- generateTrainData(n = n.train, calcOut = calcOut, C = 1.5, gamma = gamma)
    res.sQ <- sQ.pred(train, test)
    res.mQ <- mQ.pred(train, test)
    res.iQ <- iQ.pred(train, test)
    res.imQ <- imQ.pred(train, test)
    A1.perc[i, ] <- c(sum(res.sQ$A1opt == test$A1opt)/nrow(test), 
                      sum(res.mQ$A1opt == test$A1opt)/nrow(test), 
                      sum(res.iQ$A1opt == test$A1opt)/nrow(test), 
                      sum(res.imQ$A1opt == test$A1opt)/nrow(test))
  }
  res <- c(gamma, apply(A1.perc, 2, mean))
  return(res)
}
res.uncor <- data.frame()
res.cor <- data.frame()
res.heteroA1 <- data.frame()
gamma.grid <- seq(0, 4, 1)
for (i in 1:length(gamma.grid)){
  res.uncor <- rbind.data.frame(res.uncor, sim.uv(n.train = 250, n.test = 10000, 
                                                  generateTestData = generateTestData.uncorUV, 
                                                  generateTrainData = generateTrainData.uncorUV, 
                                                  calcOut = calcOut.homoA2, 
                                                  gamma = gamma.grid[i], n.sim = 100))
  res.cor <- rbind.data.frame(res.cor, sim.uv(n.train = 250, n.test = 10000, 
                                              generateTestData = generateTestData.corUV, 
                                              generateTrainData = generateTrainData.corUV, 
                                              calcOut = calcOut.homoA2, 
                                              gamma = gamma.grid[i], n.sim = 100))
  res.heteroA1 <- rbind.data.frame(res.heteroA1, sim.uv(n.train = 250, n.test = 10000, 
                                                        generateTestData = generateTestData.heteroA1, 
                                                        generateTrainData = generateTrainData.heteroA1, 
                                                        calcOut = calcOut.homoA2, 
                                                        gamma = gamma.grid[i], n.sim = 100))
  colnames(res.uncor) = colnames(res.cor) = colnames(res.heteroA1) = c("gamma", "sQ", "mQ", "IQ", "mIQ")
}

## Aim 2 (main manuscript): compare four methods for different sizes of heterogeneous treatment effects
# metric: percentage of correctly identified dhat1(X1) and dhat2(X1,A1opt,X2)
set.seed(2020)
sim.imq <- function(n.train, n.test, C, gamma, n.sim){
  test <- generateTestData.heteroA1(n = n.test, calcOut = calcOut.heteroA2, C = C, gamma = gamma)
  A1.perc <- data.frame(matrix(ncol = 4, nrow = n.sim))
  colnames(A1.perc) = c("sQ", "mQ", "IQ", "mIQ")
  for (i in 1:n.sim){
    train <- generateTrainData.heteroA1(n = n.train, calcOut = calcOut.heteroA2, C = C, gamma = gamma)
    res.sQ <- sQ.pred(train, test)
    res.mQ <- mQ.pred(train, test)
    res.iQ <- iQ.pred(train, test)
    res.imQ <- imQ.pred(train, test)
    A1.perc[i, ] <- c(sum(res.sQ$A1opt == test$A1opt)/nrow(test), 
                      sum(res.mQ$A1opt == test$A1opt)/nrow(test), 
                      sum(res.iQ$A1opt == test$A1opt)/nrow(test), 
                      sum(res.imQ$A1opt == test$A1opt)/nrow(test))
  }
  res <- c(C, gamma, apply(A1.perc, 2, mean))
  return(res)
}
res <- data.frame()
C.grid <- seq(1, 3, 0.5) # for table
gamma.grid <- c(0, 2, 4)
for (i in 1:length(gamma.grid)){
  for (j in 1:length(C.grid)){
    res <- rbind.data.frame(res, sim.imq(n.train = 250, n.test = 10000, C = C.grid[j], gamma = gamma.grid[i], n.sim = 100))
    colnames(res) <- c("C", "gamma", "sQ", "mQ", "IQ", "mIQ")
  }
}
write.csv(res, "pci_table.csv", row.names = FALSE)
