library(iqLearn)
source("helpFunctions.R")

# Q-learning
sQ.pred <- function(train, test){
  # stage 2 estimation
  mod.s2 <- glm(Y ~ X1 + A1 + X2 + A2 * (X1 + A1 + X2), data = train)
  train.opt <- predTrain(dat = train, trt.name = "A2", mod = mod.s2)
  test.opt <- predTest(dat = test, trt.name = "A2", mod = mod.s2)
  # stage 1 estimation
  mod.s1 <- glm(Y ~ X1*A1, data = train.opt)
  train.opt <- predTrain(dat = train.opt, trt.name = "A1", mod = mod.s1)
  test.opt <- predTest(dat = test.opt, trt.name = "A1", mod = mod.s1)
  
  return(list("d1opt" = test.opt$A1, "d2opt" = test.opt$A2, 
              "HTE1hat" = test.opt$HTE1hat, "Yopthat" = test.opt$Yopthat, 
              "trt2hat" = test.opt$trt2hat, "main2hat" = test.opt$main2hat))
}

# modified Q-learning
mQ.pred <- function(train, test){
  # stage 2 estimation
  mod.s2 <- glm(Y ~ X1 + A1 + X2 + A2 * (X1 + A1 + X2), data = train)
  train.opt <- predTrain(dat = train, trt.name = "A2", mod = mod.s2, modified = TRUE)
  test.opt <- predTest(dat = test, trt.name = "A2", mod = mod.s2)
  # stage 1 estimation
  mod.s1 <- glm(Y ~ X1*A1, data = train.opt)
  train.opt <- predTrain(dat = train.opt, trt.name = "A1", mod = mod.s1, modified = TRUE)
  test.opt <- predTest(dat = test.opt, trt.name = "A1", mod = mod.s1)
  
  return(list("d1opt" = test.opt$A1, "d2opt" = test.opt$A2, 
              "HTE1hat" = test.opt$HTE1hat, "Yopthat" = test.opt$Yopthat, 
              "trt2hat" = test.opt$trt2hat, "main2hat" = test.opt$main2hat))
}

# interactive Q-learning
iQ.pred <- function(train, test){
  ## stage 2 prediction
  # estimate stage 2 model using training data
  # fitIQ2 <- learnIQ2(Y ~ X1 + A1 + X2 + A2 * (X1 + A1 + X2), data = train, treatName = "A2", intNames = c("X1", "A1", "X2"))
  s2vars <- train[, c(2:4)]
  s2ints <- 1:3
  fitIQ2 <- learnIQ2(H2 = s2vars, Y = train$Y, A2 = train$A2, s2ints = s2ints)
  # use the estimated model to predict A2opt for test data
  H21.test = H20.test = cbind.data.frame("intercept" = rep(1, nrow(test)), "X1" = test$X1, "A1" = test$A1opt, "X2" = test$X2)
  test.opt <- test
  test.opt$A2 <- ifelse(as.matrix(H21.test) %*% fitIQ2$betaHat21 > 0, -1, 1)
  test.opt$trt2hat <- as.matrix(H21.test) %*% fitIQ2$betaHat21
  test.opt$main2hat <- as.matrix(H20.test) %*% fitIQ2$betaHat20
  
  # step 2: stage 1 estimation of main effect
  # fitIQ1main <- learnIQ1main.formula(~ X1 + A1 * (X1), data = dat, treatName = "A1", intNames = "X1", s2object = fitIQ2)
  s1vars <- train[, 2, drop = FALSE]
  s1mainInts <- c(1)
  fitIQ1main <- iqQ1MainEst(mainResp = fitIQ2$main, H1Main = s1vars, A1 = train$A1, s1mainInts = s1mainInts)
  testIQ1main <- predIQ1main(fitIQ1main = fitIQ1main, H1Main = test[, 2, drop = FALSE], A1 = test$A1, s1mainInts = s1mainInts)
  
  # step 3: density modeling of contrast function
  fitIQ1cm <- learnIQ1cm(object = fitIQ2, H1CMean = s1vars, A1 = train$A1, s1cmInts = s1mainInts)
  fitIQ1var <- learnIQ1var(object = fitIQ1cm, method = "homo")
  testIQ1cm <- predIQ1cm(fitIQ1cm = fitIQ1cm, H1CMean = test[, 2, drop = FALSE], A1 = test$A1, s1cmInts = s1mainInts)
  
  # step 4: stage 1 estimation
  fitIQ1 <- learnIQ1Est(mainObj = fitIQ1main, cmObj = fitIQ1cm, sigObj = fitIQ1var, dens = "norm")
  testIQ1 <- learnIQ1Est(mainObj = testIQ1main, cmObj = testIQ1cm, sigObj = fitIQ1var, dens = "norm")
  test.opt$A1 <- testIQ1$optA1
  test.opt$HTE1hat <- testIQ1$HTE1hat
  test.opt$Yopthat <- testIQ1$Yopthat
  
  return(list("d1opt" = test.opt$A1, "d2opt" = test.opt$A2, 
              "HTE1hat" = test.opt$HTE1hat, "Yopthat" = test.opt$Yopthat, 
              "trt2hat" = test.opt$trt2hat, "main2hat" = test.opt$main2hat))
}

# modified interactive Q-learning
imQ.pred <- function(train, test){
  # step 1: regress Y on H2 and A2 to obtain Q2(H2,A2)
  s2vars <- train[, c(2:4)]
  s2ints <- 1:3
  fitIQ2 <- learnIQ2(H2 = s2vars, Y = train$Y, A2 = train$A2, s2ints = s2ints)
  # use the estimated model to predict A2opt for test data
  H21.test = H20.test = cbind.data.frame("intercept" = rep(1, nrow(test)), "X1" = test$X1, "A1" = test$A1opt, "X2" = test$X2)
  test.opt <- test
  test.opt$A2 <- ifelse(as.matrix(H21.test) %*% fitIQ2$betaHat21 > 0, -1, 1)
  test.opt$trt2hat <- as.matrix(H21.test) %*% fitIQ2$betaHat21
  test.opt$main2hat <- as.matrix(H20.test) %*% fitIQ2$betaHat20
  
  # step 2: regress Y - H21betahat21 on H1 and A1 to obtain l(H1,A1)
  # fitIQ1main <- learnIQ1main.formula(~ X1 + A1*(X1),
  #                            data = dat, treatName = "A1", intNames = "X1", s2object = fitIQ2)
  s1vars <- train[, 2, drop = FALSE]
  s1mainInts <- c(1)
  fitIQ1main <- iqQ1MainEst(mainResp = train$Y - fitIQ2$A2*fitIQ2$contrast, H1Main = s1vars, A1 = train$A1, s1mainInts = s1mainInts)
  testIQ1main <- predIQ1main(fitIQ1main = fitIQ1main, H1Main = test[, 2, drop = FALSE], A1 = test$A1, s1mainInts = s1mainInts)
  
  # step 3: estimate g(H21betahat21|H1,A1)
  fitIQ1cm <- learnIQ1cm(object = fitIQ2, H1CMean = s1vars, A1 = train$A1, s1cmInts = s1mainInts)
  fitIQ1var <- learnIQ1var(object = fitIQ1cm, method = "homo")
  testIQ1cm <- predIQ1cm(fitIQ1cm = fitIQ1cm, H1CMean = test[, 2, drop = FALSE], A1 = test$A1, s1cmInts = s1mainInts)
  
  # step 4: stage 1 estimation
  fitIQ1 <- learnIQ1Est(mainObj = fitIQ1main, cmObj = fitIQ1cm, sigObj = fitIQ1var, dens = "norm")
  testIQ1 <- learnIQ1Est(mainObj = testIQ1main, cmObj = testIQ1cm, sigObj = fitIQ1var, dens = "norm")
  test.opt$A1 <- testIQ1$optA1
  test.opt$HTE1hat <- testIQ1$HTE1hat
  test.opt$Yopthat <- testIQ1$Yopthat
  
  return(list("d1opt" = test.opt$A1, "d2opt" = test.opt$A2, 
              "HTE1hat" = test.opt$HTE1hat, "Yopthat" = test.opt$Yopthat, 
              "trt2hat" = test.opt$trt2hat, "main2hat" = test.opt$main2hat))
}
