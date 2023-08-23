library(ggplot2)

source("functions.R")

###############################################
########## Data Generative Functions ##########
###############################################

# function to control stage 2 heterogeneous treatment effects
calcEY <- function(X1, A1, X2, A2, c2, HTE2 = "homo"){ # HTE2 could be "homo" or "hetero"
  H2 <- cbind.data.frame("intercept" = rep(1, length(X1)), "X1" = X1, "A1" = A1, "X2" = X2)
  beta20 <- c(3, -1, 0.1, -0.1)
  if (HTE2 == "homo"){
    beta21 <- c2 * c(-6, 0, 5, 0)
  }
  else if (HTE2 == "hetero"){
    beta21 <- c2 * c(-6, -4, 5, -0.2)
  }
  EY <- as.matrix(H2)%*%beta20 + A2*(as.matrix(H2)%*%beta21)
  return(EY)
}

# function to control stage 1 heterogeneous treatment effects (and omitted variables)
generateTrainData <- function(n = 250, 
                              UV = "HTE1", c1 = 0, 
                              c2 = 1.5, HTE2 = "homo"){ # UV could be "uncorrelated", "correlated", or "HTE1"
  X1 <- rnorm(n, -2, 1)
  A1 <- 2 * rbinom(n, 1, 0.5) - 1
  A2 <- 2 * rbinom(n, 1, 0.5) - 1
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma) # sd(phi) is too small, 7-12 is appropriate
  if (UV == "uncorrelated"){
    # unmeasured variable uncorrelated with stage 2 main effects predictors
    V <- rnorm(n, -1, 1)
  }
  else if (UV == "correlated"){
    # unmeasured variable correlated with stage 2 main effects predictors
    V <- rnorm(n, 2*X1*X2, 1)
  }
  else if (UV == "HTE1"){
    V <- X1*A1
  }
  Y <- calcEY(X1, A1, X2, A2, c2, HTE2) + c1*V + phi
  train <- data.frame("ID" = 1:n, "X1" = X1, "A1" = A1, "X2" = X2, "A2" = A2, "Y" = Y)
  return(train)
}

generateTestData <- function(n = 10000, 
                             UV = "HTE1", c1 = 0, 
                             c2 = 1.5, HTE2 = "homo"){ # UV could be "uncorrelated", "correlated", or "HTE1"
  X1 <- rnorm(n, -2, 1)
  epsilon <- rnorm(n, 0, 2)
  X2 <- X1 + epsilon
  phi <- rnorm(n, 0, sigma)
  if (UV == "uncorrelated"){
    V <- rnorm(n, -1, 1)
    Y.1.1 <- calcEY(X1, A1 = rep(1, n), X2, A2 = rep(1, n), c2, HTE2) + c1*V + phi
    Y.1.n1 <- calcEY(X1, A1 = rep(1, n), X2, A2 = rep(-1, n), c2, HTE2) + c1*V + phi
    Y.n1.1 <- calcEY(X1, A1 = rep(-1, n), X2, A2 = rep(1, n), c2, HTE2) + c1*V + phi
    Y.n1.n1 <- calcEY(X1, A1 = rep(-1, n), X2, A2 = rep(-1, n), c2, HTE2) + c1*V + phi
  }
  else if (UV == "correlated"){
    V <- rnorm(n, 2*X1*X2, 1)
    Y.1.1 <- calcEY(X1, A1 = rep(1, n), X2, A2 = rep(1, n), c2, HTE2) + c1*V + phi
    Y.1.n1 <- calcEY(X1, A1 = rep(1, n), X2, A2 = rep(-1, n), c2, HTE2) + c1*V + phi
    Y.n1.1 <- calcEY(X1, A1 = rep(-1, n), X2, A2 = rep(1, n), c2, HTE2) + c1*V + phi
    Y.n1.n1 <- calcEY(X1, A1 = rep(-1, n), X2, A2 = rep(-1, n), c2, HTE2) + c1*V + phi
  }
  else if (UV == "HTE1"){
    Y.1.1 <- calcEY(X1, A1 = rep(1, n), X2, A2 = rep(1, n), c2, HTE2) + c1*X1 + phi
    Y.1.n1 <- calcEY(X1, A1 = rep(1, n), X2, A2 = rep(-1, n), c2, HTE2) + c1*X1 + phi
    Y.n1.1 <- calcEY(X1, A1 = rep(-1, n), X2, A2 = rep(1, n), c2, HTE2) - c1*X1 + phi
    Y.n1.n1 <- calcEY(X1, A1 = rep(-1, n), X2, A2 = rep(-1, n), c2, HTE2) - c1*X1 + phi
  }
  HTE1_1 <- Y.1.1 - Y.n1.1
  HTE1_n1 <- Y.1.n1 - Y.n1.n1
  trt2_1 <- (Y.1.1 - Y.1.n1)/2
  trt2_n1 <- (Y.n1.1 - Y.n1.n1)/2
  main2_1 <- (Y.1.1 + Y.1.n1)/2 - phi
  main2_n1 <- (Y.n1.1 + Y.n1.n1)/2 - phi
  test <- data.frame("ID" = 1:n, "X1" = X1, "X2" = X2, 
                     "Y(1,1)" = Y.1.1, "Y(1,-1)" = Y.1.n1, 
                     "Y(-1,1)" = Y.n1.1, "Y(-1,-1)" = Y.n1.n1, 
                     "Yopt" = pmin(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1), 
                     check.names = FALSE)
  dtr.opt <- apply(cbind.data.frame(Y.1.1, Y.1.n1, Y.n1.1, Y.n1.n1), 1, which.min)
  # optimal treatments
  test$A1opt <- ifelse(dtr.opt == 1|dtr.opt == 2, 1, -1)
  test$A2opt <- ifelse(dtr.opt == 1|dtr.opt == 3, 1, -1)
  test$HTE1 <- ifelse(test$A2opt == 1, HTE1_1, HTE1_n1)
  test$trt2 <- ifelse(test$A1opt == 1, trt2_1, trt2_n1)
  test$main2 <- ifelse(test$A1opt == 1, main2_1, main2_n1)
  return(test)
}

############################################
########## Preliminary Simulation ##########
############################################

# Impact of an unmeasured variable on stage 1 rule identification

sigma <- 1

sim.uv <- function(n.sim, n.train, n.test, UV, c1){
  test <- generateTestData(n = n.test, UV = UV, c1 = c1, c2 = 1.5, HTE2 = "homo")
  A1.perc <- data.frame(matrix(ncol = 4, nrow = n.sim))
  A2.perc <- data.frame(matrix(ncol = 4, nrow = n.sim))
  trt2.bias <- rep(NA, n.sim)
  main2.bias <- rep(NA, n.sim)
  colnames(A1.perc) = c("sQ", "mQ", "IQ", "mIQ")
  for (i in 1:n.sim){
    train <- generateTrainData(n = n.train, UV = UV, c1 = c1, c2 = 1.5, HTE2 = "homo")
    res.sQ <- sQ.pred(train, test)
    res.mQ <- mQ.pred(train, test)
    res.IQ <- IQ.pred(train, test)
    res.mIQ <- mIQ.pred(train, test)
    A1.perc[i, ] <- c(sum(res.sQ$d1opt == test$A1opt)/n.test, 
                      sum(res.mQ$d1opt == test$A1opt)/n.test, 
                      sum(res.IQ$d1opt == test$A1opt)/n.test, 
                      sum(res.mIQ$d1opt == test$A1opt)/n.test)
    A2.perc[i, ] <- c(sum(res.sQ$d2opt == test$A2opt)/n.test, 
                      sum(res.mQ$d2opt == test$A2opt)/n.test, 
                      sum(res.IQ$d2opt == test$A2opt)/n.test, 
                      sum(res.mIQ$d2opt == test$A2opt)/n.test)
    trt2.bias[i] <- mean(res.sQ$trt2hat - test$trt2)
    main2.bias[i] <- mean(res.sQ$main2hat - test$main2)
  }
  res.perc <- c(c1, apply(A1.perc, 2, mean), apply(A2.perc, 2, mean))
  res.bias <- c(c1, sumMetric(trt2.bias), sumMetric(main2.bias))
  
  return(res = list("res.perc" = res.perc, "res.bias" = res.bias))
}

set.seed(2021)
res.perc.uncor <- data.frame()
res.perc.cor <- data.frame()
res.perc.HTE1 <- data.frame()
res.bias.uncor <- data.frame()
res.bias.cor <- data.frame()
res.bias.HTE1 <- data.frame()
c1.grid <- seq(0, 4, 1)
for (i in 1:length(c1.grid)){
  res.uncor <- sim.uv(n.sim = 100, n.train = 250, n.test = 10000, UV = "uncorrelated", c1 = c1.grid[i])
  res.cor <- sim.uv(n.sim = 100, n.train = 250, n.test = 10000, UV = "correlated", c1 = c1.grid[i])
  res.HTE1 <- sim.uv(n.sim = 100, n.train = 250, n.test = 10000, UV = "HTE1", c1 = c1.grid[i])
  res.perc.uncor <- rbind.data.frame(res.perc.uncor, res.uncor$res.perc)
  res.perc.cor <- rbind.data.frame(res.perc.cor, res.cor$res.perc)
  res.perc.HTE1 <- rbind.data.frame(res.perc.HTE1, res.HTE1$res.perc)
  res.bias.uncor <- rbind.data.frame(res.bias.uncor, res.uncor$res.bias)
  res.bias.cor <- rbind.data.frame(res.bias.cor, res.cor$res.bias)
  res.bias.HTE1 <- rbind.data.frame(res.bias.HTE1, res.HTE1$res.bias)
}
colnames(res.perc.uncor) = colnames(res.perc.cor) = colnames(res.perc.HTE1) = c("c1", "sQ1", "mQ1", "IQ1", "mIQ1", "sQ2", "mQ2", "IQ2", "mIQ2")
colnames(res.bias.uncor) = colnames(res.bias.cor) = colnames(res.bias.HTE1) = c("c1", "trt2", "main2")
res.perc1 <- rbind(round(res.perc.uncor$sQ1, digits = 3), round(res.perc.cor$sQ1, digits = 3), round(res.perc.HTE1$sQ1, digits = 3), 
                   round(res.perc.uncor$IQ1, digits = 3), round(res.perc.cor$IQ1, digits = 3), round(res.perc.HTE1$IQ1, digits = 3))
res.bias2 <- rbind(res.bias.uncor$main2, res.bias.cor$main2, res.bias.HTE1$main2, 
                   res.bias.uncor$trt2, res.bias.cor$trt2, res.bias.HTE1$trt2)

#####################################
########## Main Simulation ##########
#####################################

sigma <- 1

sim <- function(n.sim, n.train, n.test, c1, c2){
  test <- generateTestData(n = n.test, UV = "HTE1", c1 = c1, c2 = c2, HTE2 = "hetero")
  # Metric 1: percentage of correctly identified $\hat{d}_1^{\text{opt}}(X_1)$ and $\hat{d}_2^{\text{opt}}(X_1,A_1^{\text{opt}},X2)$
  A1.perc <- data.frame(matrix(ncol = 4, nrow = n.sim))
  # Metric 2: bias of the optimal value
  Yopt.bias <- data.frame(matrix(ncol = 4, nrow = n.sim))
  colnames(A1.perc) = c("sQ", "mQ", "IQ", "mIQ")
  for (iter in 1:n.sim){
    train <- generateTrainData(n = n.train, UV = "HTE1", c1 = c1, c2 = c2, HTE2 = "hetero")
    res.sQ <- sQ.pred(train, test)
    res.mQ <- mQ.pred(train, test)
    res.IQ <- IQ.pred(train, test)
    res.mIQ <- mIQ.pred(train, test)
    A1.perc[iter, ] <- c(sum(res.sQ$d1opt == test$A1opt)/n.test, 
                         sum(res.mQ$d1opt == test$A1opt)/n.test, 
                         sum(res.IQ$d1opt == test$A1opt)/n.test, 
                         sum(res.mIQ$d1opt == test$A1opt)/n.test)
    Yopt.bias[iter, ] <- c(mean(res.sQ$Yopthat-test$Yopt), 
                           mean(res.mQ$Yopthat-test$Yopt), 
                           mean(res.IQ$Yopthat-test$Yopt), 
                           mean(res.mIQ$Yopthat-test$Yopt))
  }
  return(list("res.A1" = c(c1, c2, apply(A1.perc, 2, mean)), 
              "res.Yopt" = c(c1, c2, apply(Yopt.bias, 2, sumMetric))))
}

set.seed(2023)
res.Yopt <- data.frame()
c1.grid <- c(0, 2, 4)
c2.grid <- seq(1, 3, 1) # for table
for (i in 1:length(c1.grid)){
  for (j in 1:length(c2.grid)){
    res.sim <- sim(n.sim = 100, n.train = 250, n.test = 10000, c1 = c1.grid[i], c2 = c2.grid[j])
    res.Yopt <- rbind.data.frame(res.Yopt, res.sim$res.Yopt)
  }
}

res.A1 <- data.frame()
c1.grid <- c(0, 2, 4)
c2.grid <- seq(0, 3, 0.2) # for figure
for (i in 1:length(c1.grid)){
  for (j in 1:length(c2.grid)){
    res.sim <- sim(n.sim = 1000, n.train = 250, n.test = 10000, c1 = c1.grid[i], c2 = c2.grid[j])
    res.A1 <- rbind.data.frame(res.A1, res.sim$res.A1)
  }
}

colnames(res.A1) = colnames(res.Yopt) = c("c1", "c2", "sQ", "mQ", "IQ", "mIQ")

write.csv(res.Yopt, "Yopt.csv", row.names = FALSE)
write.csv(res.A1, "A1.csv", row.names = FALSE)

res.pci <- reshape(res.A1, idvar = c("c1", "c2"), varying = list(3:6), v.names = "pci", times = names(res)[3:6], timevar = "method", direction = "long")
res.pci$c1 <- factor(res.pci$c1, levels = c("0", "2", "4"),
                     labels = c(expression(c[1]==0), expression(c[1]==2), expression(c[1]==4)))
pci.plot <- ggplot(res.pci, aes(x = c2, y = pci, shape = method, color = method, size = method)) +
  geom_line(size = 0.8) + guides(linetype = FALSE) +
  geom_point() +
  labs(title = "", x = expression(c[2]), y = "PCI") +
  ylim(0.8, 1) + scale_x_continuous(breaks = seq(0, 3, 0.5)) +
  scale_shape_manual(values = c(15, 16, 17, 18), 
                     name = "Method",
                     breaks = c("sQ", "mQ", "IQ", "mIQ"),
                     labels = c("Standard Q-learning", "Modified Q-learning",
                                "Interactive Q-learning", "Modified Interactive Q-learning")) + 
  scale_color_manual(# values = c("darkgoldenrod", "cadetblue", "lightsalmon", "darkolivegreen4"), 
                     # values = c("#D7FFF1", "#8CD790", "#77AF9C", "#285943"), 
                     values = c("#9DC8C8", "#58C9B9", "#519D9E", "#285943"), 
                     name = "Method",
                     breaks = c("sQ", "mQ", "IQ", "mIQ"),
                     labels = c("Standard Q-learning", "Modified Q-learning",
                                "Interactive Q-learning", "Modified Interactive Q-learning")) + 
  scale_size_manual(values = c(2.5, 3, 2.5, 3.5), 
                    name = "Method",
                    breaks = c("sQ", "mQ", "IQ", "mIQ"),
                    labels = c("Standard Q-learning", "Modified Q-learning",
                               "Interactive Q-learning", "Modified Interactive Q-learning")) + 
  facet_wrap(~ c1, ncol = 3, labeller = label_parsed) + 
  theme(legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), 
        strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
png(filename = "pci.png", width = 27, height = 12, units = "cm", res = 600)
pci.plot
dev.off()
