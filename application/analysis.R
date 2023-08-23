library(dplyr)
library(iqLearn)
source("HelpFunctions.R")
load("widedata.Rdata")
apply(mbridge, 2, function(x) sum(is.na(x)))

covnames <- c("nonwhite", "female", "yes_greek", 
              "norm_percent_drink", "norm_num_drinks_typical", "norm_most_drinks_typical", "norm_percent_binge", 
              "intent_freq_per_month", "intent_typical_drinks", "intent_drunk_per_month", 
              "yes_parent", "num_days0", "avg_drinks0")

analyze <- function(Yname = "max_drinks", covnames = covnames){ # Yname can start with "binge", "max_drinks", "byaacq"
  dat <- mbridge %>% select(External_Data_Reference, nonwhite, female, yes_greek, 
                            starts_with("norm_"), starts_with("intent_"), yes_parent, num_days0, avg_drinks0, 
                            A1, A2, HD, 
                            starts_with(Yname))
  dat <- dat[, -ncol(dat)] # remove follow-up 2
  dat <- dat[apply(dat, 1, function(x) sum(is.na(x))==0), ]
  
  Y <- paste0(Yname, 1, sep = "")
  Y0 <- paste0(Yname, 0, sep = "")
  covnames <- c(covnames, Y0)
  f2 <- as.formula(paste0(Y, " ~ A1 + ", paste(covnames, collapse = " + "), 
                          " + A2 * (A1 + ", paste(covnames, collapse = " + "), ")", sep = ""))
  f2.glm <- as.formula(paste0(Y, " ~ A2 * (A1 + ", paste(covnames, collapse = " + "), ")", sep = ""))
  f1.glm <- as.formula(paste0("Y ~ A1 * (", paste(covnames, collapse = " + "), ")", sep = ""))
  
  # modified interactive Q-learning - normal
  dat.s2 <- dat[dat$HD == 1, ]
  fitIQ2 <- learnIQ2(f2, data = dat.s2, treatName = "A2", intNames = c("A1", covnames))
  HTE2 <- 2 * fitIQ2$contrast
  dat.s1 <- dat
  dat.s1$contrast <- 0
  dat.s1$contrast[dat.s1$HD != 0] <- fitIQ2$contrast
  s1vars <- dat.s1[, c(2:14,18)]
  s1mainInts <- 1:14
  fitIQ1main <- iqQ1MainEst(mainResp = dat.s1[, Y] - dat.s1$A2*dat.s1$contrast, 
                            H1Main = s1vars, A1 = dat.s1$A1, s1mainInts = s1mainInts)
  fitIQ1cm <- learnIQ1cm(object = fitIQ2, H1CMean = dat.s2[, c(2:14,18)], A1 = dat.s2$A1, s1cmInts = 1:14)
  fitIQ1var <- learnIQ1var(object = fitIQ1cm, method = "homo")
  fitIQ1 <- learnIQ1Est(tailor = dat.s1$HD, mainObj = fitIQ1main, cmObj = fitIQ1cm, sigObj = fitIQ1var, dens = "norm")
  HTE1 <- fitIQ1$HTE
  d1opt <- fitIQ1$optA1
  H21 <- cbind.data.frame("intercept" = rep(1, nrow(dat.s2)), "A1" = d1opt[dat.s1$HD == 1], 
                          dat.s2[, covnames])
  d2opt <- ifelse(as.matrix(H21) %*% fitIQ2$betaHat21 > 0, -1, 1)
  dopt.mIQ <- list("d1opt" = d1opt, "d2opt" = d2opt)
  HTE.mIQ <- cbind.data.frame("Stage" = c(rep("Stage 1", length(HTE1)), rep("Stage 2", length(HTE2))), 
                              "HTE" = c(HTE1, HTE2))
  
  # modified interactive Q-learning - nonparametric
  fitIQ1 <- learnIQ1Est(tailor = dat.s1$HD, mainObj = fitIQ1main, cmObj = fitIQ1cm, sigObj = fitIQ1var, dens = "nonpar")
  HTE1 <- fitIQ1$HTE
  d1opt <- fitIQ1$optA1
  H21 <- cbind.data.frame("intercept" = rep(1, nrow(dat.s2)), "A1" = d1opt[dat.s1$HD == 1], 
                          dat.s2[, covnames])
  d2opt <- ifelse(as.matrix(H21) %*% fitIQ2$betaHat21 > 0, -1, 1)
  dopt.mIQnonpar <- list("d1opt" = d1opt, "d2opt" = d2opt)
  
  # interactive Q-learning
  dat.s2 <- dat[dat$HD == 1, ]
  fitIQ2 <- learnIQ2(f2, data = dat.s2, treatName = "A2", intNames = c("A1", covnames))
  dat.s1 <- dat
  dat.s1$main <- 0
  dat.s1$main[dat.s1$HD != 0] <- fitIQ2$main
  dat.s1$main[dat.s1$HD == 0] <- dat.s1[dat.s1$HD == 0, Y]
  s1vars <- dat.s1[, c(2:14,18)]
  s1mainInts <- 1:14
  fitIQ1main <- iqQ1MainEst(mainResp = dat.s1$main, H1Main = s1vars, A1 = dat.s1$A1, s1mainInts = s1mainInts)
  fitIQ1cm <- learnIQ1cm(object = fitIQ2, H1CMean = dat.s2[, c(2:14,18)], A1 = dat.s2$A1, s1cmInts = 1:14)
  fitIQ1var <- learnIQ1var(object = fitIQ1cm, method = "homo")
  fitIQ1 <- learnIQ1Est(tailor = dat.s1$HD, mainObj = fitIQ1main, cmObj = fitIQ1cm, sigObj = fitIQ1var, dens = "norm")
  d1opt <- fitIQ1$optA1
  H21 <- cbind.data.frame("intercept" = rep(1, nrow(dat.s2)), "A1" = d1opt[dat.s1$HD == 1], 
                          dat.s2[, covnames])
  d2opt <- ifelse(as.matrix(H21) %*% fitIQ2$betaHat21 > 0, -1, 1)
  dopt.IQ <- list("d1opt" = d1opt, "d2opt" = d2opt)
  
  # interactive Q-learning - nonparametric
  fitIQ1 <- learnIQ1Est(tailor = dat.s1$HD, mainObj = fitIQ1main, cmObj = fitIQ1cm, sigObj = fitIQ1var, dens = "nonpar")
  d1opt <- fitIQ1$optA1
  H21 <- cbind.data.frame("intercept" = rep(1, nrow(dat.s2)), "A1" = d1opt[dat.s1$HD == 1], 
                          dat.s2[, covnames])
  d2opt <- ifelse(as.matrix(H21) %*% fitIQ2$betaHat21 > 0, -1, 1)
  dopt.IQnonpar <- list("d1opt" = d1opt, "d2opt" = d2opt)
  
  # standard Q-learning
  dat.s2 <- dat[dat$HD == 1, ]
  mod.s2 <- glm(f2.glm, data = dat.s2)
  dat.opt <- predDat(dat = dat.s2, trt.name = "A2", mod = mod.s2, modified = FALSE)
  dat.s1 <- dat
  dat.s1$Y <- dat.s1[, Y]
  dat.s1$Y[dat.s1$HD != 0] <- dat.opt$Y
  mod.s1 <- glm(f1.glm, data = dat.s1)
  dat.opt <- predDat(dat = dat.s1, trt.name = "A1", mod = mod.s1, modified = FALSE)
  d1opt <- dat.opt$A1
  H21 <- cbind.data.frame("intercept" = rep(1, nrow(dat.s2)), "A1" = dat.opt$A1[dat.opt$HD == 1], 
                          dat.s2[, covnames])
  d2opt <- ifelse(as.matrix(H21) %*% coefficients(mod.s2)[c(2,18:32)] > 0, -1, 1)
  dopt.Q <- list("d1opt" = d1opt, "d2opt" = d2opt)
  
  # modified Q-learning
  dat.s2 <- dat[dat$HD == 1, ]
  mod.s2 <- glm(f2.glm, data = dat.s2)
  dat.opt <- predDat(dat = dat.s2, trt.name = "A2", mod = mod.s2, modified = TRUE, Y = dat.s2[, Y])
  dat.s1 <- dat
  dat.s1$Y <- dat.s1[, Y]
  dat.s1$Y[dat.s1$HD != 0] <- dat.opt$Y
  mod.s1 <- glm(f1.glm, data = dat.s1)
  dat.opt <- predDat(dat = dat.s1, trt.name = "A1", mod = mod.s1, modified = TRUE, Y = dat.s1$Y)
  d1opt <- dat.opt$A1
  H21 <- cbind.data.frame("intercept" = rep(1, nrow(dat.s2)), "A1" = dat.opt$A1[dat.opt$HD == 1], 
                          dat.s2[, covnames])
  d2opt <- ifelse(as.matrix(H21) %*% coefficients(mod.s2)[c(2,18:32)] > 0, -1, 1)
  dopt.mQ <- list("d1opt" = d1opt, "d2opt" = d2opt)
  
  return(rbind("Q" = c(sapply(dopt.Q, table)), "mQ" = c(sapply(dopt.mQ, table)), 
               "IQ" = c(sapply(dopt.IQ, table)), "IQnonpar" = c(sapply(dopt.IQnonpar, table)), 
               "mIQ" = c(sapply(dopt.mIQ, table)), "mIQnonpar" = c(sapply(dopt.mIQnonpar, table))))
}

analyze(Yname = "max_drinks", covnames = covnames)
analyze(Yname = "byaacq", covnames = covnames)

# model diagnostics of the g function
diagnose <- function(Yname = "max_drinks", covnames = covnames){
  dat <- mbridge %>% select(External_Data_Reference, nonwhite, female, yes_greek, 
                            starts_with("norm_"), starts_with("intent_"), yes_parent, num_days0, avg_drinks0, 
                            A1, A2, HD, 
                            starts_with(Yname))
  dat <- dat[, -ncol(dat)] # remove follow-up 2
  dat <- dat[apply(dat, 1, function(x) sum(is.na(x))==0), ]
  Y <- paste0(Yname, 1, sep = "")
  Y0 <- paste0(Yname, 0, sep = "")
  covnames <- c(covnames, Y0)
  f2 <- as.formula(paste0(Y, " ~ A1 + ", paste(covnames, collapse = " + "), 
                          " + A2 * (A1 + ", paste(covnames, collapse = " + "), ")", sep = ""))
  dat.s2 <- dat[dat$HD == 1, ]
  fitIQ2 <- learnIQ2(f2, data = dat.s2, treatName = "A2", intNames = c("A1", covnames))
  HTE2 <- 2 * fitIQ2$contrast
  mod.HTE2 <- lm(as.formula(paste0("HTE2 ~ A1 * (", paste(covnames, collapse = " + "), ")", sep = "")), data = dat.s2)
  png(filename = paste0("diag_HTE2_", Yname, ".png", sep = ""), width = 18, height = 7, units = "cm", res = 600)
  par(mfrow = c(1, 3), oma = c(0, 0, 0, 0))
  plot(mod.HTE2, which = 1:2, sub.caption = "")
  hist(mod.HTE2$residuals, xlab = "Residuals", main = "Residual Histogram")
  dev.off()
  return(shapiro.test(mod.HTE2$residuals))
}
diagnose(Yname = "max_drinks", covnames = covnames)
diagnose(Yname = "byaacq", covnames = covnames)
