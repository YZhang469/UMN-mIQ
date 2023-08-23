###############################
## Functions to Analyze Data ##
###############################

# for non-interactive Q-learning
predTrain <- function(dat, trt.name, mod, modified = FALSE){
  dat.opt <- dat
  dat.n1 <- dat.opt
  dat.n1[, trt.name] <- -1
  dat.n1$Y <- predict(mod, newdata = dat.n1)
  dat.1 <- dat.opt
  dat.1[, trt.name] <- 1
  dat.1$Y <- predict(mod, newdata = dat.1)
  dat.opt[, trt.name] <- ifelse(dat.n1$Y < dat.1$Y, -1, 1)
  if (!modified){
    dat.opt$Y <- predict(mod, newdata = dat.opt)
  }
  if (modified){
    int <- (dat.1$Y - dat.n1$Y)/2
    dat.opt$Y <- dat$Y - 2 * ifelse(dat.opt[, trt.name] == dat[, trt.name], 0, 1) * abs(int)
  }
  return(dat.opt)
}
predTest <- function(dat, trt.name, mod){
  dat.opt <- dat
  if (trt.name == "A2"){
    dat.opt$A1 <- dat.opt$A1opt # base estimation of A2 on true A1opt
  }
  dat.n1 <- dat.opt
  dat.n1[, trt.name] <- -1
  dat.n1$Y <- predict(mod, newdata = dat.n1)
  dat.1 <- dat.opt
  dat.1[, trt.name] <- 1
  dat.1$Y <- predict(mod, newdata = dat.1)
  dat.opt[, trt.name] <- ifelse(dat.n1$Y < dat.1$Y, -1, 1)
  if (trt.name == "A2"){
    dat.opt$main2hat <- (dat.n1$Y + dat.1$Y)/2
    dat.opt$trt2hat <- (dat.1$Y - dat.n1$Y)/2
  }
  if (trt.name == "A1"){
    dat.opt$HTE1hat <- dat.1$Y - dat.n1$Y
    dat.opt$Yopthat <- pmin(dat.1$Y, dat.n1$Y)
  }
  return(dat.opt)
}

# for interactive Q-learning
# use the estimated main model to predict the output "fitIQ1main" for test data
predIQ1main <- function(fitIQ1main, H1Main, A1, s1mainInts){
  s1mainVars <- as.data.frame(H1Main)
  alphaHat <- c(fitIQ1main$alphaHat0, fitIQ1main$alphaHat1)
  mainPos <- as.matrix(cbind (1, s1mainVars, 1, s1mainVars[, s1mainInts])) %*% alphaHat
  mainNeg <- as.matrix(cbind (1, s1mainVars, -1, -s1mainVars[, s1mainInts])) %*% alphaHat
  list ("alphaHat0" = fitIQ1main$alphaHat0, "alphaHat1" = fitIQ1main$alphaHat1, 
        "s1MainFit" = fitIQ1main$s1MainFit, "mainPos" = mainPos, "mainNeg" = mainNeg, 
        "s1mainInts" = s1mainInts, "A1" = A1)
}
predIQ1cm <- function(fitIQ1cm, H1CMean, A1, s1cmInts){
  s1cmVars = as.data.frame (H1CMean)
  betaHat1 <- c(fitIQ1cm$betaHat10, fitIQ1cm$betaHat11)
  cmPos = as.matrix (cbind(1, s1cmVars, 1, s1cmVars[, s1cmInts])) %*% betaHat1
  cmNeg = as.matrix (cbind(1, s1cmVars, -1, -s1cmVars[, s1cmInts])) %*% betaHat1
  list ("betaHat0" = fitIQ1cm$betaHat10, "betaHat1" = fitIQ1cm$betaHat11, 
        "s1cmFit" = fitIQ1cm$s1cmFit, 
        "cmPos" = cmPos, "cmNeg" = cmNeg, 
        "s1cmInts" = s1cmInts, "A1" = A1)
}
learnIQ1Est <- function (mainObj, cmObj, sigObj, dens){
  A1 = cmObj$A1;
  
  if (dens=="norm"){
    normal = T;
  }
  else if (dens=="nonpar"){
    normal = F;
  }
  else{
    stop ("dens must be one of {norm, nonpar}");
  }
  ## NEXT: estimate optimal treatment for each patient in training
  ## set
  
  ## contrast function mean estimate
  muPos = cmObj$cmPos;
  muNeg = cmObj$cmNeg;
  
  ## lhat is the estimate of the expectation of the main effect term
  lhatPos = mainObj$mainPos; 
  lhatNeg = mainObj$mainNeg;
  n = length (lhatPos);
  ## standard deviation estimates
  if (sigObj$homo){
    sigPos = rep (sigObj$sigPos, n);
    sigNeg = rep (sigObj$sigNeg, n);
  } else {
    sigPos = sigObj$sigPos;
    sigNeg = sigObj$sigNeg;
  }
  ## Estimate Q1 for both A1=1 and A1=-1 using either normal density
  ## or empirical estimate
  if (normal){
    q1HatPos = lhatPos - muPos*(1-2*pnorm (-muPos/sigPos)) - sqrt(2/pi)*sigPos*exp (-muPos^2/(2*sigPos^2)); 
    q1HatNeg = lhatNeg - muNeg*(1-2*pnorm (-muNeg/sigNeg)) - sqrt(2/pi)*sigNeg*exp (-muNeg^2/(2*sigNeg^2)); 
  }
  else{
    n = length (lhatPos)
    q1HatPos = rep (NA, n);
    q1HatNeg = rep (NA, n);
    for (i in 1:n){
      opPos = sum (abs (muPos[i] + sigPos[i]*sigObj$stdResids));  
      opNeg = sum (abs (muNeg[i] + sigNeg[i]*sigObj$stdResids));  
      q1HatPos[i] = lhatPos[i] - opPos;
      q1HatNeg[i] = lhatNeg[i] - opNeg;
    }
  }
  
  ## Vector of optimal first-stage txts
  optA1 = ifelse(q1HatPos < q1HatNeg, 1, -1);
  
  list ("optA1" = optA1, "HTE1hat" = q1HatPos - q1HatNeg, "Yopthat" = pmin(q1HatPos, q1HatNeg))
}

##############################
##### Modified Functions #####
##############################

# for interactive Q-learning
iqQ1MainEst <-
  function (mainResp, H1Main, A1, s1mainInts, ...){
    # mainResp = Object$main;
    s1mainVars = as.data.frame(H1Main);
    s1names = names(s1mainVars);
    s1mainInts = as.vector(s1mainInts);
    
    if (length (s1mainInts) > 0){  
      s1m. = as.matrix(cbind(1, s1mainVars, A1,
                             A1*s1mainVars[, s1mainInts])); 
      colnames(s1m.) = c ("intercept", s1names, "A1",
                          paste(s1names[s1mainInts], "A1", sep = ":")); 
      p10 = ncol(s1mainVars);
      p11 = ncol(s1mainVars[,s1mainInts]);
      p1 = ncol(s1m.);
      H10. = s1m.[, 1:(p10+1)];
      H11. = s1m.[, (p10+2):p1]*A1;
      
      ## first-stage regression for the main effect term
      s1MainFit = lm (mainResp ~ s1m. - 1, ...);
      alphaHat = s1MainFit$coefficients;
      alphaHat0 = alphaHat[1:(p10+1)];
      alphaHat1 = alphaHat[(p10+2):p1];
      
      mainPos = as.matrix(cbind(1, s1mainVars, 1,
                                s1mainVars[,s1mainInts])) %*% alphaHat;
      mainNeg = as.matrix(cbind(1, s1mainVars, -1,
                                -s1mainVars[,s1mainInts])) %*% alphaHat;
    }
    else{
      s1m. = as.matrix(cbind(1, s1mainVars, A1)); 
      colnames(s1m.) = c("intercept", s1names, "A1"); 
      p10 = ncol(s1mainVars);
      p1 = ncol(s1m.);
      H10. = s1m.[, 1:(p10+1)];
      H11. = 1;
      
      ## first-stage regression for the main effect term
      s1MainFit = lm(mainResp ~ s1m. - 1, ...);
      alphaHat = s1MainFit$coefficients;
      alphaHat0 = alphaHat[1:(p10+1)];
      alphaHat1 = alphaHat[p1];
      
      mainPos = as.matrix (cbind (1, s1mainVars, 1)) %*% alphaHat;
      mainNeg = as.matrix (cbind (1, s1mainVars, -1)) %*% alphaHat;
    }
    
    list ("alphaHat0" = alphaHat0, "alphaHat1" = alphaHat1,
          "s1MainFit" = s1MainFit, "mainPos" = mainPos, "mainNeg" = mainNeg,
          "s1mainInts" = s1mainInts, "A1" = A1);        
  }
