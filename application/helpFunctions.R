########################################
########## Modified Functions ##########
########################################
# for interactive Q-learning
iqQ1MainEst <-
  function (mainResp, H1Main, A1, s1mainInts, ...){
    # mainResp = Object$main;
    s1mainVars = as.data.frame (H1Main);
    s1names = names (s1mainVars);
    s1mainInts = as.vector (s1mainInts);
    
    if (length (s1mainInts) > 0){  
      s1m. = as.matrix (cbind (1, s1mainVars, A1,
                               A1*s1mainVars[,s1mainInts])); 
      colnames (s1m.) = c ("intercept", s1names, "A1",
                           paste(s1names[s1mainInts], "A1", sep=":")); 
      p10 = ncol (s1mainVars);
      p11 = ncol (s1mainVars[,s1mainInts]);
      p1 = ncol (s1m.);
      H10. = s1m.[, 1:(p10+1)];
      H11. = s1m.[, (p10+2):p1]*A1;
      
      ## first-stage regression for the main effect term
      s1MainFit = lm (mainResp ~ s1m. - 1, ...);
      alphaHat = s1MainFit$coefficients;
      alphaHat0 = alphaHat[1:(p10+1)];
      alphaHat1 = alphaHat[(p10+2):p1];
      
      mainPos = as.matrix (cbind (1, s1mainVars, 1,
                                  s1mainVars[,s1mainInts])) %*% alphaHat;
      mainNeg = as.matrix (cbind (1, s1mainVars, -1,
                                  -s1mainVars[,s1mainInts])) %*% alphaHat;
    }
    else{
      s1m. = as.matrix (cbind (1, s1mainVars, A1)); 
      colnames (s1m.) = c ("intercept", s1names, "A1"); 
      p10 = ncol (s1mainVars);
      p1 = ncol (s1m.);
      H10. = s1m.[, 1:(p10+1)];
      H11. = 1;
      
      ## first-stage regression for the main effect term
      s1MainFit = lm (mainResp ~ s1m. - 1, ...);
      alphaHat = s1MainFit$coefficients;
      alphaHat0 = alphaHat[1:(p10+1)];
      alphaHat1 = alphaHat[p1];
      
      mainPos = as.matrix (cbind (1, s1mainVars, 1)) %*% alphaHat;
      mainNeg = as.matrix (cbind (1, s1mainVars, -1)) %*% alphaHat;
    }
    
    list ("alphaHat0"=alphaHat0, "alphaHat1"=alphaHat1,
          "s1MainFit"=s1MainFit, "mainPos"=mainPos, "mainNeg"=mainNeg,
          "s1mainInts"=s1mainInts, "A1"=A1);        
  }

# modify function for SMARTs with embedded tailoring
learnIQ1Est <- function (tailor, mainObj, cmObj, sigObj, dens){
  A1 = mainObj$A1;
  
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
  muPos = rep(0, length(A1));
  muPos[tailor == 1] = cmObj$cmPos;
  muNeg = rep(0, length(A1));
  muNeg[tailor == 1] = cmObj$cmNeg;
  
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
    q1HatPos = lhatPos + muPos*(1-2*pnorm (-muPos/sigPos)) +
      sqrt(2/pi)*sigPos*exp (-muPos^2/(2*sigPos^2)); 
    q1HatNeg = lhatNeg + muNeg*(1-2*pnorm (-muNeg/sigNeg)) +
      sqrt(2/pi)*sigNeg*exp (-muNeg^2/(2*sigNeg^2)); 
  }
  else{
    n = length (lhatPos)
    q1HatPos = rep (NA, n);
    q1HatNeg = rep (NA, n);
    for (i in 1:n){
      opPos = sum (abs (muPos[i] +
                          sigPos[i]*sigObj$stdResids));  
      opNeg = sum (abs (muNeg[i] +
                          sigNeg[i]*sigObj$stdResids));  
      q1HatPos[i] = lhatPos[i] + opPos;
      q1HatNeg[i] = lhatNeg[i] + opNeg;
    }
  }
  
  ## Vector of optimal first-stage txts
  optA1 = - sign (q1HatPos - q1HatNeg);
  
  list ("optA1" = optA1, "HTE" = q1HatPos - q1HatNeg)
}

predDat <- function(dat, trt.name, mod, modified = FALSE, Y){
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
    dat.opt$Y <- Y - 2 * ifelse(dat.opt[, trt.name] == dat[, trt.name], 0, 1) * abs(int)
  }
  return(dat.opt)
}
