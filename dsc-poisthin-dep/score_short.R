#This file should define your score function
library(AUC)

score = function(data, output){
  if (class(output)=="list"){
    if (is.null(output$svalue)){
      output$svalue = NA
    }
    
    pi0 = data$meta$pi0
    pi0.est = output$pi0
    
    DP_005 = mean(output$qvalue<=0.05,na.rm=TRUE)
    FDP_005 = sum(data$meta$null==1 & output$qvalue<=0.05,na.rm=TRUE)/sum(output$qvalue<=0.05,na.rm=TRUE)
    FSP_005 = sum(data$meta$beta*output$beta.est<=0 & output$svalue<=0.05,na.rm=TRUE)/sum(output$svalue<=0.05,na.rm=TRUE)
    
    AUC = auc(roc(output$qvalue,factor(data$meta$null)))
    rmse.beta = sqrt(mean((data$meta$beta-output$beta.est)^2))
    
    res = c(pi0, pi0.est, DP_005, FDP_005, FSP_005, AUC, rmse.beta)
  }else{
    res = rep(NA, 7)
  } 
  names(res) = c("pi0","pi0.est","DP_005","FDP_005","FSP_005","AUC","rmse.beta") 
  return(res)
}

score3 = function(data, output){
  if (class(output)=="list"){
    return(c(S=pcdf_post(output$fit$fitted_g, data$meta$beta, output$fit$data)))
  }else{
    return(c(S=rep(NA,length(data$meta$beta))))
  }
  
}

score_neg = function(data, output){
  if (class(output)=="list"){
    return(c(S=output$fit$result$NegativeProb))
  }else{
    return(c(S=rep(NA,length(data$meta$beta)))) 
  }
}

score_pos = function(data, output){
  if (class(output)=="list"){
    return(c(S=output$fit$result$PositiveProb))
  }else{
    return(c(S=rep(NA,length(data$meta$beta))))
  }
}

score_loglike = function(data, output){
  if (class(output)=="list"){
    return(output$fit$loglik)
  }else{
    return(NA)
  }
}

score_lambda = function(data, output){
  if (class(output)=="list"){
    res = c(output$lambda1, output$lambda2)
  }else{
    res = c(NA, NA)
  }
  names(res) = c("lambda1", "lambda2")
  return(res)
}

