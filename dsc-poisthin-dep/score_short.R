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
    FSP_005 = sum(data$meta$beta*output$beta.est<0 & output$svalue<=0.05,na.rm=TRUE)/sum(output$svalue<=0.05,na.rm=TRUE)
    
    AUC = auc(roc(output$qvalue,factor(data$meta$null)))
    rmse.beta = sqrt(mean((data$meta$beta-output$beta.est)^2))
    
    res = c(pi0, pi0.est, DP_005, FDP_005, FSP_005, AUC, rmse.beta)
    names(res) = c("pi0","pi0.est","DP_005","FDP_005","FSP_005","AUC","rmse.beta") 
    return(res)
  }else{
    res = rep(NA, 7)
    names(res) = c("pi0","pi0.est","DP_005","FDP_005","FSP_005","AUC","rmse.beta")  
    return(res)
  } 
}


score3 = function(data, output){
  return(c(S=pcdf_post(output$fit$fitted.g, data$meta$beta,
                       output$fit$data$betahat,output$fit$data$sebetahat,v=output$fit$data$df)))
}

score_neg = function(data, output){
  return(c(S=output$fit$NegativeProb))
}

score_pos = function(data, output){
  return(c(S=output$fit$PositiveProb))
}