library(vicar)
library(sva)
vruv4ash.wrapper = function(input, args){
  Y = log(t(input$counts+1))
  X = model.matrix(~input$condition)
  
  # use the first nctl genes as control genes
  ctl = as.logical(c(rep(1,args$nctl),rep(0,nrow(input$counts)-args$nctl)))
  
  # estimate k (# of confounders)
  k = num.sv(as.matrix(log(input$counts+1)), 
             mod = stats::model.matrix(~input$condition), method = "be")
  k = min(k, nrow(X)-ncol(X)-1)
  
  if (k>0){
    # run vicarius_ruv4, cov_of_interest are all covariates but the intercept
    fit = ash_ruv4(Y, X, ctl, k=k, cov_of_interest = (1:ncol(X))[-1], 
                       likelihood = "t",
                       limmashrink = TRUE, include_intercept = FALSE,
                       scale_var = TRUE)
    return(list(fit=fit))
    
  }else{
    return(NA)
  }
  
}