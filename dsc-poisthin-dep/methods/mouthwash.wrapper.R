library(vicar)
library(sva)
mouthwash.wrapper = function(input, args){
  Y = log2(t(input$counts+1))
  X = model.matrix(~input$condition)
  colnames(X) = c("Intercept", "Treatment")
  
  # estimate k (# of confounders)
  k = num_sv = max(sva::num.sv(t(Y), mod = X, method = "be"), 1)
  #k = min(k, nrow(X)-ncol(X)-1)
  

  # run mouthwash, cov_of_interest are all covariates but the intercept
  fit = mouthwash(Y, X, k=k, scale_var=TRUE)
  return(list(fit=fit))

}