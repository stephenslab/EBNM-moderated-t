library(ashr)
nullash.wrapper = function(input, args){
  voom = voom_transform(input$counts, input$condition)
  betahat = voom$betahat
  sebetahat = voom$sebetahat

  # use the first nctl genes as control genes
  ctls = rep(0,length(input$null))  # Use nctl true nulls to do supervised RUV/SVA
  ctls[which(input$null==1)[1:args$nctl]] = 1
  
  # use controls to fit g0
  nullfit = ash(betahat[ctls==1], sebetahat[ctls==1], mixcompdist = "normal",prior="uniform")
  
  if (nullfit$fitted_g$pi[nullfit$fitted_g$sd==0]!=1){
    # if the fitted g0 is not point mass on 0
    pilik = nullfit$fitted_g$pi
    sdlik = sqrt(outer(sebetahat^2,(nullfit$fitted_g$sd)^2,FUN="+"))
    fit = ash(betahat, 1, 
              lik=normalmix_lik(pilik,sdlik),mixcompdist = "uniform")
  }else{
    fit = ash(betahat, sebetahat, mixcompdist = "uniform")
  }
  
  return(list(fit=fit, nullfit=nullfit))
}