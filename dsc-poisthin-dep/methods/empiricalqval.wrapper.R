empiricalqval.wrapper = function(input,args){
  fit = lmFit(input$v)
  fit = eBayes(fit)
  pvalue = fit$p.value[,2]
  
  # use the first nctl genes as control genes
  ctls = rep(0,length(input$null))  # Use nctl true nulls to do supervised RUV/SVA
  ctls[which(input$null==1)[1:args$nctl]] = 1
  o = order(pvalue)
  qvalue=rep(NA,length(pvalue))
  qvalue[o] = cumsum(ctls[o])/(1:sum(!is.na(pvalue)))
  return(list(beta.est=fit$coefficients[,2],
              qvalue = qvalue,
              pi0 = NA))
}