# Just use the log-cpm, but don't use the voom weights

library(ashr)
library(limma)
noweightash.wrapper = function(input,args){
  
  lim = lmFit(input$v$E, design=model.matrix(~input$condition))
  lim = eBayes(lim)
  fit = ash(lim$coefficients[,2], sqrt(lim$s2.post)*lim$stdev.unscaled[,2],
            df = lim$df.total[1], mixcompdist = "uniform")
  
  return(list(fit=fit))
}