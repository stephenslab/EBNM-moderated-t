# Use voom+limma+ash to test differential expression
# input$counts is the G*N expression count matrix (G: # of genes, N: # of samples)
# input$condition is the length N condition indicator vector 
#   e.g. condition=factor(c(1,1,2,2)) means that the first two columns of input$counts
#   belong to group "1" and the last two columns of input$counts belong to group "2"

library(limma)
library(edgeR)
library(ashr)

limmaash.wrapper = function(input,args=NULL){
  # calculate normalization factors
  dgecounts = calcNormFactors(DGEList(counts=input$counts,group=input$condition))
  
  # voom transformation
  design = model.matrix(~input$condition)
  v = voom(dgecounts,design,plot=FALSE)
  
  # use limma to estimate effect sizes and (shrunk) standard errors
  lim = lmFit(v)
  lim = eBayes(lim)
  betahat = lim$coefficients[,2]
  sebetahat = lim$stdev.unscaled[,2]*sqrt(lim$s2.post) # EB shrunk s.e.
  
  # fit ash with the shrunk s.e. and moderated d.f
  fit = ash(betahat,sebetahat,df=lim$df.total[1])
  
  return(list(fit=fit))
}