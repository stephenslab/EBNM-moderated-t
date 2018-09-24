# ad-hoc ash:
# use limma's betahat and p-value of moderated t-test
# to recover the "sebetahat" which makes the z-score
# produce the same p-value
# then feed betahat and "sebetahat" into normal ash

library(limma)
library(ashr)
adhocash.wrapper = function(input,args){
  lim = lmFit(input$v)
  lim = eBayes(lim)
  betahat = lim$coefficients[,2]
  pval = lim$p.value[,2]
  zscore = qnorm(1-pval/2)
  sebetahat = abs(betahat/zscore)

  fit = ash(betahat,sebetahat,method="fdr")

  return(list(fit=fit))
}