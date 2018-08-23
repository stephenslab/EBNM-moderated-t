nullash.simp.wrapper = function(input, args){
  # variance shrinkage
  ash.params <- choose.ash.parameters(input, args)
  va = vash(ash.params$sebetahat/ash.params$scale,
            df=ash.params$df,singlecomp=TRUE)
  sebetahat.new = ash.params$scale*c(sqrt(va$PosteriorRate/va$PosteriorShape))
  df.new = 2*va$PosteriorShape[1]
  betahat = ash.params$betahat
  
  # use the first nctl genes as control genes
  ctls = rep(0,length(input$null))  # Use nctl true nulls to do supervised RUV/SVA
  ctls[which(input$null==1)[1:args$nctl]] = 1
  
  # use controls to fit lambda, such that t-score|H0 ~ lambda*T(df)
  # MoM estimator: match 2nd moment
  lambda = sqrt(sum(((betahat/sebetahat.new)[ctls==1])^2)/
                  ((sum(ctls==1)-1)*df.new/(df.new-2)))
  
  fit = ash(betahat, lambda*sebetahat.new, df=df.new, mixcompdist = "uniform", alpha=args$alpha)
  
  return(list(fit=fit, lambda=lambda))
}