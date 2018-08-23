library(ashr)
library(gaussquad)
nullash.approx.wrapper = function(input, args){
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
  
  # use controls to fit g0
  nullfit = ash(betahat[ctls==1], sebetahat.new[ctls==1], df=df.new, 
                mixcompdist = "uniform", prior="uniform", alpha=args$alpha)
  
  # approximate unimix distn nullfit$fitted_g by a normal distn (with same sd)
  if (args$pointmass==TRUE){
    nullg.pi0 = sum(nullfit$fitted_g$pi[nullfit$fitted_g$b==0])
    normsd = sqrt(sum(nullfit$fitted_g$pi*(nullfit$fitted_g$b)^2/3)/(1-nullg.pi0))
  }else{
    normsd = sqrt(sum(nullfit$fitted_g$pi*(nullfit$fitted_g$b)^2/3))
  }
  
  if (nullfit$fitted_g$pi[nullfit$fitted_g$a==0]!=1){
    # if the fitted g0 is not point mass on 0
    nonzerog0 = 1
    
    app = approxt(df.new, 10)
    if (args$pointmass==TRUE){
      pilik = c(outer(app$weight, c(nullg.pi0,1-nullg.pi0)))
      sdlik = cbind(outer(sebetahat.new^(1-args$alpha), app$sigma),
                    sqrt(outer(sebetahat.new^(2*(1-args$alpha)), app$sigma^2)+normsd^2))
    }else{
      pilik = app$weight
      sdlik = sqrt(outer(sebetahat.new^(2*(1-args$alpha)), app$sigma^2)+normsd^2)
    }
    betahat.new = betahat/sebetahat.new^(args$alpha)
    fit = ash(betahat.new, 1, lik=normalmix_lik(pilik,sdlik), 
              mixcompdist = "uniform")
  }else{
    nonzerog0 = 0
    fit = ash(betahat, sebetahat.new, df=df.new, 
              mixcompdist = "uniform", alpha=args$alpha)
    
  }
  
  return(list(fit=fit, nullfit=nullfit, nonzerog0=nonzerog0))
}

# Approximate t-distribution (of df) by r-components normal mixture
approxt = function(df, r){
  alpha = df/2-1
  rules = glaguerre.quadrature.rules(r,alpha,normalized=TRUE)
  sigma = sqrt(df/(2*rules[[r]]$x))
  weight = rules[[r]]$w/sum(rules[[r]]$w)
  return(list(sigma=sigma, weight=weight))
}