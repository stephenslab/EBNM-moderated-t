ctlinflate.wrapper = function(input, args){
  # variance shrinkage
  ash.params = choose.ash.parameters(input, args)
  va = vash(ash.params$sebetahat/ash.params$scale,
            df=ash.params$df,singlecomp=TRUE)
  sebetahat.new = ash.params$scale*c(sqrt(va$PosteriorRate/va$PosteriorShape))
  df.new = 2*va$PosteriorShape[1]
  betahat = ash.params$betahat
  
  # use nctl genes as control genes
  ctls = rep(0,length(input$null))  # Use nctl true nulls to do supervised RUV/SVA
  idx = sample(which(input$null==1), min(args$nctl, sum(input$null)))
  ctls[idx] = 1
  
  if(args$normapprox==TRUE){
    # convert t scores to z scores by quantile matching
    sebetahat.norm = betahat/qnorm(pt(betahat/sebetahat.new, df=df.new))
    
    # suppose betahat~N(0,lambda1*sebetahat^2+lambda2)
    # use control genes to estimate lambda1, lambda2 (matching 2nd moment)
    betahat.ctl = betahat[ctls==1]
    sebetahat.ctl = sebetahat.norm[ctls==1]
    s2 = sebetahat.ctl^2
    b2 = betahat.ctl^2
    
    if (args$mode == "additive"){
      lambda1 = 1
      tmp = b2-s2
      x = rep(1, length(tmp))
      reg = lm(tmp~0+x)
      lambda2 = reg$coefficients[1]
    }else if (args$mode == "multiplicative"){
      lambda2 = 0
      reg = lm(b2~0+s2)
      lambda1 = reg$coefficients[1]
    }else{
      reg = lm(b2~s2)
      lambda2 = reg$coefficients[1]
      lambda1 = reg$coefficients[2]
    }
    sebetahat.new = sebetahat.norm
    df.new = NULL
    
  }else{
    # directly use t likelihood
    # Suppose betahat/(lambda1*sebetahat+lambda2) ~ t(df)
    # MoM: betahat*sqrt((df-2)/df) = lambda1*sebetahat+lambda2 
    
    betahat.ctl = betahat[ctls==1]
    sebetahat.ctl = sebetahat.new[ctls==1]
    b = betahat.ctl*sqrt((df.new-2)/df.new)
    
    if (args$mode == "additive"){
      lambda1 = 1
      tmp = b-sebetahat.ctl
      x = rep(1, length(tmp))
      reg = lm(tmp~0+x)
      lambda2 = reg$coefficients[1]
    }else if (args$mode == "multiplicative"){
      lambda2 = 0
      reg = lm(b~0+sebetahat.ctl)
      lambda1 = reg$coefficients[1]
    }else{
      reg = lm(b~sebetahat.ctl)
      lambda2 = reg$coefficients[1]
      lambda1 = reg$coefficients[2]
    }
  }
  
  if(args$conservative==TRUE){
    lambda1 = max(1, lambda1)
    lambda2 = max(0, lambda2)
  }
  
  fit = ash(betahat, sqrt(pmax(lambda1*sebetahat.new^2+lambda2,0)), df=df.new)
  
  return(list(fit=fit, lambda1=lambda1, lambda2=lambda2))
}