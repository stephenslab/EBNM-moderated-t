library(locfdr)
inflate.wrapper = function(input, args){
  # variance shrinkage
  ash.params = choose.ash.parameters(input, args)
  va = vash(ash.params$sebetahat/ash.params$scale,
            df=ash.params$df,singlecomp=TRUE)
  sebetahat.new = ash.params$scale*c(sqrt(va$PosteriorRate/va$PosteriorShape))
  df.new = 2*va$PosteriorShape[1]
  betahat = ash.params$betahat
  
  # convert t scores to z scores by quantile matching
  # use Efron's empirical method to estimate lambda1
  sebetahat.norm = betahat/qnorm(pt(betahat/sebetahat.new, df=df.new))
  z = betahat/sebetahat.norm
  z[abs(z)==Inf] = (betahat/(sebetahat.new+1e-50))[abs(z)==Inf]
  empnull = locfdr(z)
  lambda1 = empnull$fp0[3,2]
  
  # lambda1 fixed, directly use t likelihood to estimate lambda2
  # Suppose betahat/(lambda1*sebetahat+lambda2) ~ t(df)
  # MoM: betahat*sqrt((df-2)/df) = lambda1*sebetahat+lambda2 
  b = betahat*sqrt((df.new-2)/df.new)
  tmp = b-sebetahat.new*lambda1
  x = rep(1, length(tmp))
  reg = lm(tmp~0+x)
  lambda2 = reg$coefficients[1]
  
  if(args$conservative==TRUE){
    lambda1 = max(1, lambda1)
    lambda2 = max(0, lambda2)
  }
  
  fit = ash(betahat, sqrt(pmax(lambda1*sebetahat.new^2+lambda2,0)), df=df.new)
  
  return(list(fit=fit, lambda1=lambda1, lambda2=lambda2))
}