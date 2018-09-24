sourceDir("datamakers")

nsamp = c(2,4,10)
scale = c(0.25,1,3)
for (i in 1:length(nsamp)){
  addScenario(dsc_gtex,name=paste0("spiky,nsamp=",nsamp[i]), 
              fn=datamaker,
              args=list(ngene=10000, nsamp=nsamp[i], pi0="random",
                        betaargs=list(betapi=c(.4,.2,.2,.2),betamu=c(0,0,0,0),
                                      betasd=c(.25,.5,1,2)/scale[i]),
                        breaksample=TRUE),
              seed=1:1)
}




