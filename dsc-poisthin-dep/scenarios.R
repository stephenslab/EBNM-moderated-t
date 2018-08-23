sourceDir("datamakers")


nsamp = c(2,4,10)
#scale = c(0.25,1,3)
scale = c(0.25,1,3)/2

for (i in 1:length(nsamp)){
  addScenario(dsc_gtex,name=paste0("spiky,nsamp=",nsamp[i]), 
              fn=datamaker,
              args=list(ngene=10000, nsamp=nsamp[i], pi0="random",
                        betaargs=list(betapi=c(.4,.2,.2,.2),
                                      betamu=c(0,0,0,0),betasd=c(.25,.5,1,2)/scale[i]),
                        breaksample=FALSE, nctl=100),
              seed=1:50)
  
  addScenario(dsc_gtex,name=paste0("near_normal,nsamp=",nsamp[i]), 
              fn=datamaker,
              args=list(ngene=10000, nsamp=nsamp[i], pi0="random",
                        betaargs=list(betapi=c(2/3,1/3),
                                      betamu=c(0,0),betasd=c(1,2)/scale[i]),
                        breaksample=FALSE, nctl=100),
              seed=1:50)
  
  addScenario(dsc_gtex,name=paste0("flat_top,nsamp=",nsamp[i]), 
              fn=datamaker,
              args=list(ngene=10000, nsamp=nsamp[i], pi0="random",
                        betaargs=list(betapi=rep(1/7,7),
                                      betamu=c(-1.5,-1,-0.5,0,0.5,1,1.5),betasd=rep(0.5,7)/scale[i]),
                        breaksample=FALSE, nctl=100),
              seed=1:50)
  
  addScenario(dsc_gtex,name=paste0("big-normal,nsamp=",nsamp[i]), 
              fn=datamaker,
              args=list(ngene=10000, nsamp=nsamp[i], pi0="random",
                        betaargs=list(betapi=c(1),
                                      betamu=c(0),betasd=c(4)/scale[i]),
                        breaksample=FALSE, nctl=100),
              seed=1:50)
  
  addScenario(dsc_gtex,name=paste0("bimodal,nsamp=",nsamp[i]), 
              fn=datamaker,
              args=list(ngene=10000, nsamp=nsamp[i], pi0="random",
                        betaargs=list(betapi=c(0.5,0.5),
                                      betamu=c(-2,2),betasd=c(1,1)/scale[i]),
                        breaksample=FALSE, nctl=100),
              seed=1:50) 
}