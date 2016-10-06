library(dscr)

dsc_gtex = new.dsc("gtex","dsc-gtex-files")
source("scenarios.R")
source("methods.R")
source("score_short.R")
source("choose_ash_parameters.R")

jointash2qval_est =function(output){
  if (class(output)=="list"){
    qvalue = output$fit$result$qvalue
    svalue = qval.from.lfdr(output$fit$result$lfsr)
    pi0 = output$fit$fitted_g$pi[1]
    return(list(qvalue=qvalue, 
                svalue=svalue,
                pi0=pi0,
                beta.est=output$fit$result$PosteriorMean))
  }else{
    return(list(qvalue=NA, qvalue.fsr=NA, pi0=NA, beta.est=NA))
  }
} 
pval2qval_est =function(output){
  if (class(output)=="list"){
    qq = qvalue(output$pvalue[!is.na(output$pvalue)])
    qvalue = rep(NA,length(output$pvalue))
    qvalue[!is.na(output$pvalue)] = qq$qval
    return(list(qvalue=qvalue, pi0=qq$pi0, beta.est=output$beta.est))
  }else{
    return(list(qvalue=NA, pi0=NA, beta.est=NA))
  }
}
addOutputParser(dsc_gtex,"pval2qval",pval2qval_est,"pval_output","qval_output")
addOutputParser(dsc_gtex,"jointash2qval",jointash2qval_est,"jointash_output","qval_output")

addScore(dsc_gtex,score,name="score",outputtype="qval_output")
addScore(dsc_gtex,score3,"cdf_score","jointash_output")
addScore(dsc_gtex,score_neg,"negprob","jointash_output") #just extract the negativeprobs
addScore(dsc_gtex,score_pos,"posprob","jointash_output") #just extracts the positiveprobs


res=run_dsc(dsc_gtex)

save(res,file="res.Rdata")

# library(ggplot2)
# res = separate(res,scenario,c("scenario","nsamp"),",")
# res$nsamp = factor(res$nsamp, levels=c("nsamp=2","nsamp=10","nsamp=50"))
# 
# ggplot(res, aes(pi0,pi0.est,colour=method))+
#   facet_grid(nsamp~scenario) + geom_point(shape=1) +xlim(0,1) +ylim(0,1) + 
#   xlab("true pi0") +ylab("estimated pi0") +geom_abline(slope=1,intercept=0,color=1) + coord_equal()
# 
# ggplot(res, aes(pi0,FDP_005,colour=method)) + 
#   facet_grid(nsamp~scenario) + geom_point(shape=1) +xlim(0,1) +ylim(0,0.2) + 
#   xlab("true pi0") +ylab("false discovery proportion when q=0.05") +geom_abline(slope=0,intercept=0.05,color=1)