library("dscr")
library("ashr")
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

#######################
# Make plots (dsc with independent genes)
load("../dsc-poisthin-indep/res.RData")
res$score = separate(res$score,scenario,c("scenario","nsamp"),",nsamp=")
res$score$nsamp = paste0("N=",res$score$nsamp)
res$score$nsamp = factor(res$score$nsamp, levels=c("N=2","N=10","N=50"))
res$score$scenario = factor(res$score$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))

setEPS()
postscript('../paper/figures/pi0est_indep.eps',width=6,height=6)
pi0_plot=ggplot(res$score,
                aes(pi0,pi0.est,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(colour = "black") +
  xlab("True pi0") +
  ylab("Estimated pi0") 
print(pi0_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/fdp_indep.eps',width=6,height=6)
res$score$FDP_005[is.na(res$score$FDP_005)]=0
fdp_plot=ggplot(res$score,
                aes(pi0,FDP_005,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0.05,colour = "black") +
  xlab("True pi0") +
  ylab("False discovery proportion when q=0.05") 
print(fdp_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/fsp_indep.eps',width=6,height=6)
res$score$FSP_005[is.na(res$score$FSP_005)]=0
fsp_plot=ggplot(res$score[res$score$method=="voom+vash+ash",],
                aes(pi0,FSP_005,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0.05,colour = "black") +
  xlab("True pi0") +
  ylab("False sign proportion when s=0.05") 
print(fsp_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/dp_indep.eps',width=6,height=6)
dp_plot=ggplot(res$score,
                aes(pi0,DP_005,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=-1,intercept=1,colour = "black") +
  xlab("True pi0") +
  ylab("Discovery proportion when q=0.05") 
print(dp_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

####
newres = res$score
newres = newres[newres$method=="voom+limma",]
newres = cbind(newres[,2:4],newres$rmse.beta,newres$DP_005,newres$AUC)
names(newres)[4:6] = c("rmse.beta_voomlimma","DP_005_voomlimma","AUC_voomlimma")
test = left_join(res$score,newres,by=c("seed","scenario","nsamp"))
test$rmse.beta_rel = test$rmse.beta/test$rmse.beta_voomlimma 
test$AUC_rel = test$AUC-test$AUC_voomlimma
test$DP_005_rel = test$DP_005-test$DP_005_voomlimma

setEPS()
postscript('../paper/figures/rmse_indep.eps',width=6,height=6)
rmse_plot=ggplot(test,
                aes(pi0,rmse.beta_rel,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=1,colour = "black") +
  xlab("True pi0") +
  ylab("Relative RMSE of effect estimates") 
print(rmse_plot +scale_y_continuous(limits=c(0,1.2)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/auc_indep.eps',width=6,height=6)
auc_plot=ggplot(test,
                 aes(pi0,AUC_rel,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0,colour = "black") +
  xlab("True pi0") +
  ylab("Relative AUC") 
print(auc_plot +scale_y_continuous(limits=c(-0.03,0.03)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=10) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/dprel_indep.eps',width=6,height=6)
dprel_plot=ggplot(test,
                aes(pi0,DP_005_rel,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0,colour = "black") +
  xlab("True pi0") +
  ylab("Relative discovery proportion when q=0.05") 
print(dprel_plot +scale_y_continuous(limits=c(-0.05,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()


# Coverage tables
coverthresh = 0.05 # threshold at which we look at coverage
findthresh=0.95 #threshold at we define a discovery significant

neglong = 
  res$negprob %>% 
  select(-user.self,-sys.self,-elapsed,-user.child,-sys.child) %>%
  melt(id.vars=c("method","scenario","seed",".id"),value.name="negprob") %>%
  filter(negprob > findthresh) %>%
  filter(method=="voom+vash+ash")

poslong = 
  res$posprob %>% 
  select(-user.self,-sys.self,-elapsed,-user.child,-sys.child) %>%
  melt(id.vars=c("method","scenario","seed",".id"),value.name="posprob") %>%
  filter(posprob > findthresh) %>%
  filter(method=="voom+vash+ash")

reslong = 
  res$cdf_score %>% 
  select(-user.self,-sys.self,-elapsed,-user.child,-sys.child) %>%    
  melt(id.vars=c("method","scenario","seed",".id")) %>%
  filter(method=="voom+vash+ash")

reslong.pos = inner_join(reslong,poslong)
reslong.neg = inner_join(reslong,neglong)
dim(reslong.pos)
dim(reslong.neg)

xtabs(lt~as.numeric(nsamp)+scenario,reslong %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp=")) %>% round(2)
xtabs(lt~as.numeric(nsamp)+scenario,reslong.pos %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp=")) %>% round(2)
xtabs(lt~as.numeric(nsamp)+scenario,reslong.neg %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp=")) %>% round(2)


save_latex_coverage_table=function(df,methodnames,filename,SCENARIONAMES=c("spiky","near_normal","flat_top","big-normal","bimodal"),switch=TRUE){
  library(xtable)
  df$method = factor(df$method,levels=methodnames)
  df$scenario=factor(df$scenario,levels=SCENARIONAMES)
  df$scenario = gsub("_","-",df$scenario)
  mat <- as.matrix(xtabs(lt~as.numeric(nsamp)+scenario,df))
  if(switch){
    mat=1-mat
  }
  mat <- xtable(mat,digits=rep(2,ncol(mat)+1))
  rownames(mat)=paste0("N=",rownames(mat))
  write(print(mat, 
              sanitize.text.function = function(x){x},
              floating=FALSE, 
              hline.after=NULL, 
              add.to.row=list(pos=list(-1,0, nrow(mat)), 
                              command=c('\\toprule ',
                                        '\\midrule ',
                                        '\\bottomrule '))
  ),file=filename)
  
}

save_latex_coverage_table(reslong.neg %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_neg.tex")
save_latex_coverage_table(reslong.pos %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_pos.tex")
save_latex_coverage_table(reslong %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_all.tex")


#####################
# Make plots (dsc with dependent genes)
rm(list=ls())
load("../dsc-poisthin-dep/res.RData")

res$score = separate(res$score,scenario,c("scenario","nsamp"),",nsamp=")
res$score$nsamp = paste0("N=",res$score$nsamp)
res$score$nsamp = factor(res$score$nsamp, levels=c("N=2","N=10","N=50"))
res$score$scenario = factor(res$score$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))

res$score = filter(res$score,nsamp %in% c("N=2","N=10","N=50"))
res$score = filter(res$score, scenario %in% c("spiky","near_normal","flat_top","big-normal","bimodal"))
res$score1 = filter(res$score, method %in% c("DESeq2","edgeR","voom+limma","voom+vash+ash"))
#res$score2 = filter(res$score, method %in% c("RUV+DESeq2","RUV+edgeR","RUV+voom+limma","RUV+voom+vash+ash"))
res$score2 = filter(res$score, method %in% c("voom+vash+ash","vruv4+ash"))

setEPS()
postscript('../paper/figures/pi0est_dep.eps',width=6,height=6)
pi0_plot=ggplot(res$score1,
                aes(pi0,pi0.est,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(colour = "black") +
  xlab("True pi0") +
  ylab("Estimated pi0") 
print(pi0_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/pi0est2_dep.eps',width=6,height=6)
pi0_plot2=ggplot(res$score2,
                aes(pi0,pi0.est,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(colour = "black") +
  xlab("True pi0") +
  ylab("Estimated pi0") 
print(pi0_plot2 +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/fdp_dep.eps',width=6,height=6)
res$score1$FDP_005[is.na(res$score1$FDP_005)]=0
fdp_plot=ggplot(res$score1,
                aes(pi0,FDP_005,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0.05,colour = "black") +
  xlab("True pi0") +
  ylab("False discovery proportion when q=0.05") 
print(fdp_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/fdp2_dep.eps',width=6,height=6)
res$score2$FDP_005[is.na(res$score2$FDP_005)]=0
fdp_plot=ggplot(res$score2,
                aes(pi0,FDP_005,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0.05,colour = "black") +
  xlab("True pi0") +
  ylab("False discovery proportion when q=0.05") 
print(fdp_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/fsp_dep.eps',width=6,height=6)
res$score$FSP_005[is.na(res$score$FSP_005)]=0
fsp_plot=ggplot(res$score[res$score$method=="voom+vash+ash",],
                aes(pi0,FSP_005,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0.05,colour = "black") +
  xlab("True pi0") +
  ylab("False sign proportion when s=0.05") 
print(fsp_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

setEPS()
postscript('../paper/figures/fsp2_dep.eps',width=6,height=6)
res$score2$FSP_005[is.na(res$score2$FSP_005)]=0
fsp_plot=ggplot(res$score2[res$score2$method %in% c("voom+vash+ash","vruv4+ash"),],
                aes(pi0,FSP_005,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0.05,colour = "black") +
  xlab("True pi0") +
  ylab("False sign proportion when s=0.05") 
print(fsp_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

###
newres = res$score1
newres = newres[newres$method=="voom+limma",]
newres = cbind(newres[,2:4],newres$rmse.beta,newres$DP_005,newres$AUC)
names(newres)[4:6] = c("rmse.beta_voomlimma","DP_005_voomlimma","AUC_voomlimma")
test = left_join(res$score1,newres,by=c("seed","scenario","nsamp"))
test$rmse.beta_rel = test$rmse.beta/test$rmse.beta_voomlimma 
test$AUC_rel = test$AUC-test$AUC_voomlimma
test$DP_005_rel = test$DP_005-test$DP_005_voomlimma

setEPS()
postscript('../paper/figures/rmse_dep.eps',width=6,height=6)
rmse_plot=ggplot(test,
                 aes(pi0,rmse.beta_rel,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=1,colour = "black") +
  xlab("True pi0") +
  ylab("Relative RMSE of effect estimates") 
print(rmse_plot +scale_y_continuous(limits=c(0,1.3)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

# test2 = left_join(res$score2,newres,by=c("seed","scenario","nsamp"))
# test2$rmse.beta_rel = test2$rmse.beta/test2$rmse.beta_voomlimma 
# test2$AUC_rel = test2$AUC-test2$AUC_voomlimma
# test2$DP_005_rel = test2$DP_005-test2$DP_005_voomlimma
# 
# setEPS()
# postscript('../paper/figures/rmse2_dep.eps',width=6,height=6)
# rmse_plot=ggplot(test2,
#                  aes(pi0,rmse.beta_rel,colour=method)) +geom_point(shape=1) +
#   facet_grid(nsamp ~ scenario) + 
#   guides(alpha=FALSE) +
#   geom_abline(slope=0,intercept=1,colour = "black") +
#   xlab("True pi0") +
#   ylab("Relative RMSE of effect estimates") 
# print(rmse_plot +scale_y_continuous(limits=c(0,1.3)) +
#         scale_x_continuous(limits=c(0,1)) +
#         coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
# dev.off()


# Coverage tables
coverthresh = 0.05 # threshold at which we look at coverage
findthresh=0.95 #threshold at we define a discovery significant

neglong = 
  res$negprob %>% 
  select(-user.self,-sys.self,-elapsed,-user.child,-sys.child) %>%
  melt(id.vars=c("method","scenario","seed",".id"),value.name="negprob") %>%
  filter(negprob > findthresh) %>%
  filter(method=="voom+vash+ash")

poslong = 
  res$posprob %>% 
  select(-user.self,-sys.self,-elapsed,-user.child,-sys.child) %>%
  melt(id.vars=c("method","scenario","seed",".id"),value.name="posprob") %>%
  filter(posprob > findthresh) %>%
  filter(method=="voom+vash+ash")

reslong = 
  res$cdf_score %>% 
  select(-user.self,-sys.self,-elapsed,-user.child,-sys.child) %>%    
  melt(id.vars=c("method","scenario","seed",".id")) %>%
  filter(method=="voom+vash+ash")

reslong.pos = inner_join(reslong,poslong)
reslong.neg = inner_join(reslong,neglong)
dim(reslong.pos)
dim(reslong.neg)

xtabs(lt~as.numeric(nsamp)+scenario,reslong %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp=")) %>% round(2)
xtabs(lt~as.numeric(nsamp)+scenario,reslong.pos %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp=")) %>% round(2)
xtabs(lt~as.numeric(nsamp)+scenario,reslong.neg %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp=")) %>% round(2)


save_latex_coverage_table=function(df,methodnames,filename,SCENARIONAMES=c("spiky","near_normal","flat_top","big-normal","bimodal"),switch=TRUE){
  library(xtable)
  df$method = factor(df$method,levels=methodnames)
  df$scenario=factor(df$scenario,levels=SCENARIONAMES)
  df$scenario = gsub("_","-",df$scenario)
  mat <- as.matrix(xtabs(lt~as.numeric(nsamp)+scenario,df))
  if(switch){
    mat=1-mat
  }
  mat <- xtable(mat,digits=rep(2,ncol(mat)+1))
  rownames(mat)=paste0("N=",rownames(mat))
  write(print(mat, 
              sanitize.text.function = function(x){x},
              floating=FALSE, 
              hline.after=NULL, 
              add.to.row=list(pos=list(-1,0, nrow(mat)), 
                              command=c('\\toprule ',
                                        '\\midrule ',
                                        '\\bottomrule '))
  ),file=filename)
  
}

save_latex_coverage_table(reslong.neg %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_neg_dep.tex")
save_latex_coverage_table(reslong.pos %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_pos_dep.tex")
save_latex_coverage_table(reslong %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_all_dep.tex")
