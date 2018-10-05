library("dscr")
library("ashr")
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

#######################
# Make plots (dsc with independent genes)
load("../dsc-poisthin-indep/res.RData")
resscore = separate(res$score,scenario,c("scenario","nsamp"),",nsamp=")

resscore$method = replace(resscore$method,resscore$method=="voom+vash+ash","VL+eBayes+ash")
resscore$method = replace(resscore$method,resscore$method=="voom+vash+ash.alpha=1",
                          "VL+eBayes+ash.alpha=1")
resscore$method = replace(resscore$method,resscore$method=="voom+ash","VL+ash")
resscore$method = replace(resscore$method,resscore$method=="voom+limma","VL+eBayes+qval")
resscore$method[resscore$method %in% c("DESeq2","edgeR")] = 
  paste0(resscore$method[resscore$method %in% c("DESeq2","edgeR")],"+qval")

resscore$nsamp = paste0("N=",resscore$nsamp)
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="N=2","2 vs 2")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="N=4","4 vs 4")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="N=10","10 vs 10")
resscore$nsamp = factor(resscore$nsamp, levels=c("2 vs 2","4 vs 4","10 vs 10"))

resscore$scenario = factor(resscore$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))
resscore2 = filter(resscore, method %in% c("DESeq2+qval","VL+eBayes+qval","VL+pval2se+ash",
                                           "VL+eBayes+ash", "VL+eBayes+ash.alpha=1"))
resscore3 = filter(resscore, method %in% c("DESeq2+qval","edgeR+qval","VL+eBayes+qval",
                                           "VL+pval2se+ash"))
resscore = filter(resscore, method %in% c("VL+eBayes+qval", "VL+ash","VL+pval2se+ash",
                                          "VL+eBayes+ash", "VL+eBayes+ash.alpha=1"))


setEPS()
postscript('../paper/figures/pi0est_indep.eps',width=8,height=8)
pi0_plot=ggplot(resscore,
                aes(pi0,pi0.est,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(colour = "black") +
  xlab("True pi0") +
  ylab("Estimated pi0") 
print(pi0_plot +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + 
        guides(fill=guide_legend(nrow=2))+
        theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
        
dev.off()

setEPS()
postscript('../paper/figures/fdp_indep.eps',width=8,height=8)
fdp_plot=ggplot(resscore,
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
postscript('../paper/figures/fdp_qval_indep.eps',width=8,height=8)
fdp_plot=ggplot(resscore3,
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
postscript('../paper/figures/fsp_indep.eps',width=8,height=8)
fsp_plot=ggplot(resscore[resscore$method %in% c("VL+eBayes+ash","VL+eBayes+ash.alpha=1"),],
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
postscript('../paper/figures/dp_indep.eps',width=8,height=8)
dp_plot=ggplot(resscore,
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
newres = resscore2
newres = newres[newres$method=="VL+eBayes+qval",]
newres = cbind(newres[,2:4],newres$rmse.beta,newres$DP_005,newres$AUC)
names(newres)[4:6] = c("rmse.beta_voomlimma","DP_005_voomlimma","AUC_voomlimma")
test = left_join(resscore2,newres,by=c("seed","scenario","nsamp"))
test$rmse.beta_rel = test$rmse.beta/test$rmse.beta_voomlimma 
test$AUC_rel = test$AUC-test$AUC_voomlimma
test$DP_005_rel = test$DP_005-test$DP_005_voomlimma
test$method[test$method=="VL+eBayes+qval"] = "VL"
test$method[test$method=="DESeq2+qval"] = "DESeq2"

setEPS()
postscript('../paper/figures/rmse_indep.eps',width=8,height=8)
rmse_plot=ggplot(test,
                aes(pi0,rmse.beta_rel,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0,colour = "black") +
  xlab("True pi0") +
  ylab("Relative RMSE of effect estimates") 
print(rmse_plot +scale_y_continuous(limits=c(0,1.2)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

# setEPS()
# postscript('../paper/figures/auc_indep.eps',width=8,height=8)
# auc_plot=ggplot(test,
#                  aes(pi0,AUC_rel,colour=method)) +
#   #geom_point(shape=1) +
#   geom_smooth() +
#   facet_grid(nsamp ~ scenario) + 
#   guides(alpha=FALSE) +
#   #geom_abline(slope=0,intercept=0,colour = "black",lty=2) +
#   xlab("True pi0") +
#   ylab("Relative AUC") 
# print(auc_plot +scale_y_continuous(limits=c(-0.03,0.03)) +
#         scale_x_continuous(limits=c(0,1)) +
#         coord_equal(ratio=10) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
# dev.off()

# setEPS()
# postscript('../paper/figures/dprel_indep.eps',width=8,height=8)
# dprel_plot=ggplot(test,
#                 aes(pi0,DP_005_rel,colour=method)) +geom_point(shape=1) +
#   facet_grid(nsamp ~ scenario) + 
#   guides(alpha=FALSE) +
#   geom_abline(slope=0,intercept=0,colour = "black") +
#   xlab("True pi0") +
#   ylab("Relative discovery proportion when q=0.05") 
# print(dprel_plot +scale_y_continuous(limits=c(-0.05,1)) +
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
                          c("voom+vash+ash"),"../paper/tables/coverage_neg.tex")
save_latex_coverage_table(reslong.pos %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_pos.tex")
save_latex_coverage_table(reslong %>% group_by(scenario,method) %>% summarize(lt = mean(value<coverthresh)) %>% separate(scenario,c("scenario","nsamp"),",nsamp="),
                          c("voom+vash+ash"),"../paper/tables/coverage_all.tex")


#####################
# Make plots (dsc with dependent genes)
rm(list=ls())
load("../dsc-poisthin-dep/res.RData")

resscore = res$score
resscore$method = replace(resscore$method,resscore$method=="voom+vash+ash","VL+eBayes+ash")
resscore$method = replace(resscore$method,resscore$method=="voom+vash+ash.alpha=1",
                          "VL+eBayes+ash.alpha=1")
resscore$method = replace(resscore$method,resscore$method=="voom+ash","VL+ash")
resscore$method = replace(resscore$method,resscore$method=="voom+limma","VL+eBayes+qval")
resscore$method[resscore$method %in% c("DESeq2","edgeR")] = 
  paste0(resscore$method[resscore$method %in% c("DESeq2","edgeR")],"+qval")
resscore$method[resscore$method=="ctlinflate.both.cons"] = "VL+eBayes+ash+inflate"
resscore$method[resscore$method=="inflate.cons"] = "VL+eBayes+ash+inflate.ctl"

resscore = separate(resscore,scenario,c("scenario","nsamp"),",nsamp=")
resscore$nsamp = paste0("N=",resscore$nsamp)
resscore$scenario = factor(resscore$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))

resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="N=2","2 vs 2")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="N=4","4 vs 4")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="N=10","10 vs 10")
resscore$nsamp = factor(resscore$nsamp, levels=c("2 vs 2","4 vs 4","10 vs 10"))

resscore2 = filter(resscore, method %in% c("VL+eBayes+ash","VL+eBayes+ash+inflate","VL+eBayes+ash+inflate.ctl","mouthwash"))
resscore1 = filter(resscore, method %in% c("DESeq2+qval","edgeR+qval","VL+eBayes+qval","VL+eBayes+ash","VL+eBayes+ash.alpha=1"))



setEPS()
postscript('../paper/figures/pi0est_dep.eps',width=10,height=10)
pi0_plot=ggplot(resscore1,
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
postscript('../paper/figures/pi0est2_dep.eps',width=10,height=10)
pi0_plot2=ggplot(resscore2,
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
postscript('../paper/figures/fdp_dep.eps',width=10,height=10)
fdp_plot=ggplot(resscore1,
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
postscript('../paper/figures/fdp2_dep.eps',width=10,height=10)
fdp_plot=ggplot(resscore2,
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
postscript('../paper/figures/fsp_dep.eps',width=10,height=10)
fsp_plot=ggplot(resscore[resscore$method %in% c("voom+limma+ash","voom+limma+ash.alpha=1"),],
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
postscript('../paper/figures/fsp2_dep.eps',width=10,height=10)
fsp_plot=ggplot(resscore2,
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

# setEPS()
# postscript('../paper/figures/auc_dep.eps',width=8,height=8)
# fdp_plot=ggplot(resscore1,
#                 aes(pi0,AUC,colour=method)) +geom_point(shape=1) +
#   geom_smooth(method="loess") +
#   facet_grid(nsamp ~ scenario) + 
#   guides(alpha=FALSE) +
#   xlab("True pi0") +
#   ylab("AUC") 
# print(fdp_plot +scale_y_continuous(limits=c(0.4,1)) +
#         scale_x_continuous(limits=c(0,1)) +
#         coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
# dev.off()

# setEPS()
# postscript('../paper/figures/auc2_dep.eps',width=8,height=8)
# fdp_plot=ggplot(res$score2,
#                 aes(pi0,AUC,colour=method)) +geom_point(shape=1) +
#   facet_grid(nsamp ~ scenario) +
#   guides(alpha=FALSE) +
#   xlab("True pi0") +
#   ylab("AUC")
# print(fdp_plot +scale_y_continuous(limits=c(0.4,1)) +
#         scale_x_continuous(limits=c(0,1)) +
#         coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
# dev.off()

###
newres = resscore1ß
newres = newres[newres$method=="VL+eBayes+qval",]
newres = cbind(newres[,2:4],newres$rmse.beta,newres$DP_005,newres$AUC)
names(newres)[4:6] = c("rmse.beta_voomlimma","DP_005_voomlimma","AUC_voomlimma")
test = left_join(resscore1,newres,by=c("seed","scenario","nsamp"))
test$rmse.beta_rel = test$rmse.beta/test$rmse.beta_voomlimma 
test$AUC_rel = test$AUC-test$AUC_voomlimma
test$DP_005_rel = test$DP_005-test$DP_005_voomlimma
test$method[test$method=="VL+eBayes+qval"] = "VL"
test$method[test$method=="DESeq2+qval"] = "DESeq2"
test$method[test$method=="edgeR+qval"] = "edgeR"

setEPS()
postscript('../paper/figures/rmse_dep.eps',width=8,height=8)
rmse_plot=ggplot(test,
                 aes(pi0,rmse.beta_rel,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  geom_abline(slope=0,intercept=0,colour = "black") +
  xlab("True pi0") +
  ylab("Relative RMSE of effect estimates") 
print(rmse_plot +scale_y_continuous(limits=c(0,1.3)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

# setEPS()
# postscript('../paper/figures/auc_dep.eps',width=8,height=8)
# auc_plot=ggplot(test,
#                 aes(pi0,AUC_rel,colour=method)) +
#   #geom_point(shape=1) +
#   geom_smooth(method="loess") +
#   facet_grid(nsamp ~ scenario) + 
#   guides(alpha=FALSE) +
#   #geom_abline(slope=0,intercept=0,colour = "black") +
#   xlab("True pi0") +
#   ylab("Relative AUC") 
# print(auc_plot +scale_y_continuous(limits=c(-0.03,0.03)) +
#         scale_x_continuous(limits=c(0,1)) +
#         coord_equal(ratio=10) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
# dev.off()


# test2 = left_join(res$score2,newres,by=c("seed","scenario","nsamp"))
# test2$rmse.beta_rel = test2$rmse.beta/test2$rmse.beta_voomlimma 
# test2$AUC_rel = test2$AUC-test2$AUC_voomlimma
# test2$DP_005_rel = test2$DP_005-test2$DP_005_voomlimma
# 
# setEPS()
# postscript('../paper/figures/rmse2_dep.eps',width=8,height=8)
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


##################################
load("../../dscr-gtex-indep-null/res_indep.RData")
resscore = separate(res$score_ks,scenario,c("scenario","nsamp"),",")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="2vs2","2 vs 2")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="4vs4","4 vs 4")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="10vs10","10 vs 10")
resscore$nsamp = factor(resscore$nsamp, levels=c("2 vs 2","4 vs 4","10 vs 10"))
resscore$method = replace(res$score$method,res$score$method=="voom+limma","VL+eBayes")


# setEPS()
# postscript('../paper/figures/pvalks_indep_null.eps',width=8,height=8/3)
# pval_ks_plot = ggplot(resscore, aes(method,pval.ks)) +geom_boxplot() +
#   facet_grid(~nsamp) + 
#   geom_abline(slope=0,intercept=0.05,colour="red") +
#   xlab("Methods") +
#   ylab("P-value of KS test") 
# print(pval_ks_plot +scale_y_continuous(limits=c(0,0.5)))
# dev.off()

setEPS()
postscript('../paper/figures/npval005_indep_null.eps',width=8,height=8/3)
npval_005_plot = ggplot(resscore, aes(method,npval.005)) +geom_boxplot() +
  facet_grid(~nsamp) + 
  geom_abline(slope=0,intercept=0.05,colour="red") +
  xlab("Methods") +
  ylab("Proportion of genes with p<0.05") 
print(npval_005_plot +scale_y_continuous(limits=c(0,0.2)))
dev.off()

# setEPS()
# postscript('../paper/figures/npval_indep_null.eps',width=8,height=8)
# tmp = gather(resscore, "thres", "npval", 7:11)
# thres = unlist(strsplit(tmp$thres,"[.]"))
# thres = as.numeric(paste0('0.',thres[seq(2,length(thres),by=2)]))*10
# tmp$thres = thres
# tmp = tmp[(tmp$thres!=0.0001 & tmp$thres!=0.05),]
# 
# npval_plot = ggplot(tmp, aes(thres,npval,group=thres)) +geom_boxplot(width=0.02) +
#   facet_grid(nsamp~method) + geom_abline(slope=1,intercept=0,colour="red",linetype=2)+
#   scale_x_continuous(breaks=c(0.01,0.05,0.1))+
#   xlab("threshold") +
#   ylab("Proportion of genes with p<threshold") 
# print(npval_plot)
# dev.off()



npvaltab = resscore[,c(1:5,7:11)] %>% gather(thres, npval, 6:10)
npvaltab$thres = as.numeric(substring(npvaltab$thres, 6))*10
npvalsum = npvaltab %>% group_by(nsamp, method, thres) %>% 
  summarise(mean = mean(npval), upp=quantile(npval,0.975), low=quantile(npval,0.025))
npvalsum = npvalsum[(npvalsum$thres!=0.05 & npvalsum$thres!=0.0001),]

setEPS()
postscript('../paper/figures/npval_indep_null.eps',width=6,height=6)
npval_plot = ggplot(npvalsum, aes(x=log10(thres), y=mean/thres, colour=method)) +
  facet_grid(~nsamp) +
  geom_abline(slope=0,intercept=1,linetype="dashed") +
  geom_errorbar(width=0.3,aes(ymin=low/thres, ymax=upp/thres)) +
  geom_point()+
  xlab("p-value threshold") +
  ylab("actual proportion / theoretical proportion")
print(npval_plot +
        scale_x_continuous(breaks=log10(c(0.001,0.01,0.1)), labels=c(0.001,0.01,0.1)) +
        scale_y_continuous(breaks=c(0,1,5,10,20,30))+
        theme(legend.position = "top"))
dev.off()

  
resscore2 = separate(res$score,scenario,c("scenario","nsamp"),",")
resscore2$nsamp = replace(resscore2$nsamp,resscore2$nsamp=="2vs2","2 vs 2")
resscore2$nsamp = replace(resscore2$nsamp,resscore2$nsamp=="4vs4","4 vs 4")
resscore2$nsamp = replace(resscore2$nsamp,resscore2$nsamp=="10vs10","10 vs 10")
resscore2$nsamp = factor(resscore2$nsamp, levels=c("2 vs 2","4 vs 4","10 vs 10"))

setEPS()
postscript('../paper/figures/dp_indep_null.eps',width=8,height=8/3)
dp_plot = ggplot(resscore2, aes(method,DP_005)) +geom_boxplot() +
  facet_grid(~nsamp) + 
  xlab("Methods") +
  ylab("Proportion of genes with q<0.05") 
print(dp_plot +scale_y_continuous(limits=c(0,0.05)))
dev.off()

setEPS()
postscript('../paper/figures/pi0est_indep_null.eps',width=8,height=8/3)
pi0_plot = ggplot(resscore2, aes(method,pi0.est)) +geom_boxplot() +
  facet_grid(~nsamp) + 
  xlab("Methods") +
  ylab("Estimated pi0") 
print(pi0_plot +scale_y_continuous(limits=c(0,1)))
dev.off()


#########################
##########################
load("../../dscr-gtex-indep-null/res_dep.RData")
resscore = separate(res$score_ks,scenario,c("scenario","nsamp"),",")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="2vs2","2 vs 2")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="4vs4","4 vs 4")
resscore$nsamp = replace(resscore$nsamp,resscore$nsamp=="10vs10","10 vs 10")
resscore$nsamp = factor(resscore$nsamp, levels=c("2 vs 2","4 vs 4","10 vs 10"))
resscore$method = replace(res$score$method,res$score$method=="voom+limma","VL+eBayes")

npvaltab = resscore[,c(1:5,7:11)] %>% gather(thres, npval, 6:10)
npvaltab$thres = as.numeric(substring(npvaltab$thres, 6))*10
npvalsum = npvaltab %>% group_by(nsamp, method, thres) %>% 
  summarise(mean = mean(npval), upp=quantile(npval,0.975), low=quantile(npval,0.025))
npvalsum = npvalsum[(npvalsum$thres!=0.05 & npvalsum$thres!=0.0001),]


setEPS()
postscript('../paper/figures/npval_dep_null.eps',width=6,height=6)
npval_plot = ggplot(npvalsum, aes(x=log10(thres), y=mean/thres, colour=method)) +
  facet_grid(~nsamp) +
  geom_abline(slope=0,intercept=1,linetype="dashed") +
  geom_errorbar(width=0.3,aes(ymin=low/thres, ymax=upp/thres)) +
  geom_point()+
  xlab("p-value threshold") +
  ylab("actual proportion / theoretical proportion")
print(npval_plot +
        scale_x_continuous(breaks=log10(c(0.001,0.01,0.1)), labels=c(0.001,0.01,0.1)) +
        scale_y_continuous(breaks=c(0,1,5,10,20,30))+
        theme(legend.position = "top"))
dev.off()

setEPS()
postscript('../paper/figures/pvalks_dep_null.eps',width=8,height=8/3)
pval_ks_plot = ggplot(resscore, aes(method,pval.ks)) +geom_boxplot() +
  facet_grid(~nsamp) + 
  geom_abline(slope=0,intercept=0.05,colour="red") +
  xlab("Methods") +
  ylab("P-value of KS test") 
print(pval_ks_plot +scale_y_continuous(limits=c(0,0.5)))
dev.off()

setEPS()
postscript('../paper/figures/npval005_dep_null.eps',width=8,height=8/3)
npval_005_plot = ggplot(resscore, aes(method,npval.005)) +geom_boxplot() +
  facet_grid(~nsamp) + 
  geom_abline(slope=0,intercept=0.05,colour="red") +
  xlab("Methods") +
  ylab("Proportion of genes with p<0.05") 
print(npval_005_plot +scale_y_continuous(limits=c(0,0.5)))
dev.off()

setEPS()
postscript('../paper/figures/npval_dep_null.eps',width=8,height=8)
tmp = gather(resscore, "thres", "npval", 7:9)
thres = unlist(strsplit(tmp$thres,"[.]"))
thres = as.numeric(paste0('0.',thres[seq(2,length(thres),by=2)]))*10
tmp$thres = thres

npval_plot = ggplot(tmp, aes(thres,npval,group=thres)) +geom_boxplot(width=0.02) +
  facet_grid(nsamp~method) + geom_abline(slope=1,intercept=0,colour="red",linetype=2)+
  scale_x_continuous(breaks=c(0.01,0.05,0.1))+
  xlab("threshold") +
  ylab("Proportion of genes with p<threshold") 
print(npval_plot)
dev.off()

npvaltab = resscore[,c(1:5,7:9)] %>% gather(thres, npval, 6:8)
npvaltab$thres = as.numeric(substring(npvaltab$thres, 6))*10
npvalsum = npvaltab %>% group_by(nsamp, method, thres) %>% 
  summarise(mean = mean(npval), upp=quantile(npval,0.975), low=quantile(npval,0.025))

setEPS()
postscript('../paper/figures/npval_dep_null.eps',width=8,height=4)
npval_plot = ggplot(npvalsum, aes(x=thres, y=mean, colour=method)) +
  facet_grid(~nsamp) +
  geom_abline(slope=1,intercept=0,linetype="dashed") +
  geom_errorbar(width=0.015,aes(ymin=low, ymax=upp)) +
  geom_point()+
  xlab("Threshold") +
  ylab("Proportion of p-values under threshold")
print(npval_plot +
        scale_y_continuous(limits=c(0,0.35)) +
        scale_x_continuous(limits=c(0,0.12), breaks=c(0.01,0.05,0.1), labels=c(0.01,0.05,0.1)) +
        coord_equal(ratio=1) + 
        theme(legend.position = "top"))
dev.off()


resscore2 = separate(res$score,scenario,c("scenario","nsamp"),",")
resscore2$nsamp = replace(resscore2$nsamp,resscore2$nsamp=="2vs2","2 vs 2")
resscore2$nsamp = replace(resscore2$nsamp,resscore2$nsamp=="4vs4","4 vs 4")
resscore2$nsamp = replace(resscore2$nsamp,resscore2$nsamp=="10vs10","10 vs 10")
resscore2$nsamp = factor(resscore2$nsamp, levels=c("2 vs 2","4 vs 4","10 vs 10"))

setEPS()
postscript('../paper/figures/dp_dep_null.eps',width=8,height=8/3)
dp_plot = ggplot(resscore2, aes(method,DP_005)) +geom_boxplot() +
  facet_grid(~nsamp) + 
  xlab("Methods") +
  ylab("Proportion of genes with q<0.05") 
print(dp_plot +scale_y_continuous(limits=c(0,0.2)))
dev.off()

setEPS()
postscript('../paper/figures/pi0est_dep_null.eps',width=8,height=8/3)
pi0_plot = ggplot(resscore2, aes(method,pi0.est)) +geom_boxplot() +
  facet_grid(~nsamp) + 
  xlab("Methods") +
  ylab("Estimated pi0") 
print(pi0_plot +scale_y_continuous(limits=c(0,1)))
dev.off()


############
df = data.frame()
g = normalmix(pi=c(.4,.2,.2,.2),mean=c(0,0,0,0),sd=c(.25,.5,1,2))
x  = seq(-6,6,length = 100)
y  = as.numeric(dens(g,x))
df = rbind(df,data.frame(x = x,y = y,scenario = "spiky"))

g = normalmix(c(2/3,1/3),c(0,0),c(1,2))
x  = seq(-6,6,length = 100)
y  = as.numeric(dens(g,x))
df = rbind(df,data.frame(x = x,y = y,scenario = "near_normal"))

g = normalmix(rep(1/7,7),c(-1.5,-1,-0.5,0,0.5,1,1.5),rep(0.5,7))
x  = seq(-6,6,length = 100)
y  = as.numeric(dens(g,x))
df = rbind(df,data.frame(x = x,y = y,scenario = "flat_top"))

g = normalmix(c(1),c(0),c(4))
x  = seq(-6,6,length = 100)
y  = as.numeric(dens(g,x))
df = rbind(df,data.frame(x = x,y = y,scenario = "big-normal"))

g = normalmix(c(0.5,0.5),c(-2,2),c(1,1))
x  = seq(-6,6,length = 100)
y  = as.numeric(dens(g,x))
df = rbind(df,data.frame(x = x,y = y,scenario = "bimodal"))

setEPS()
postscript('../paper/figures/scenario_dens.eps',width=10,height=2)
ggplot(df,aes(x = x,y = y)) + 
  geom_line(size = 1.2,linetype = 1) + 
  facet_grid(.~scenario) + 
  ylab("density")
dev.off()
