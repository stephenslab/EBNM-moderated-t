library("dscr")
library("ashr")
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

load("../dsc-poisthin-dep/res_nullash.RData")

res = separate(res$score,scenario,c("scenario","nsamp"),",nsamp=")
res$nsamp = paste0("N=",res$nsamp)
res$nsamp = factor(res$nsamp, levels=c("N=2","N=10","N=50",
                                       "N=2,alpha=1","N=10,alpha=1","N=50,alpha=1"))
res$scenario = factor(res$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))

#res = filter(res,nsamp %in% c("N=2","N=10","N=50","N=2,alpha=1","N=10,alpha=1","N=50,alpha=1"))
res = filter(res,nsamp %in% c("N=2","N=10","N=50"))
res = filter(res, scenario %in% c("spiky","near_normal","flat_top","big-normal","bimodal"))
#res1 = filter(res, method %in% c("DESeq2","edgeR","voom+limma","voom+vash+ash","logcpm+limma+ash"))
#res2 = filter(res, method %in% c("RUV+DESeq2","RUV+edgeR","RUV+voom+limma","RUV+voom+vash+ash"))
res2 = filter(res, method %in% c("voom+vash+ash","vruv4+ash","ctlinflate.both","ctlinflate.add",
                                 "ctlinflate.multi","ctlinflate.both.1000"))
res1 = filter(res, method %in% c("ctlinflate.both.cons","ctlinflate.add.cons","ctlinflate.multi.cons"))
res3 = filter(res, method %in% c("inflate","inflate.cons"))

setEPS()
postscript('../paper/figures/pi0est_dep_inflate.eps',width=8,height=8)
pi0_plot2=ggplot(res3,
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
postscript('../paper/figures/pi0est_dep_ctlinflate.eps',width=8,height=8)
pi0_plot2=ggplot(res2,
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
postscript('../paper/figures/pi0est_dep_ctlinflate_cons.eps',width=8,height=8)
pi0_plot2=ggplot(res1,
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
postscript('../paper/figures/fdp_dep_ctlinflate.eps',width=8,height=8)
res2$FDP_005[is.na(res2$FDP_005)]=0
fdp_plot=ggplot(res2,
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
postscript('../paper/figures/fsp_dep_ctlinflate.eps',width=8,height=8)
res2$FSP_005[is.na(res2$FSP_005)]=0
fsp_plot=ggplot(res2,
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
postscript('../paper/figures/auc_dep_ctlinflate.eps',width=8,height=8)
fdp_plot=ggplot(res2,
                aes(pi0,AUC,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  xlab("True pi0") +
  ylab("AUC") 
print(fdp_plot +scale_y_continuous(limits=c(0.4,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

################
filter(res2,scenario=="near_normal" & method=="vruv4+ash" & FDP_005>0.15)
filter(res2,scenario=="near_normal" & method=="voom+vash+ash" & FDP_005>0.5)
filter(res2[,1:9], seed==19 & scenario=="near_normal")

nullash.out=readRDS("./dsc-gtex-files/output/near_normal,nsamp=50/nullash/jointash_output/output.19.rds")
nullash.out$nullfit$fitted_g
p.null = 2*(1-pnorm(abs(nullash.out$nullfit$data$x/nullash.out$nullfit$data$s)))
hist(p.null,20)

ash.out = readRDS("./dsc-gtex-files/output/near_normal,nsamp=50/voom+vash+ash/jointash_output/output.19.rds")
ash.out$fit$fitted_g
nullash.out$fit$fitted_g

sum(ash.out$fit$result$qvalue<0.05)
sum(nullash.out$fit$result$qvalue<0.05)

vruvout = readRDS("./dsc-gtex-files/output/near_normal,nsamp=50/vruv4+ash/jointash_output/output.19.rds")
vruvout$fit$ruv4$multiplier
#######################
load("../dsc-poisthin-dep/res.RData")

res = separate(res,scenario,c("scenario","nsamp"),",nsamp=")
res$nsamp = paste0("N=",res$nsamp)
res$nsamp = factor(res$nsamp, levels=c("N=2","N=10","N=50"))
res$scenario = factor(res$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))

res = filter(res,nsamp %in% c("N=2","N=10","N=50"))
res = filter(res, scenario %in% c("spiky","near_normal","flat_top","big-normal","bimodal"))
#res1 = filter(res, method %in% c("DESeq2","edgeR","voom+limma","voom+vash+ash","logcpm+limma+ash"))
#res2 = filter(res, method %in% c("RUV+DESeq2","RUV+edgeR","RUV+voom+limma","RUV+voom+vash+ash"))
res2 = filter(res, method %in% c("voom+vash+ash","vruv4+ash","nullash.approx"))

setEPS()
postscript('../paper/figures/pi0est_dep_nullash_approx.eps',width=6,height=6)
pi0_plot2=ggplot(res2,
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
postscript('../paper/figures/fdp_dep_nullash_approx.eps',width=6,height=6)
res2$FDP_005[is.na(res2$FDP_005)]=0
fdp_plot=ggplot(res2,
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
postscript('../paper/figures/fsp_dep_nullash_approx.eps',width=6,height=6)
res2$FSP_005[is.na(res2$FSP_005)]=0
fsp_plot=ggplot(res2,
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
postscript('../paper/figures/auc_dep_nullash_approx.eps',width=6,height=6)
fdp_plot=ggplot(res2,
                aes(pi0,AUC,colour=method)) +geom_point(shape=1) +
  facet_grid(nsamp ~ scenario) + 
  guides(alpha=FALSE) +
  xlab("True pi0") +
  ylab("AUC") 
print(fdp_plot +scale_y_continuous(limits=c(0.4,1)) +
        scale_x_continuous(limits=c(0,1)) +
        coord_equal(ratio=1) + theme(legend.position = "top",axis.text.x = element_text(size = 8,angle=45)))
dev.off()

##########
load("../dsc-poisthin-dep/res_nullash_tmp.RData")

res = separate(res,scenario,c("scenario","nsamp"),",nsamp=")
res$nsamp = paste0("N=",res$nsamp)
res$nsamp = factor(res$nsamp, levels=c("N=2","N=10","N=50",
                                       "N=2,alpha=1","N=10,alpha=1","N=50,alpha=1"))
res$scenario = factor(res$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))

res = filter(res,nsamp %in% c("N=2","N=10","N=50","N=2,alpha=1","N=10,alpha=1","N=50,alpha=1"))
res = filter(res, scenario %in% c("spiky","near_normal","flat_top","big-normal","bimodal"))
#res1 = filter(res, method %in% c("DESeq2","edgeR","voom+limma","voom+vash+ash","logcpm+limma+ash"))
#res2 = filter(res, method %in% c("RUV+DESeq2","RUV+edgeR","RUV+voom+limma","RUV+voom+vash+ash"))
#res2 = filter(res, method %in% c("voom+vash+ash","vruv4+ash","nullash.approx","nullash.approx.alpha=1"))
res2 = res

setEPS()
postscript('../paper/figures/pi0est_dep_nullash_tmp.eps',width=6,height=6)
pi0_plot2=ggplot(res2,
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


#############
load("../dsc-poisthin-dep/res_nullash.RData")
res = separate(res$score_lambda,scenario,c("scenario","nsamp"),",nsamp=")
res$nsamp = paste0("N=",res$nsamp)
res$nsamp = factor(res$nsamp, levels=c("N=2","N=10","N=50"))
res$scenario = factor(res$scenario, levels=c("spiky"))
res = filter(res, method %in% c("ctlinflate.both","ctlinflate.add","ctlinflate.multi"))

setEPS()
postscript('../paper/figures/hist_lambda1.eps',width=6,height=6)
ggplot(res[res$method!="ctlinflate.add",], aes(lambda1,colour=method)) + 
  geom_histogram(bins=50) + facet_grid(nsamp ~ method) 
dev.off()

setEPS()
postscript('../paper/figures/hist_lambda2.eps',width=6,height=6)
ggplot(res[res$method!="ctlinflate.multi",], aes(lambda2,colour=method)) + 
  geom_histogram(bins=100) + facet_grid(nsamp ~ method) + xlim(-3,3) 
dev.off()


#############
load("../dsc-poisthin-dep/res_nullash.RData")
res = separate(res$score_lambda,scenario,c("scenario","nsamp"),",nsamp=")
res$nsamp = paste0("N=",res$nsamp)
res$nsamp = factor(res$nsamp, levels=c("N=2","N=10","N=50"))
res$scenario = factor(res$scenario, levels=c("spiky","near_normal","flat_top","big-normal","bimodal"))
res = filter(res, method %in% c("inflate"))

setEPS()
postscript('../paper/figures/hist_inflate_lambda1.eps',width=10,height=6)
ggplot(res, aes(lambda1,colour=method)) + 
  geom_histogram(bins=50) + geom_vline(xintercept=1) + 
  facet_grid(nsamp ~ scenario) 
dev.off()

setEPS()
postscript('../paper/figures/hist_inflate_lambda2.eps',width=10,height=6)
ggplot(res, aes(lambda2,colour=method)) + 
  geom_histogram(bins=100) + geom_vline(xintercept=0) + 
  facet_grid(nsamp ~ scenario) + xlim(-3,3) 
dev.off()
