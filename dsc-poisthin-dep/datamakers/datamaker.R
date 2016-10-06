library(data.table)
library(ashr)
library(vashr)
library(limma)
library(qvalue)
library(edgeR)
library(DESeq2)
library(RUVSeq)
library(sva)
library(ggplot2)
  
datamaker = function(args){  
  #rawdata = read.table("/Volumes/PERSONAL/MS/gtex/tissues/liver.txt",header=TRUE)
  rawdata = as.matrix(fread("/Volumes/PERSONAL/MS/gtex/tissues/Muscle-Skeletal.txt",
                            header=FALSE,drop=1,skip=1))
  colnames(rawdata) = as.character(read.table("/Volumes/PERSONAL/MS/gtex/tissues/Muscle-Skeletal.txt",
                                              nrows=1,stringsAsFactors=FALSE))
  
  #rawdata = rawdata[,-(1:2)]  # remove the first 2 samples (outliers)
  counts = selectsample(rawdata, 2*args$nsamp,args$breaksample)
  
  # Remove genes without any reads
  counts = counts[apply(counts,1,sum)>0,]
  
  # Take the top Ngene high-expressed genes
  counts = counts[sort(order(rowSums(counts),decreasing=TRUE)[1:args$ngene]),]
  
  # Model's design: nsamp samples for group A and nsamp samples for group B
  condition = factor(rep(1:2,each=args$nsamp))
  design = model.matrix(~condition)
  
  # Thinned effect size generated from normal mixture
  be = gen_normalmix(args$ngene, 
                     args$betaargs$betapi, args$betaargs$betamu, args$betaargs$betasd, 
                     args$pi0)
  null = be$null
  
  # Poisson thinning (optional)
  counts = pois_thinning(counts, be$beta)
  
  if (is.null(args$num.sv)){
    args$num.sv = num.sv(as.matrix(log(counts+1)), 
                           mod = stats::model.matrix(~condition), method = "be")
    args$num.sv = min(args$num.sv, 2*args$nsamp-3)
  }
  
  # Voom transformation
  voom = voom_transform(counts, condition)
  
  # RUV & voom
  ctls = rep(0,length(null))  # Use nctl true nulls to do supervised RUV/SVA
  #ctls[which(null==1)[1:floor(length(which(null==1))/2)]] = 1
  ctls[which(null==1)[1:args$nctl]] = 1
  W.RUV = RUV_factor(counts, args, ctls)
  RUVvoom = voom_transform(counts, condition, W=W.RUV)
  
  # SVA & voom
  W.SVA = SVA_factor(counts, condition, args, ctls)
  SVAvoom = voom_transform(counts, condition, W=W.SVA)
  
  meta = list(pi0=be$pi0, null=be$null, beta=be$beta,
              betaprior=list(pi=args$betaargs$betapi,mu=args$betaargs$betamu,sd=args$betaargs$betasd), 
              args=args)
  input = list(counts=counts, condition=condition,
               v=voom$v, RUVv=RUVvoom$v, SVAv=SVAvoom$v, 
               betahat.voom=voom$betahat, sebetahat.voom=voom$sebetahat, df.voom=voom$df,
               betahat.RUVvoom=RUVvoom$betahat, sebetahat.RUVvoom=RUVvoom$sebetahat, df.RUVvoom=RUVvoom$df, W.RUV=W.RUV,
               betahat.SVAvoom=SVAvoom$betahat, sebetahat.SVAvoom=SVAvoom$sebetahat, df.SVAvoom=SVAvoom$df, W.SVA=W.SVA)
  data = list(meta=meta,input=input)
  return(data)
}

# Generate beta from normal mixture prior
gen_normalmix = function(ngene, pi, mu, sd, pi0){
  if (pi0=="random"){
    pi0 = runif(1,0,1) #generate the proportion of true nulls randomly
  }
  k = length(pi) # number of components
  comp = sample(1:k,ngene,pi,replace=TRUE) #randomly draw a component
  isnull = (runif(ngene,0,1) < pi0)
  beta = ifelse(isnull, 0, rnorm(ngene,mu[comp],sd[comp]))
  return(list(beta=beta, pi0=pi0, null=isnull))
}

# Poisson thinning
pois_thinning = function(counts, log2foldchanges){
  nsamp = dim(counts)[2]/2
  null = (log2foldchanges==0)
  log2foldchanges = log2foldchanges[!null]
  foldchanges = 2^log2foldchanges
  
  # thin group A
  counts[which(!null)[log2foldchanges>0],1:nsamp]=matrix(rbinom(sum(log2foldchanges>0)*nsamp, 
                                                                size=c(as.matrix(counts[which(!null)[log2foldchanges>0],1:nsamp])),
                                                                prob=rep(1/foldchanges[log2foldchanges>0],nsamp)),ncol=nsamp)
  # thin group B
  counts[which(!null)[log2foldchanges<0],(nsamp+1):(2*nsamp)]=matrix(rbinom(sum(log2foldchanges<0)*nsamp, 
                                                                            size=c(as.matrix(counts[which(!null)[log2foldchanges<0],
                                                                                                    (nsamp+1):(2*nsamp)])),
                                                                            prob=rep(foldchanges[log2foldchanges<0],nsamp)),
                                                                     ncol=nsamp)
  
  return(counts)
}

# Voom transformation
voom_transform = function(counts, condition, W=NULL){
  dgecounts = calcNormFactors(DGEList(counts=counts,group=condition))
  
  if (is.null(W)){
    design = model.matrix(~condition)
  }else{
    design = model.matrix(~condition+W)
  }
  
  v = voom(dgecounts,design,plot=FALSE)
  lim = lmFit(v)
  
  betahat.voom = lim$coefficients[,2]
  sebetahat.voom = lim$stdev.unscaled[,2]*lim$sigma
  if (!is.null(W)){
    df.voom = length(condition)-2-dim(W)[2]
  }else{
    df.voom = length(condition)-2
  }
  
  scale.voom = lim$stdev.unscaled[,2]
  
  return(list(betahat=betahat.voom, sebetahat=sebetahat.voom, df=df.voom, 
              v=v, scale=scale.voom))
}

# randomly subsample data for each gene
# gene: a vector of reads for one gene
# nsamp: # of samples wanted
sampleingene = function(gene, nsamp){
  sample = sample(length(gene),nsamp)
  return(c(gene[sample]))
}

# Randomly select samples
# counts: full count matrix
# nsamp: # of samples wanted
# breaksample: flag, if select different samples for each gene
selectsample = function(counts, nsamp, breaksample){
  if (breaksample==FALSE){
    subsample = sample(1:dim(counts)[2],nsamp)
    counts = counts[,subsample]
  }else{
    counts = t(apply(counts, 1, sampleingene, nsamp=nsamp))
  }
  return(counts)
}

## Use RUV to estimate confounding factor
RUV_factor = function(counts, args, null) {
  seq = EDASeq::newSeqExpressionSet(as.matrix(counts[as.logical(null), ]))
  if (sum(null)>0 & args$num.sv>0) {
    controls = rownames(seq)
    differences = matrix(data = c(1:args$nsamp, (args$nsamp + 1):(2 * args$nsamp)), 
                         byrow = TRUE, nrow = 2)
    seqRUV = RUVSeq::RUVg(seq, controls, k = min(args$num.sv,sum(null)))
    return(W = as.matrix(Biobase::pData(seqRUV)))
  } else {
    return(W = NULL)
  }
}

# Use SVA to estimate confounding factor
SVA_factor = function(counts, condition, args, null = NULL) {
  mod1 = model.matrix(~condition)
  mod0 = cbind(mod1[, 1]) 
  if (!is.null(null)) {
    svseq_out = sva::svaseq(counts, mod1, mod0, control = null, n.sv = min(args$num.sv,sum(null)))
  } else {
    svseq_out = sva::svaseq(as.matrix(counts), mod1, mod0, n.sv = args$num.sv)
  }
  if (svseq_out$n.sv > 0) {
    W = svseq_out$sv
    W = matrix(W, ncol=svseq_out$n.sv)
    return(W)
  } else {
    return(W = NULL)
  }
}