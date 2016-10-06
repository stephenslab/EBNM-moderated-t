library(DESeq2)
library(qvalue)
deseq2.wrapper = function(input,args){

  if(is.null(args$pseudocounts)){
    args$pseudocounts = 1
  }
  cond = input$condition
  
  if(is.null(args$RUV)){
    args$RUV = FALSE
  }
  if(is.null(args$SVA)){
    args$SVA = FALSE
  }
  if(args$RUV==TRUE & !is.null(input$W.RUV)){
    vars = data.frame(input$W.RUV,cond)
    dds = DESeqDataSetFromMatrix(input$counts+args$pseudocounts, colData=vars, 
                                 as.formula(paste0("~",paste0(names(vars),collapse = "+"))))
  }else if(args$SVA==TRUE & !is.null(input$W.SVA)){
    vars = data.frame(input$W.SVA,cond)
    dds = DESeqDataSetFromMatrix(input$counts+args$pseudocounts, colData=vars, 
                                   as.formula(paste0("~",paste0(names(vars),collapse = "+"))))        
  }else{
    dds = DESeqDataSetFromMatrix(input$counts+args$pseudocounts, DataFrame(cond), ~cond)
  }
  
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  dds = nbinomWaldTest(dds)
  res = results(dds,cooksCutoff=FALSE,independentFiltering=FALSE)
  beta.est = res$log2FoldChange
  pvalue = res$pvalue
  
  return(list(pvalue=pvalue, beta.est=beta.est))
}