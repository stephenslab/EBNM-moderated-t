library(DESeq2)
library(qvalue)
deseq2.wrapper = function(input,args){

  if(is.null(args$pseudocounts)){
    args$pseudocounts = 1
  }
  cond = input$condition
  dds = DESeqDataSetFromMatrix(input$counts+args$pseudocounts, DataFrame(cond), ~cond)
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  dds = nbinomWaldTest(dds)
  res = results(dds,cooksCutoff=FALSE,independentFiltering=FALSE)
  beta.est = res$log2FoldChange
  pvalue = res$pvalue
  
  return(list(pvalue=pvalue, beta.est=beta.est))
}