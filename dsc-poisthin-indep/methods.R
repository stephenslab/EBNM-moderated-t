# Methods: 
sourceDir("methods")

addMethod(dsc_gtex,name="DESeq2",fn=deseq2.wrapper,outputtype="pval_output")

addMethod(dsc_gtex,name="edgeR",fn=edger.wrapper,outputtype="pval_output",args=list(exacttest=TRUE))

addMethod(dsc_gtex,name="voom+limma",fn=limma.wrapper,outputtype="pval_output",
          args=list(transform="voom",robust=FALSE))

addMethod(dsc_gtex,name="voom+vash+ash",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="voom",singlecomp=TRUE))