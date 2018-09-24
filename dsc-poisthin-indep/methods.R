# Methods: 
sourceDir("methods")

addMethod(dsc_gtex,name="DESeq2",fn=deseq2.wrapper,outputtype="pval_output")

addMethod(dsc_gtex,name="edgeR",fn=edger.wrapper,outputtype="pval_output",args=list(exacttest=TRUE))

addMethod(dsc_gtex,name="VL+eBayes+qval",fn=limma.wrapper,outputtype="pval_output",
          args=list(transform="voom",robust=FALSE))

addMethod(dsc_gtex,name="VL+ash",fn=ash.wrapper,outputtype="jointash_output",
          args=list(transform="voom"))

addMethod(dsc_gtex,name="VL+eBayes+ash",fn=limmaash.wrapper,outputtype="jointash_output",
          args=list(alpha=0))

addMethod(dsc_gtex,name="VL+eBayes+ash.alpha=1",fn=limmaash.wrapper,outputtype="jointash_output",
          args=list(alpha=1))

addMethod(dsc_gtex,name="VL+pval2se+ash",fn=adhocash.wrapper,outputtype="jointash_output",
          args=list(transform="voom"))

# addMethod(dsc_gtex,name="logcpm+limma+ash",fn=noweightash.wrapper,outputtype="jointash_output")