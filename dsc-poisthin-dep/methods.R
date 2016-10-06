# Methods: 
sourceDir("methods")

addMethod(dsc_gtex,name="DESeq2",fn=deseq2.wrapper,outputtype="pval_output")

addMethod(dsc_gtex,name="RUV+DESeq2",fn=deseq2.wrapper,outputtype="pval_output",
          args=list(RUV=TRUE))

addMethod(dsc_gtex,name="edgeR",fn=edger.wrapper,outputtype="pval_output",args=list(exacttest=TRUE))

addMethod(dsc_gtex,name="RUV+edgeR",fn=edger.wrapper,outputtype="pval_output",
          args=list(exacttest=FALSE, RUV=TRUE))

addMethod(dsc_gtex,name="SVA+edgeR",fn=edger.wrapper,outputtype="pval_output",
          args=list(exacttest=FALSE, SVA=TRUE))

addMethod(dsc_gtex,name="voom+limma",fn=limma.wrapper,outputtype="pval_output",
          args=list(transform="voom",robust=FALSE))

addMethod(dsc_gtex,name="RUV+voom+limma",fn=limma.wrapper,outputtype="pval_output",
          args=list(transform="RUVvoom",robust=FALSE))

addMethod(dsc_gtex,name="SVA+voom+limma",fn=limma.wrapper,outputtype="pval_output",
          args=list(transform="SVAvoom",robust=FALSE))

addMethod(dsc_gtex,name="voom+vash+ash",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="voom",singlecomp=TRUE))

addMethod(dsc_gtex,name="RUV+voom+vash+ash",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="RUVvoom",singlecomp=TRUE))

addMethod(dsc_gtex,name="SVA+voom+vash+ash",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="SVAvoom",singlecomp=TRUE))

addMethod(dsc_gtex,name="vruv4+ash",fn=vruv4ash.wrapper,outputtype="jointash_output",
          args=list(nctl=100))