# Methods: 
sourceDir("methods")

addMethod(dsc_gtex,name="DESeq2",fn=deseq2.wrapper,outputtype="pval_output")

# addMethod(dsc_gtex,name="RUV+DESeq2",fn=deseq2.wrapper,outputtype="pval_output",
#           args=list(RUV=TRUE))

addMethod(dsc_gtex,name="edgeR",fn=edger.wrapper,outputtype="pval_output",args=list(exacttest=TRUE))

# addMethod(dsc_gtex,name="RUV+edgeR",fn=edger.wrapper,outputtype="pval_output",
#           args=list(exacttest=FALSE, RUV=TRUE))
# 
# addMethod(dsc_gtex,name="SVA+edgeR",fn=edger.wrapper,outputtype="pval_output",
#           args=list(exacttest=FALSE, SVA=TRUE))

addMethod(dsc_gtex,name="VL+eBayes+qval",fn=limma.wrapper,outputtype="pval_output",
          args=list(transform="voom",robust=FALSE))

# addMethod(dsc_gtex,name="RUV+voom+limma",fn=limma.wrapper,outputtype="pval_output",
#           args=list(transform="RUVvoom",robust=FALSE))
# 
# addMethod(dsc_gtex,name="SVA+voom+limma",fn=limma.wrapper,outputtype="pval_output",
#           args=list(transform="SVAvoom",robust=FALSE))

addMethod(dsc_gtex,name="VL+ash",fn=ash.wrapper,outputtype="jointash_output",
          args=list(transform="voom"))

addMethod(dsc_gtex,name="VL+eBayes+ash",fn=limmaash.wrapper,outputtype="jointash_output",
          args=list(alpha=0))

addMethod(dsc_gtex,name="VL+eBayes+ash.alpha=1",fn=limmaash.wrapper,outputtype="jointash_output",
          args=list(alpha=1))

# addMethod(dsc_gtex,name="RUV+voom+vash+ash",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="RUVvoom",singlecomp=TRUE))
# 
# addMethod(dsc_gtex,name="SVA+voom+vash+ash",fn=jointash.wrapper,outputtype="jointash_output",
#           args=list(transform="SVAvoom",singlecomp=TRUE))

# addMethod(dsc_gtex,name="vruv4+ash",fn=vruv4ash.wrapper,outputtype="jointash_output",
#           args=list(nctl=100))

# addMethod(dsc_gtex,name="vruv4+ash1000",fn=vruv4ash.wrapper,outputtype="jointash_output",
#           args=list(nctl=1000))

# addMethod(dsc_gtex,name="logcpm+limma+ash",fn=noweightash.wrapper,outputtype="jointash_output")

# addMethod(dsc_gtex,name="nullash",fn=nullash.wrapper,outputtype="jointash_output",
#           args=list(nctl=100))

# addMethod(dsc_gtex,name="nullash.approx",fn=nullash.approx.wrapper,outputtype="jointash_output",
#           args=list(nctl=100, transform="voom"))

# addMethod(dsc_gtex,name="empiricalqval",fn=empiricalqval.wrapper,outputtype="qval_output",
#           args=list(nctl=100))

addMethod(dsc_gtex,name="VL+eBayes+ash+inflate",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=100, transform="voom", mode="both", 
                    conservative=TRUE, normapprox=FALSE))

addMethod(dsc_gtex,name="VL+eBayes+ash+inflate.ctl",fn=inflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(transform="voom", conservative=TRUE))

addMethod(dsc_gtex,name="mouthwash",fn=mouthwash.wrapper,outputtype="jointash_output",
          args=list())