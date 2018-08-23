# Methods: 
sourceDir("methods")

addMethod(dsc_gtex,name="voom+vash+ash",fn=jointash.wrapper,outputtype="jointash_output",
          args=list(transform="voom",singlecomp=TRUE))

addMethod(dsc_gtex,name="vruv4+ash",fn=vruv4ash.wrapper,outputtype="jointash_output",
          args=list(nctl=100))

# addMethod(dsc_gtex,name="nullash.approx",fn=nullash.approx.wrapper,
#           outputtype="jointash_output",
#           args=list(nctl=100, transform="voom", alpha=0, pointmass=TRUE))

# addMethod(dsc_gtex,name="nullash.approx1",fn=nullash.approx.wrapper,
#           outputtype="jointash_output",
#           args=list(nctl=100, transform="voom", alpha=0, pointmass=FALSE))

# addMethod(dsc_gtex,name="nullash.approx.alpha1",fn=nullash.approx.wrapper,
#           outputtype="jointash_output", args=list(nctl=100, transform="voom", alpha=1))

# addMethod(dsc_gtex,name="nullash.simp",fn=nullash.simp.wrapper,
#           outputtype="jointash_output",args=list(nctl=100, transform="voom", alpha=0))

addMethod(dsc_gtex,name="ctlinflate.both",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=100, transform="voom", mode="both", 
                    conservative=FALSE, normapprox=TRUE))

addMethod(dsc_gtex,name="ctlinflate.add",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=100, transform="voom", mode="additive", 
                    conservative=FALSE, normapprox=TRUE))

addMethod(dsc_gtex,name="ctlinflate.multi",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=100, transform="voom", mode="multiplicative", 
                    conservative=FALSE, normapprox=TRUE))

addMethod(dsc_gtex,name="ctlinflate.both.cons",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=100, transform="voom", mode="both", 
                    conservative=TRUE, normapprox=TRUE))

addMethod(dsc_gtex,name="ctlinflate.add.cons",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=100, transform="voom", mode="additive", 
                    conservative=TRUE, normapprox=TRUE))

addMethod(dsc_gtex,name="ctlinflate.multi.cons",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=100, transform="voom", mode="multiplicative", 
                    conservative=TRUE, normapprox=TRUE))

addMethod(dsc_gtex,name="ctlinflate.both.1000",fn=ctlinflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(nctl=1000, transform="voom", mode="both", 
                    conservative=FALSE, normapprox=FALSE))

addMethod(dsc_gtex,name="inflate",fn=inflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(transform="voom", conservative=FALSE))

addMethod(dsc_gtex,name="inflate.cons",fn=inflate.wrapper,
          outputtype="ctlinflate_output",
          args=list(transform="voom", conservative=TRUE))