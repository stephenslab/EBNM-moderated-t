# EBNM

This repository contains the codes and analyses for the "EBNM" project (Empirical Bayes estimation of normal means, accounting for uncertainty in estimated standard errors).

Vignette for our proposed method *VL+eBayes+ash*: `analysis/method.Rmd`.

For simulation studies, we used the [dscr](https://github.com/stephens999/dscr) package to perform dynamic comparisons for different estimation methods on simulated RNA-seq data (simulation based on [GTEx](https://gtexportal.org/home/) RNA-seq dataset). The codes for the simulation studies are in directory `dsc-poisthin-indep` (simulate RNA-seq data with independent genes) and `dsc-poisthin-dep` (simulate RNA-seq data with dependent genes). The comparison results are visualized in `analysis/simulation.Rmd`.

