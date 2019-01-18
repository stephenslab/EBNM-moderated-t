# Empirical Bayes Normal Means with Estimated Standard errors

This repository contains the codes and analyses for 
work by Lu and Stephens, Empirical Bayes (EB) estimation of normal means, accounting for uncertainty in estimated standard errors.

## The bottom line, and example pipeline

The bottom line is that we recommend a two-stage procedure:

  1. Perform EB shrinkage on estimated standard errors; eg as in the limma R package.

  2. Apply the EB T means model (ie the EB normal means, but with the normal likelihood replaced with a t likelihood); eg as in the ashr R package.
 
An example of this pipeline, which we call *VL+eBayes+ash*, is illustrated in `analysis/method.Rmd`.

## Simulation studies

For simulation studies, we used the [dscr](https://github.com/stephens999/dscr) package to perform dynamic comparisons for different estimation methods on simulated RNA-seq data (simulation based on [GTEx](https://gtexportal.org/home/) RNA-seq dataset). The codes for the simulation studies are in directory `dsc-poisthin-indep` (simulate RNA-seq data with independent genes) and `dsc-poisthin-dep` (simulate RNA-seq data with dependent genes). 

To run the simulation and plot the results see the code in `analysis/simulation.Rmd`.
(There is also code to produce `.eps` files of the figures in `analysis/plots.R`).
