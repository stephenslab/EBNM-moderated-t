# DSC of RNA-seq differential expression analysis methods 

This directory contains the R code of DSC (Dynamic Statistical Comparison) of DE analysis methods on simulated RNA-Seq data where genes are independent. The simulations are based on real [GTEx](https://gtexportal.org/home/) RNA-seq data of human liver ([raw RNA-seq data](https://github.com/mengyin/EBNM/blob/master/data/Liver.txt)).

To run the DSC, first install the [`dscr`](https://github.com/stephens999/dscr) package (and its required packages dplyr, shiny):

```{r}
install.packages(c("dplyr", "shiny"))

library(devtools)
install_github("stephens999/dscr")
```

Also make sure you already installed the other required R packages: 

* From Bioconductor: limma, DESeq2, edgeR, qvalue, RUVSeq, sva
* From CRAN: ashr, ggplot2, data.table, locfdr, gaussquad, AUC, SQUAREM, foreach, pscl
* From Github: vashr, vicar

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma", "DESeq2", "edgeR", "qvalue", "RUVSeq", "sva"))

install.packages(c("ashr", "ggplot2", "data.table", "locfdr", "gaussquad", "AUC", "SQUAREM", "foreach", "pscl"))

library(devtools)
install_github("mengyin/vashr")
install_github("dcgerard/vicar")
```

Then run the DSC by

```{r}
setwd(".")
source("run_dsc.R")
```

It will generate the simulated RNA-seq datasets for different simulation scenarios, apply different DE methods, and record the final results & scores in the file `res.Rdata`. The whole procedure may take quite some time (at least couple of hours on my machine), and you can reduce the number of scenarios/methods/replicates by editting `methods.R` and `scenarios.R`.

Note: normally I first run the testing script (with restricted number of scenarios and methods) to check if the DSC is runnable:

```{r}
source("run_dsc_test.R")
```

