
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coNBMTX

<!-- badges: start -->
<!-- badges: end -->

The goal of coNBMTX is to implement Differential Expression (DE) for
metatranscriptomics with sample-paired metagenomics data.

## Installation

To install coNBMTX, you need to firstly install edgeR and metagenomeSeq:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("metagenomeSeq")
```

You can install the development version of coNBMTX like so:

``` r
library(devtools)
install_github('Lizz647/coNBMTX')
```

## Example

This is a basic example which shows you how to use coNBMTX:

``` r
library(coNBMTX)

# Obtain the sample-paired data, with genes on rows and samples on columns.
mtxtable <- matrix(rnbinom(5000,size=10,mu=50),nrow=100)
mgxtable <- matrix(rnbinom(5000,size=5,mu=30),nrow=100)
metatable <- data.frame(ID=paste("S",seq(1:50),sep=""),group=c(rep("a",25),rep("b",25)))

# Packaging data using a list.
datamat = list(mtx=mtxtable,mgx=mgxtable,metadata=metatable)

# Specify the formula in regression.
group_formula = "~group"

output = coNBMTX(datamat=datamat,group_formula=group_formula)
#> Preprocessing...
#> Caculating normalize factors...
#> Default value being used.
#> Default value being used.
#> Estimating genewise-dispersion...
#> Estimating non-conditional genewise-dispersion...
#> Estimating conditional genewise-dispersion...
#> Estimating posterior dispersion...
#> Estimating fold changes...
#> Inference...
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
