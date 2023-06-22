# R CODE

This folder contains the R code to fit the models and reproduce similar figures and tables to the ones described in the paper _"A scalable approach for short-term disease forecasting in high spatial resolution areal data"_ [(Orozco-Acosta et al., 2023)](https://arxiv.org/abs/2303.16549). Due to confidentiality issues simulated data sets are provided, so is not expected to get the same results given in the paper.

All computations were performed using version 0.5.1 of the R package [**bigDM**](https://cran.r-project.org/web/packages/bigDM/index.html), which was released on CRAN on February 23, 2023. This version includes adaptations to the `STCAR_INLA()`  function, making it compatible for short-term forecasting. The package also includes several univariate and multivariate spatial and spatio-temporal Bayesian models for high-dimensional areal count data based on the integrated nested Laplace approximation (INLA) estimation technique.

See [https://github.com/spatialstatisticsupna/bigDM](https://github.com/spatialstatisticsupna/bigDM) for details about installation and access to the vignettes accompanying this package.

The code of this paper is organized in self-contained folders, which are named according to the corresponding sections of the paper.


## Table of contents

- [Data](#Lung-cancer-mortality-data)
- [Section 4. Predictive validation study](#Section-4.-Predictive-validation-study)
- [Section 5. Illustration: projections of cancer mortality in Spain](#Section-5.-Illustration:-projections-of-cancer-mortality-in-Spain-study)
- [Version info](#Version-info)


## Cancer mortality data

The proposed methodology is applied to project male lung cancer and overall cancer (all cites) mortality data in the 7907 municipalities of continental Spain by considering three-year ahead predictions using the period 1991-2012 as reference.

The `Data_LungCancer` object available at the **bigDM** package, contains simulated data of male lung cancer mortality counts (modified in order to preserve the confidentiality of the original data). Specifically, the data set contains the following variables:
- ```ID```: character vector of geographic identifiers
- ```year```: numeric vector of years identifiers
- ```obs```: observed number of cases
- ```exp```: expected number of cases
- ```SMR```: standardized mortality ratios
- ```pop```: population at risk

Use the following commands to load the data
```r 
> library(bigDM)
> data("Data_LungCancer")

> head(Data_LungCancer)
     ID year obs        exp      SMR     pop
1 01001 1991   0 0.27554975 0.000000  483.77
2 01002 1991   3 2.72124997 1.102435 4948.98
3 01003 1991   0 0.49612295 0.000000  667.91
4 01004 1991   0 0.37040942 0.000000  591.11
5 01006 1991   0 0.06999827 0.000000   62.81
6 01008 1991   0 0.29240328 0.000000  354.80
```

Similarly, an Rdata containing simulated data of male overall cancer mortality counts can be loaded as:
```r 
> load(url("https://emi-sstcdapp.unavarra.es/bigDM/inst/Rdata/Data_OverallCancer.Rdata"))

> head(Data_OverallCancer)
     ID year obs       exp       SMR        pop
1 01001 1991   1 1.0125164 0.9876384  483.77402
2 01002 1991  10 9.8199947 1.0183305 4948.98255
3 01003 1991   1 1.8756597 0.5331458  667.90829
4 01004 1991   1 1.4252320 0.7016401  591.11323
5 01006 1991   0 0.2518388 0.0000000   62.81084
6 01008 1991   2 1.0323942 1.9372445  354.79780
```


## [Section 4. Predictive validation study](https://github.com/spatialstatisticsupna/Scable_Prediction/tree/main/R/Section4_PredictiveValidationStudy)


## [Section 5. Illustration: projections of cancer mortality in Spain](https://github.com/spatialstatisticsupna/Scable_Prediction/tree/main/R/Section5_Illustration)


## Version info

``` {.r}
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17763)

  
```