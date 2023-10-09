# A scalable approach for short-term disease forecasting in high spatial resolution areal data

This repository contains the R code to reproduce the analyses presented in the paper entitled _"A scalable approach for short-term disease forecasting in high spatial resolution areal data"_ [(Orozco-Acosta et al., 2023)](https://arxiv.org/abs/2303.16549). Specifically, it contains several scripts to reproduce the **_Predictive validation study_** and **_Illustration_** sections of the paper. Note that due to confidentiality issues, simulated data sets are provided, and hence, results are not fully reproducible.

Model fitting is performed using the `STCAR_INLA()` function of the R package [**bigDM**](https://cran.r-project.org/web/packages/bigDM/index.html). The package also includes several univariate and multivariate spatial and spatio-temporal Bayesian models for high-dimensional areal count data fitted using the integrated nested Laplace approximation (INLA) estimation technique (http://www.r-inla.org/).

See [https://github.com/spatialstatisticsupna/bigDM](https://github.com/spatialstatisticsupna/bigDM) for details about installation and access to the vignettes accompanying this package.


## Table of contents

- [Data](#Data)
- [R code](#R-code)
- [Acknowledgements](#Acknowledgements)
- [References](#References)


## Data

Our methodology is employed to forecast male lung cancer and overall cancer (all sites) mortality data across the 7,907 municipalities in continental Spain. This forecast extends three years into the future, using the reference period of 1991-2012.

Please, note that data used in this paper are subject to confidentiality agreements with the Spanish Statistical Office (INE), as we analyze cancer mortality data at the municipal level in Spain (NUTS4 level from the European nomenclature of territorial units for statistics).

Then, synthetic data comparable to the original data in size and structure have been included. The data can be directly loaded in R through our package **bigDM** by using the command
```r
> library(bigDM)
> data(Data_LungCancer, package="bigDM")

> head(Data_LungCancer)
     ID year obs        exp      SMR     pop
1 01001 1991   0 0.27554975 0.000000  483.77
2 01002 1991   3 2.72124997 1.102435 4948.98
3 01003 1991   0 0.49612295 0.000000  667.91
4 01004 1991   0 0.37040942 0.000000  591.11
5 01006 1991   0 0.06999827 0.000000   62.81
6 01008 1991   0 0.29240328 0.000000  354.80
```
A copy of this data has also been included in the `data/` folder.

The `Data_LungCancer` object contains the following variables:
- ```ID```: character vector of geographic identifiers
- ```year```: numeric vector containing year identifiers
- ```obs```: observed number of cases
- ```exp```: expected number of cases
- ```SMR```: standardized mortality ratios
- ```pop```: population at risk


Similarly, an `.Rdata` file containing simulated data of male overall cancer mortality counts can be loaded as:
```r
> load("./data/Data_OverallCancer.Rdata")

> head(Data_OverallCancer)
     ID year obs       exp       SMR     pop
1 01001 1991   1 1.0125164 0.9876384  483.77
2 01002 1991  10 9.8199947 1.0183305 4948.98
3 01003 1991   1 1.8756597 0.5331458  667.91
4 01004 1991   1 1.4252320 0.7016401  591.11
5 01006 1991   0 0.2518388 0.0000000   62.81
6 01008 1991   2 1.0323942 1.9372445  354.80
```


## R code

The [R code](./Rcode/) of this repository is organized according to the corresponding sections of the paper.

### Section 4 - Predictive validation study

The script [Section4_PredictiveValidationStudy.R](./Rcode/Section4_PredictiveValidationStudy.R) enables the replication of the predictive validation study presented in Section 4 of Orozco-Acosta et al. (2023) using synthetic data. Similar results to those presented in Table 2, Table A1 and Table A2 will be obtained.

The script is structured in four main steps:

1. **Generate the 8 validation configurations**

    As described in the paper, each configuration uses $T=15$ years of data to fit the model and predict at time points $T+1$, $T+2$ and $T+3$. The first configuration uses data from 1991 to 2005, the second configuration from 1992 to 2006, while the last configuration uses data from 1998 to 2012.
  
2. **Fit the models with INLA using the `STCAR_INLA()` function**

    The script allows to select the model to be fitted using the `model` (one of either `"Classical"`, `"Disjoint"` or `"1st-order nb"`) and `interaction` (one of either `"TypeI"`, `"TypeII"`, `"TypeIII"` or `"TypeIV"`) arguments. Each combination of these arguments will reproduce a line of Table 2.

    Please note that we have set an alternative approximation strategy for INLA (which was not discussed in the paper) as the default option to expedite computations.

3. **Compute model assessment criteria**

    After fitting the models for each configuration, several measures (MAE - mean absolute error, RMSE - root mean square error, IS - interval score) are computed to assess their predictive performance. The results will reproduce similar values to those presented in Table 2, as well as Tables A1 and A2 of the Appendix section. The results are stored in the `results/` folder.
  
    **CAUTION!** These computations might be very time consuming in Windows OS.

4. **Compute Figure 2**
   
    Please note that only the left column graphs of the figure (Disjoint model - Type IV) will be generated. To fully reproduce Figure 2, rerun the script with the arguments `model <- "1st-order nb"` and `interaction <- "TypeIV"` (lines 61-62).
    

### Section 5 - Illustration: projections of cancer mortality in Spain

- The script [Section5_FitModels.R](./Rcode/Section5_FitModels.R) allows to replicate the fit of spatio-temporal *classical* and *partition* models described in the illustration section of Orozco-Acosta et al. (2023) using the [bigDM](https://github.com/spatialstatisticsupna/bigDM) package. The code can be used with any other data sets with similar structure. It also computes the Logarithmic Score (based on both LOOCV and LGOCV approaches) and model selection criteria as DIC and WAIC.

    To reproduce the results shown in **Table 3**, all the models must be previously fitted following the code described in the script. To compute the Logarithmic Score measures for the partition models (i.e., *Disjoint* and *1st-order neighbourhood* models) their corresponding sub-models must be also saved by setting the argument `STCAR_INLA(..., save.models=TRUE)`. Please, note that this function will automatically create a "temp" folder in the current working directory.

    **CAUTION!** These computations are very time consuming (several hours) and generate rather large `.Rdata` files (several Gb) for each of the models.

- The script [Section5_Figures_and_Tables.R](./Rcode/Section5_Figures_and_Tables.R) contains the necessary functions to replicate the figures and tables of Section 5.1 (*Lung cancer mortality*) and Section 5.2 (*Overall cancer mortality*). Note that slightly different results are obtained since we are using simulated counts to preserve the confidentiality of the original data.

    In order to avoid fitting the models (see [Section5_FitModels.R](./Rcode/Section5_FitModels.R) for details), the final INLA models (1st-order neighbourhood + Type IV interaction) can be downloaded from [https://figshare.com/articles/dataset/Scalable_prediction/24229951](https://figshare.com/articles/dataset/Scalable_prediction/24229951).

### Auxiliary functions

The script [Auxiliary_functions.R](./Rcode/Auxiliary_functions.R) contains several additional functions to perform the analyses described in the paper, so that the results are reproducible when using an alternative (but similarly structured) data.


### R and R packages. Version info
``` {r}
R version 4.2.0 (2022-04-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /users/r/R-4.2.0/lib64/R/lib/libRblas.so
LAPACK: /users/r/R-4.2.0/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] ggpubr_0.6.0       ggplot2_3.4.2      tmap_3.3-3         RColorBrewer_1.1-3
 [5] INLA_22.12.16      sp_2.0-0           foreach_1.5.2      Matrix_1.5-4
 [9] Hmisc_5.1-0        bigDM_0.5.1

loaded via a namespace (and not attached):
  [1] leafem_0.2.0       colorspace_2.0-3   ggsignif_0.6.4
  [4] deldir_1.0-6       class_7.3-20       leaflet_2.1.1
  [7] htmlTable_2.4.1    base64enc_0.1-3    dichromat_2.0-0.1
 [10] rstudioapi_0.13    proxy_0.4-26       farver_2.1.0
 [13] listenv_0.8.0      fansi_1.0.3        codetools_0.2-18
 [16] splines_4.2.0      doParallel_1.0.17  knitr_1.39
 [19] rlist_0.4.6.2      Formula_1.2-4      tmaptools_3.1-1
 [22] broom_1.0.5        cluster_2.1.3      png_0.1-8
 [25] compiler_4.2.0     backports_1.4.1    fastmap_1.1.0
 [28] cli_3.6.1          s2_1.1.4           htmltools_0.5.2
 [31] tools_4.2.0        coda_0.19-4        gtable_0.3.0
 [34] glue_1.6.2         dplyr_1.1.2        wk_0.6.0
 [37] gmodels_2.18.1.1   Rcpp_1.0.8.3       carData_3.0-5
 [40] raster_3.5-15      vctrs_0.6.3        spdep_1.2-4
 [43] gdata_2.18.0.1     nlme_3.1-157       iterators_1.0.14
 [46] leafsync_0.1.0     crosstalk_1.2.0    lwgeom_0.2-13
 [49] xfun_0.31          stringr_1.5.0      globals_0.15.0
 [52] rbibutils_2.2.8    lifecycle_1.0.3    gtools_3.9.2.1
 [55] rstatix_0.7.2      XML_3.99-0.10      future_1.26.1
 [58] terra_1.5-21       LearnBayes_2.15.1  MASS_7.3-56
 [61] scales_1.2.0       expm_0.999-6       spatialreg_1.2-3
 [64] gridExtra_2.3      rpart_4.1.16       stringi_1.7.6
 [67] fastDummies_1.6.3  e1071_1.7-9        checkmate_2.2.0
 [70] boot_1.3-28        spData_2.0.1       Rdpack_2.3
 [73] rlang_1.1.1        pkgconfig_2.0.3    evaluate_0.15
 [76] lattice_0.20-45    purrr_1.0.1        sf_1.0-13
 [79] labeling_0.4.2     htmlwidgets_1.5.4  cowplot_1.1.1
 [82] tidyselect_1.2.0   parallelly_1.34.0  magrittr_2.0.3
 [85] R6_2.5.1           generics_0.1.2     DBI_1.1.2
 [88] pillar_1.9.0       foreign_0.8-82     withr_2.5.0
 [91] units_0.8-0        stars_0.6-1        abind_1.4-5
 [94] nnet_7.3-17        tibble_3.2.1       future.apply_1.9.0
 [97] car_3.1-2          crayon_1.5.1       KernSmooth_2.23-20
[100] utf8_1.2.2         rmarkdown_2.14     grid_4.2.0
[103] data.table_1.14.2  digest_0.6.29      classInt_0.4-3
[106] tidyr_1.3.0        munsell_0.5.0      viridisLite_0.4.0
```


# Acknowledgements

This research has been supported by the project PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033.

![](./micin-aei.jpg)


# References

[Orozco-Acosta, E., Riebler, A., Adin, A., and Ugarte, M.D. (2023). A scalable approach for short-term disease forecasting in high spatial resolution areal data. *Biometrical Journal (accepted)*.](https://arxiv.org/abs/2303.16549)
