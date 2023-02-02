# A scalable methodology for forecasting short-term high spatial resolution areal count data

This repository contains the R code to fit the models described in the paper entitled _"A scalable methodology for forecasting short-term high spatial resolution areal count data"_ (Orozco-Acosta et al., 2023).

All the computations are made using the R package [**bigDM**](https://cran.r-project.org/web/packages/bigDM/index.html), which also includes several univariate and multivariate spatial and spatio-temporal Bayesian models for high-dimensional areal count data based on the integrated nested Laplace approximation (INLA) estimation technique. Visit [https://github.com/spatialstatisticsupna/bigDM](https://github.com/spatialstatisticsupna/bigDM) for details about installation and access to the vignettes accompanying this package.


## Table of contents

- [Data](#Lung-cancer-mortality-data)
- [R code](#R-code)
- [References](#References)


# Lung cancer mortality data

We illustrate our proposed methodology by projecting lung cancer mortality data in the 7907 municipalities of continental Spain by considering 3-year ahead predictions usign the period 1991-2012 as reference.

For that purpose we use the `Data_LungCancer` object available at the bigDM package, which contains simulated data of male lung cancer mortality counts (modified in order to preserve the confidentiality of the original data). Specifically, the data set contains the following variables:
- ```ID```: character vector of geographic identifiers
- ```year```: numeric vector of year’s identifiers
- ```obs```: observed number of cases
- ```exp```: expected number of cases
- ```SMR```:  standardized mortality ratios
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


# R code

This section includes the R scripts to fit with [R-INLA](https://www.r-inla.org/) the spatio-temporal projection models described in the present paper.

- [**Fit_model_projection.R**](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/R/Fit_model_projection.R)

  R code to fit the *classical* and *partition* models described in Orozco-Acosta et al. (2023) using the [bigDM](https://github.com/spatialstatisticsupna/bigDM) package.

- [**Figures_and_Tables.R**](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/R/Figures_and_Tables.R)

  R code that contains the necessary functions to replicate most of the figures and tables of the present paper. Note that slightly different results are obtained since we are using simulated counts to preserve the confidentiality of the original data. The final model (1st-order neighbourhood + Type IV interaction) fitted using INLA can be downloaded from [https://emi-sstcdapp.unavarra.es/bigDM/k1_typeIV_simulated.Rdata](https://emi-sstcdapp.unavarra.es/bigDM/k1_typeIV_simulated.Rdata).

- [**CV_measures.R**](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/R/CV_measures.R)

  R code to compute the logarithmic score (sum of the log-predictive densities computed over each area-time point) using both leave-one-out cross-validation (LOOCV) and leave-group-out cross validation (LGOCV) techniques. See [Liu, Z., and Rue, H. (2022)](https://arxiv.org/pdf/2210.04482.pdf) for further details. 
  
  **Note**: In order to compute these measures, the `inla.mode="compact"` argument must be set in the call to the `STCAR_INLA()` function (see the [reference manual](https://cran.r-project.org/web/packages/bigDM/bigDM.pdf) for details).


# Acknowledgements

This research has been supported by the project PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033. It has also been partially funded by the Public University of Navarre (project PJUPNA2001).

![image](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/micin-aei.jpg)


# References

[Liu, Z., and Rue, H. (2022). Leave-Group-Out Cross-Validation For Latent Gaussian Models. _arXiv preprint_.](https://doi.org/10.48550/arXiv.2210.04482)

Orozco-Acosta, E., Riebler, A., Adin, A., and Ugarte, M.D. (2023). A scalable methodology for forecasting short-term high spatial resolution areal count data.
