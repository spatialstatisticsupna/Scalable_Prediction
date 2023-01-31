# Predicting short-term cancer mortality in 7907 municipalities of continental Spain

This repository contains the R code to fit the models described in the paper entitled _"Predicting short-term cancer mortality in 7907 municipalities of continental Spain"_ (Orozco-Acosta et al., 2023).

## Table of contents

- [Installation](#Installation)

- [Data](#Data)

- [R code](#R-code)

- [References](#References)

# Installation

Before replicate the data analysis, it necessary install [bigDM](https://cran.r-project.org/web/packages/bigDM) package. 

## Install from CRAN

``` r
## CRAN version ##
install.packages("bigDM")

```
## Install from Github (development version)

``` r
# Install devtools package from CRAN repository
install.packages("devtools")

# Load devtools library
library(devtools)

# Install the R-INLA package
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# In some Linux OS, it might be necessary to first install the following packages
install.packages(c("cpp11","proxy","progress","tzdb","vroom"))

# Install bigDM from GitHub repositoy
install_github("spatialstatisticsupna/bigDM")
```
**IMPORTANT NOTE: At least the stable version of INLA 22.11.22 (or newest) must be installed for the correct use of the bigDM package.**

# Data

Spanish lung cancer mortality data are used to replicate the analysis of real data in the 7907 municipalities of continental Spain during the period 1991-2015 available [here](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/data). Due to confidentiality issues with the data provider, simulated dataset is used. The [Data_LungCancer.txt](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/data/Data_LungCancer.txt) file contains the following variables:

- ```ID```: character vector of geographic identifiers
- ```year```: numeric vector of year’s identifiers
- ```obs```: observed number of cases
- ```exp```: expected number of cases
- ```SMR```:  standardized mortality ratios
- ```pob```: population per area-time unit

Use the following command to read the dataset:

```r 
data <- read.table(file = "data/Data_LungCancer.txt", header = TRUE, 
                   colClasses = c("character",rep("numeric",5)))
```

# R code

R code to fit the spatio-temporal areal models described in the paper, and to reproduce the results, has been included [here](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/R).

- [fit_model_projection.R](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/R/fit_model_projection.R) is a R script to fit (and recover submodels) the classics and partitioned (disjoint and 1st-order neighbourhood) models.

  - Note: R_INLA version 22.12.16 by default use "compact" mode (in previous versions was "experimental" mode) explained in [Van Niekerk, J. et. al. (2023)](https://www.sciencedirect.com/science/article/pii/S0167947323000038). For change to classical mode, use the following command:
  
  ``` r
  inla.setOption(inla.mode="classic")
  ```

- [cv_measures.R](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/R/cv_measures.R) is a R function to compute the cross validation measures implemented in [Liu, Z., & Rue, H. (2022)](https://arxiv.org/pdf/2210.04482.pdf) called Leave-One-Out Cross-Validation (LOOCV) and Leave-Group-Out Cross-Validation (LGOCV) for latent Gaussian models. For example:

```r
## loading packages ##
library(INLA)
library(bigDM)

## reading data and cartography files ##
data <- read.table(file = "data/Data_LungCancer.txt", header = TRUE, 
                   colClasses = c("character",rep("numeric",5)))
Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID, 1, 2)

## NA for values to be predicted ##
pred <- 2013:2015
data$obs[data$year %in% pred] <- NA

## filtering for álava province ##
data <- data[substr(data$ID,1,2) %in% "01",]
carto <- Carto_SpainMUN[Carto_SpainMUN$ID.prov=="01",]

## fit INLA model ##
model <- STCAR_INLA(carto = carto, data = data, ID.area = "ID", 
                    ID.year = "year", O = "obs", E = "pob", spatial = "BYM2",
                    temporal = "rw1", interaction = "TypeI", model = "global")

## compute cross validation measures ##
source("R/cv_measures.R")
cvmeasures <- CV(model)
                           
```


# Acknowledgements

This research has been supported by the project PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033. It has also been partially funded by the Public University of Navarre (project PJUPNA20001).

![image](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/micin-aei.jpg)

# References

Liu, Z., & Rue, H. (2022). Leave-Group-Out Cross-Validation For Latent Gaussian Models. _arXiv preprint arXiv:2210.04482_. [https://doi.org/10.48550/arXiv.2210.04482](https://doi.org/10.48550/arXiv.2210.04482)

Orozco-Acosta, E., Riebler, A., Adin, A., and Ugarte, M.D. (2023). Predicting short-term cancer mortality in 7907 municipalities of continental Spain. _preprint_.

Van Niekerk, J., Krainski, E., Rustand, D., & Rue, H. (2023). A new avenue for Bayesian inference with INLA. _Computational Statistics & Data Analysis_, 107692. [https://doi.org/10.1016/j.csda.2023.107692](https://doi.org/10.1016/j.csda.2023.107692)
