# A scalable methodology for forecasting short-term high-spatial resolution cancer mortality data

This repository contains the R code to fit the models described in the paper entitled _"A scalable methodology for forecasting short-term high-spatial resolution cancer mortality data"_ (Orozco-Acosta et al., 2023).

## Table of contents

- [Data](#Data)

- [R code](#R-code)

- [References](#References)

# Data

Spanish lung cancer mortality data are used to replicate the analysis of real data in the 7907 municipalities of continental Spain during the period 1991-2015. Due to confidentiality issues with the data provider, simulated dataset is used. To analyse the data, you must have the `bigDM` package installed available [here](https://github.com/spatialstatisticsupna/bigDM).

The [Data_LungCancer](https://github.com/spatialstatisticsupna/bigDM/blob/master/data/Data_LungCancer.rda) file contains the following variables:

- ```ID```: character vector of geographic identifiers
- ```year```: numeric vector of year’s identifiers
- ```obs```: observed number of cases
- ```exp```: expected number of cases
- ```SMR```:  standardized mortality ratios
- ```pob```: population at risk

Use the following command to read the dataset:

```r 
library(bigDM)
data("Data_LungCancer"")
```

# R code

R code to fit the spatio-temporal areal models described in the paper, and to reproduce the results, has been included [here](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/R).

- [fit_model_projection.R](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/R/fit_model_projection.R) is a R script to fit (and recover submodels) the classics and partitioned (disjoint and 1st-order neighbourhood) models.

  - Note: R_INLA version 22.12.16 by default use `"compact"` mode (`"experimental"` mode in previous versions) explained in [Van Niekerk, J. et. al. (2023)](https://www.sciencedirect.com/science/article/pii/S0167947323000038). With `inla.mode` argument of `STCAR_INLA()` function it can be set. For example:
  
  ``` r
  STCAR_INLA(..., inla.mode = "compact")
  ```

- [cv_measures.R](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/R/cv_measures.R) is a R function to compute the cross validation measures implemented in [Liu, Z., & Rue, H. (2022)](https://arxiv.org/pdf/2210.04482.pdf) called Leave-One-Out Cross-Validation (LOOCV) and Leave-Group-Out Cross-Validation (LGOCV) for latent Gaussian models. For example:

```r
  ## loading packages ##
  library(bigDM)
  library(INLA)
  
  ## data preparation ##
  Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID, 1, 2)
  carto <- Carto_SpainMUN[Carto_SpainMUN$ID.prov == "01", ] # 01 Álava
  data <- Data_LungCancer
  data <- data[substr(data$ID, 1, 2) %in% "01", ]
  
  # NA values for predict ##
  data$obs[data$year %in% 2013:2015] <- NA
  
  ## fit INLA model ##
  model <- STCAR_INLA(carto = carto, data = data, ID.area = "ID", ID.year = "year", 
                      O = "obs", E = "pop", spatial = "BYM2", temporal = "rw1", 
                      interaction = "TypeI", model = "global", inla.mode = "compact")
  
  ## compute cross validation measures ##
  source("R/cv_measures.R")
  cvmeasures <- CV(model)
                           
```


# Acknowledgements

This research has been supported by the project PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033. It has also been partially funded by the Public University of Navarre (project PJUPNA20001).

![image](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/micin-aei.jpg)

# References

Liu, Z., & Rue, H. (2022). Leave-Group-Out Cross-Validation For Latent Gaussian Models. _arXiv preprint arXiv:2210.04482_. [https://doi.org/10.48550/arXiv.2210.04482](https://doi.org/10.48550/arXiv.2210.04482)

Orozco-Acosta, E., Riebler, A., Adin, A., and Ugarte, M.D. (2023). A scalable methodology for forecasting short-term high-spatial resolution cancer mortality data. _preprint_.

Van Niekerk, J., Krainski, E., Rustand, D., & Rue, H. (2023). A new avenue for Bayesian inference with INLA. _Computational Statistics & Data Analysis_, 107692. [https://doi.org/10.1016/j.csda.2023.107692](https://doi.org/10.1016/j.csda.2023.107692)
