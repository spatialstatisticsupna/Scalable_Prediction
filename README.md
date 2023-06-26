# A scalable approach for short-term disease forecasting in high spatial resolution areal data

This repository contains the R code to fit the models described in the paper entitled _"A scalable approach for short-term disease forecasting in high spatial resolution areal data"_ [(Orozco-Acosta et al., 2023)](https://arxiv.org/abs/2303.16549), as well as to reproduce similar figures and tables to the ones presented in the paper. Note that due to confidentiality issues, simulated data sets has are provided and then, results are not fully reproducible.

All computations were performed using version 0.5.1 of the R package [**bigDM**](https://cran.r-project.org/web/packages/bigDM/index.html), which was released on CRAN on February 23, 2023. This version includes adaptations to the `STCAR_INLA()` function, making it compatible for short-term forecasting. 
The package also includes several univariate and multivariate spatial and spatio-temporal Bayesian models for high-dimensional areal count data based on the integrated nested Laplace approximation (INLA) estimation technique (http://www.r-inla.org/).

See [https://github.com/spatialstatisticsupna/bigDM](https://github.com/spatialstatisticsupna/bigDM) for details about installation and access to the vignettes accompanying this package.


## Table of contents

- [Cancer mortality data](#Cancer-mortality-data)
- [R code](#R-code)
- [Acknowledgements](#Acknowledgements)
- [References](#References)


## Cancer mortality data

The proposed methodology is applied to project male lung cancer and overall cancer (all sites) mortality data in the 7907 municipalities of continental Spain by considering three-year ahead predictions using the period 1991-2012 as reference.

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
     ID year obs       exp       SMR     pop
1 01001 1991   1 1.0125164 0.9876384  483.77
2 01002 1991  10 9.8199947 1.0183305 4948.98
3 01003 1991   1 1.8756597 0.5331458  667.91
4 01004 1991   1 1.4252320 0.7016401  591.11
5 01006 1991   0 0.2518388 0.0000000   62.81
6 01008 1991   2 1.0323942 1.9372445  354.80
```


## R code

The code of this paper is organized in self-contained folders, which are named according to the corresponding sections of the paper.


[**Section 4. Predictive validation study**](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/Section4_PredictiveValidationStudy)

- **IMPORTANT NOTE**: To reproduce Table 2, it is necessary to fit all the models for each configuration set of the validation study, which entails performing computations for over 18 days using the computational architecture outlined in Section 4.2

- The following code can be used to replicate the definition of data configurations for the validation setup described in Section 4.
```r
> t.min <- min(Data_LungCancer$year)
> t.max <- max(Data_LungCancer$year)
> 
> t.length <- 15
> t.periods <- 3
> 
> Data_pred <- lapply(1:8, function(x){
+   cc <- seq(t.min+(x-1), t.min+t.length+(x+1))
+   data <- Data_LungCancer[Data_LungCancer$year %in% cc,]
+   data$obs[data$year %in% tail(unique(data$year), n=t.periods)] <- NA
+   
+   return(data)
+ })
> 
> names(Data_pred) <- paste("config",1:8,sep=".")
> str(Data_pred,2)
List of 8
 $ config.1:'data.frame':	142326 obs. of  6 variables:
  ..$ ID  : chr [1:142326] "01001" "01002" "01003" "01004" ...
  ..$ year: int [1:142326] 1991 1991 1991 1991 1991 1991 1991 1991 1991 1991 ...
  ..$ obs : int [1:142326] 0 3 0 0 0 0 1 2 0 0 ...
  ..$ exp : num [1:142326] 0.276 2.721 0.496 0.37 0.07 ...
  ..$ SMR : num [1:142326] 0 1.1 0 0 0 ...
  ..$ pop : num [1:142326] 483.8 4949 667.9 591.1 62.8 ...
 $ config.2:'data.frame':	142326 obs. of  6 variables:
  ..$ ID  : chr [1:142326] "01001" "01002" "01003" "01004" ...
  ..$ year: int [1:142326] 1992 1992 1992 1992 1992 1992 1992 1992 1992 1992 ...
  ..$ obs : int [1:142326] 0 2 0 0 0 0 0 1 1 1 ...
  ..$ exp : num [1:142326] 0.2865 2.7962 0.5091 0.3791 0.0709 ...
  ..$ SMR : num [1:142326] 0 0.715 0 0 0 ...
  ..$ pop : num [1:142326] 500.9 4909 673.2 593.5 63.2 ...
...
```  
  
- [**Validation_Study.R**](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/R/Section4_PredictiveValidationStudy/Validation_Study.R)

  R code to reproduce **Figure 2**: *One, two and three-year ahead predictions for the municipalities of Madrid, Palencia and Ávila using the disjoint model (left column) and 1st-order neighbourhood model (right column) with Type IV interactions.*
  
&nbsp;

[**Section 5. Illustration: projections of cancer mortality in Spain**](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/R/Section5_Illustration)

- [**Fit_models.R**](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/R/Section5_Illustration/Fit_models.R)

  This R script contains the necessary functions to replicate the fit of spatio-temporal *classical* and *partition* models considered in the illustration section of Orozco-Acosta et al. (2023) using the [bigDM](https://github.com/spatialstatisticsupna/bigDM) package. The code can be used with any other data sets with similar structure.
  
  It also computes the Logarithmic Score (based on both LOOCV and LGOCV approaches) and model selection criteria (DIC and WAIC) to reproduce the results shown in **Table 3**.
  
  **IMPORTANT NOTE**: In order to compute Logarithmic Score measures under the partition models (Disjoint or k-order neighbourhood models) the sub-models must be previously saved by setting the argument `STCAR_INLA(..., save.models=TRUE)`.


- [**LungCancer_Results.R**](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/R/Section5_Illustration/LungCancer_Results.R)

  R code that contains the necessary functions to replicate the figures and tables of Section 5.1 (*Lung cancer mortality*). Note that slightly different results are obtained since we are using simulated counts to preserve the confidentiality of the original data. The final model (1st-order neighbourhood + Type IV interaction) fitted using INLA can be downloaded from [https://emi-sstcdapp.unavarra.es/bigDM/inst/Rdata/INLAmodel_LungCancer.Rdata](https://emi-sstcdapp.unavarra.es/bigDM/inst/Rdata/INLAmodel_LungCancer.Rdata).
  
  - **Figure 3a:** Maps of posterior median estimates of lung cancer mortality rates per 100,000 males
  - **Figure 4a:** Posterior predictive median estimates of lung cancer mortality rates and its corresponding 95% credible intervals per 100,000 males for the municipalities of Gerona, Madrid and Bilbao.
  - **Table 4:** Posterior median estimates of the predicted lung cancer mortality rates per 100,000 males, its corresponding 95% credible intervals (CI) and width of the CIS for years 2013 and 2015 for the 47 municipalities that form the provincial capitals.
  

- [**OverallCancer_Results.R**](https://github.com/spatialstatisticsupna/Scalable_Prediction/tree/main/R/Section5_Illustration/OverallCancer_Results.R)

  R code that contains the necessary functions to replicate the figures and tables of Section 5.2 (*Overall cancer mortality*). Note that slightly different results are obtained since we are using simulated counts to preserve the confidentiality of the original data. The final model (1st-order neighbourhood + Type IV interaction) fitted using INLA can be downloaded from [https://emi-sstcdapp.unavarra.es/bigDM/inst/Rdata/INLAmodel_OverallCancer.Rdata](https://emi-sstcdapp.unavarra.es/bigDM/inst/Rdata/INLAmodel_OverallCancer.Rdata).
  
  - **Figure 3b:** Posterior median estimates of overall cancer mortality rates per 100,000 males
  - **Figure 4a:** Posterior predictive median estimates of overall cancer mortality rates and its corresponding 95% credible intervals per 100,000 males for the municipalities of Gerona, Madrid and Bilbao.
  - **Table 5:** Posterior median estimates of the predicted overall cancer mortality rates per 100,000 males, its corresponding 95% credible intervals (CI) and width of the CIS for years 2013 and 2015 for the 47 municipalities that form the provincial capitals.


### Version info
``` {.r}

```


# Acknowledgements

This research has been supported by the project PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033.

![image](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/micin-aei.jpg)


# References

[Liu, Z., and Rue, H. (2022). Leave-Group-Out Cross-Validation For Latent Gaussian Models. _arXiv preprint_.](https://doi.org/10.48550/arXiv.2210.04482)

[Orozco-Acosta, E., Riebler, A., Adin, A., and Ugarte, M.D. (2023). A scalable approach for short-term disease forecasting in high spatial resolution areal data. _arXiv preprint_.](https://arxiv.org/abs/2303.16549)
