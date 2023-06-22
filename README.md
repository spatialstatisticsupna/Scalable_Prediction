# A scalable approach for short-term disease forecasting in high spatial resolution areal data

This repository contains the R code to fit the models described in the paper entitled _"A scalable approach for short-term disease forecasting in high spatial resolution areal data"_ [(Orozco-Acosta et al., 2023)](https://arxiv.org/abs/2303.16549), as well as to reproduce similar figures and tables to the ones presented in the paper. Note that due to confidentiality issues with the real data, a simulated dataset has been used and then, results are not fully reproducible.

All computations were performed using version 0.5.1 of the R package [**bigDM**](https://cran.r-project.org/web/packages/bigDM/index.html), which was released on CRAN on February 23, 2023. This version includes adaptations to the `STCAR_INLA()`  function, making it compatible for short-term forecasting.
The package also includes several univariate and multivariate spatial and spatio-temporal Bayesian models for high-dimensional areal count data based on the integrated nested Laplace approximation (INLA) estimation technique.

See [https://github.com/spatialstatisticsupna/bigDM](https://github.com/spatialstatisticsupna/bigDM) for details about installation and access to the vignettes accompanying this package.


## Table of contents

- [R code](#R-code)
- [Acknowledgements](#Acknowledgements)
- [References](#References)


# R code

R code to fit the models and reproduce similar figures and tables to the ones presented in the paper is included [here](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/master/R/)


# Acknowledgements

This research has been supported by the project PID2020-113125RBI00/MCIN/AEI/10.13039/501100011033.

![image](https://github.com/spatialstatisticsupna/Scalable_Prediction/blob/main/micin-aei.jpg)


# References

[Liu, Z., and Rue, H. (2022). Leave-Group-Out Cross-Validation For Latent Gaussian Models. _arXiv preprint_.](https://doi.org/10.48550/arXiv.2210.04482)

[Orozco-Acosta, E., Riebler, A., Adin, A., and Ugarte, M.D. (2023). A scalable approach for short-term disease forecasting in high spatial resolution areal data. _arXiv preprint_.](https://arxiv.org/abs/2303.16549)
