## This R-package constructs valid inference for penalized G-estimation under UPoSI approach.

This repository contains the following folders:

Folder | Description
--- | ---
R | Contains the source codes
man | Contains the documentation of each function used

The R folder contains the following files:

File | Description
--- | ---
[uposi_confint.R](https://github.com/ajmeryjaman/UPoSIPeG/blob/main/R/uposi_confint.R) | Contains the main function uposi_confint() which implements our method
[bivar_quantile.R](https://github.com/ajmeryjaman/UPoSIPeG/blob/main/R/bivar_quantile.R) | Contains the function bivar_quantile() that calculates the required bivariate quantiles used to construct UPoSI confidence intervals

Please see the example given in [uposi_confint.R](https://github.com/ajmeryjaman/UPoSIPeG/blob/main/R/uposi_confint.R) to generate a longitudinal data set and implement our method. Or, do the following:

#### R commands for installing and using our package

library(devtools) # if already installed, otherwise need to install it first

install_github("ajmeryjaman/UPoSIPeg")

library(UPoSIPeg)

?uposi_confint # To follow the example
