# Code to estimate causal effects for COPD Emergency Department visits
This repository contains code for a Bayesian spatial causal inference model. Our model estimates the causal effects of county-level exposures on the number of chronic obstructive pulmonary disease (COPD) emergency department (ED) visits. The dataset comes from the counties in North Carolina in 2023. 

## Packages Required

In order to get this code to run, you must have a version of R >= 3.5.0 installed. It would be ideal to have the most recent version of R. You must install the R packages that appear in the lines `library(...)` to run these R scripts. Most of these can be installed using the R command line or the RStudio console or GUI. You also need to install `rstan`. To install `rstan`, you must first configure your R installation to be able to compile C++ code. See [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for instructions. Once the C++ toolchain is configured, you can then install `rstan` using 

```
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
```

## Income as exposure
This folder contains code to implement the model with median income (log-transformed) as the exposure. All files in this folder must be saved in the SAME local directory for the code to run. The main scripts are as follows:
* `Car_Stage_1_income.R`: Implements the Stage 1 generalized propensity score (GPS) model for median income
* `Car_Stage_2_income.R`: Implements the Stage 2 outcome regression model for estimating the causal relative risk of median income on the number of COPD ED visits
* `get_exposure_response_curve_income.R`: This estimates the causal exposure-response curve with the average causal effect of median income on the number of COPD ED visits
* `GPS_CAR_model.stan`: MCMC sampler for fitting the Stage 1 GPS model
* `outcome_CAR_mode.stan`: MCMC sampler for fitting the Stage 2 outcome regression 

## Vape shop density as exposure
This folder contains code to implement the model with vape shop density (zero-truncated and log-transformed) as the exposure. All files in this folder must be saved in the SAME local directory for the code to run. The main scripts are as follows:
* `Car_Stage_1_vape_shops.R`: Implements the Stage 1 generalized propensity score (GPS) model for vape shop density
* `Car_Stage_2_vape_shops.R`: Implements the Stage 2 outcome regression model for estimating the causal relative risk of vape shpo density on the number of COPD ED visits
* `get_exposure_response_curve_vape_shops.R`: This estimates the causal exposure-response curve with the average causal effect of vape shop density on the number of COPD ED visits
* `GPS_CAR_model.stan`: MCMC sampler for fitting the Stage 1 GPS model
* `outcome_CAR_mode.stan`: MCMC sampler for fitting the Stage 2 outcome regression model

## Visualization and test for spatial autocorrelation
This folder contains code to make heatmaps of the COPD ED visit rate, median income (log-transformed), and vape shop density (zero-truncated and log-transformed) for all 100 North Carolina counties. The code also conducts the global Moran's I test for positive spatial autocorrelation for each of these variables. The main script is:
* `heatmaps_and_Morans_I.R`: Makes the heatmaps and conducts the global Moran's I tests
