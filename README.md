# FMRHH

## Overview
FMRHH is an R package designed to identify homogeneous and heterogeneous variables in a mixture of Cox regression model and linear regression model. It uses EM and ADMM algorithms to fit the models and provide insights into the data.

## Installation
To install FMRHH, use the following R command:

```R
devtools::install_github("zongminliu9/FMRHH")

## Usage
To get started, load the package in R:
library(FMRHH)

## Functions
admm_cox: Iterates M-step in the EM algorithm for Cox proportional hazard model using the ADMM algorithm.R
admm_linear2: Iterates M-step in the EM algorithm for linear regression model with two types of heterogeneity using ADMM algorithm.
admm_linear3: Iterates M-step in the EM algorithm for linear regression model with three types of heterogeneity using ADMM algorithm.
coxHH: Fits a mixture of Cox regression models using EM and ADMM algorithm. Includes functions update_theta_cox and update_phi_cox for parameter updates.
linearHH: Fits a mixture of linear regression models using EM and ADMM algorithm. Includes functions update_theta_cox, gd_theta, gd_phi, and glasso for parameter updates.
update_admm_cox: Updates iteration parameters and returns coefficients estimates of beta and alpha for Cox regression iteration.
update_admm_linear2: Updates iteration parameters and returns coefficients estimates of beta and alpha for mixture linear regression model with two types of heterogeneity.
update_admm_linear3: Updates iteration parameters and returns coefficients estimates of beta and alpha for mixture linear regression model with three types of heterogeneity.
