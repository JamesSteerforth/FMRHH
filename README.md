# FMRHH

## Overview
FMRHH is an R package designed to identify homogeneous and heterogeneous variables in a mixture of Cox regression model and linear regression model. It uses EM and ADMM algorithms to fit the models and provide insights into the data.

## Installation
To install FMRHH, use the following R command:

## Usage
To get started, load the package in R:
library(FMRHH)

## Functions
# Iterates M-step in the EM algorithm for Cox proportional hazard model using the ADMM algorithm.
admm_cox <- admm_cox(t, status, x, index, beta_0, alpha1_0, alpha2_0, alpha3_0, vlambda1, vlambda2, eps, max)

# Iterates M-step in the EM algorithm for linear regression model with two types of heterogeneity using ADMM algorithm.
admm_linear2 <- admm_linear2(y,x,w,beta_0,alpha1_0,alpha2_0,vlambda1,vlambda2,eps=1e-7,max=50)

# Iterates M-step in the EM algorithm for linear regression model with three types of heterogeneity using ADMM algorithm.
admm_linear3 <- admm_linear3(y,x,w,beta_0,alpha1_0,alpha2_0,alpha3_0,vlambda1,vlambda2,eps=1e-6,max=50)

# Fits a mixture of Cox regression models using EM and ADMM algorithm. Includes functions update_theta_cox and update_phi_cox for parameter updates.
update_theta_cox <- update_theta_cox(t,status,x,index,beta,alpha1,alpha2,alpha3,phi,beta_t,alpha1_t,alpha2_t,alpha3_t)
update_thetat_cox <- update_thetat_cox(beta,alpha1,alpha2,alpha3,phi,vlambda1,vlambda2)
update_phi_cox <- update_phi_cox(phi_0,beta,alpha1,alpha2,alpha3,beta_t,alpha1_t,alpha2_t,alpha3_t)
coxHH <- coxHH(formula,data,lambda,max=20,eps=1e-6)

# Fits a mixture of linear regression models using EM and ADMM algorithm. Includes functions update_theta_cox, gd_theta, gd_phi, and glasso for parameter updates.
gd_theta <- gd_theta(y,x,w,beta,alpha1,alpha2,phi,beta_t,alpha1_t,alpha2_t)
glasso <- glasso(beta,alpha1,alpha2,phi,vlambda1,vlambda2)
gd_phi <- gd_phi(phi_0,beta,alpha1,alpha2,beta_t,alpha1_t,alpha2_t)
admm <- admm(y,x,w,beta_0,alpha1_0,alpha2_0,vlambda1,vlambda2,eps=1e-7,max=50)
linearHH <- linearHH(formula,data,lambda,emmax=50,eps=1e-6)
tuning <- tuning(vlambda1=1:5,vlambda2=1:5,formula,data)

# Updates iteration parameters and returns coefficients estimates of beta and alpha for Cox regression iteration.
update_theta_cox <- update_admm_cox(t,status,x,index,beta,alpha1,alpha2,alpha3,phi,beta_t,alpha1_t,alpha2_t,alpha3_t)

# Updates iteration parameters and returns coefficients estimates of beta and alpha for mixture linear regression model with two types of heterogeneity.
gd_theta2 <- gd_theta2(y,x,w,beta,alpha1,alpha2,phi,beta_t,alpha1_t,alpha2_t)
glasso2 <- glasso2(beta,alpha1,alpha2,phi,vlambda1,vlambda2)
gd_phi2 <- gd_phi2(phi_0,beta,alpha1,alpha2,beta_t,alpha1_t,alpha2_t)

# Updates iteration parameters and returns coefficients estimates of beta and alpha for mixture linear regression model with three types of heterogeneity.
gd_theta3 <- gd_theta3(y,x,w,beta,alpha1,alpha2,alpha3,phi,beta_t,alpha1_t,alpha2_t,alpha3_t)
glasso3 <- glasso3(beta,alpha1,alpha2,alpha3,phi,vlambda1,vlambda2)
gd_phi3 <- gd_phi3(phi_0,beta,alpha1,alpha2,alpha3,beta_t,alpha1_t,alpha2_t,alpha3_t)

