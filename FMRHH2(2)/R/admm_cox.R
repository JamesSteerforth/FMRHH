#' ADMM algorithm for a 3-component mixture Cox regression model
#'
#' Update the parameter estimates using ADMM algorithm in EM algorithm
#' @param t observed survival time.
#' @param status observed censoring indicator.
#' @param x covariates.
#' @param index component index.
#' @param beta_0 a vector of the average effect beta.
#' @param alpha1_0 a vector of the deviance of effects in the 1st component from the average effect.
#' @param alpha2_0 a vector of the deviance of effects in the 2nd component from the average effect.
#' @param alpha3_0 a vector of the deviance of effects in the 3rd component from the average effect.
#' @param vlambda1 a vector of tuning parameters for beta.
#' @param vlambda2 a vector of tuning parameters for alpha.
#' @param eps the convergence criteria.
#' @param max maximum number of iterations.
#'
#' @return coefficients estimates of beta and alpha
#' @import survival
#' @export


admm_cox <- function(t,status,x,index,beta_0,alpha1_0,alpha2_0,alpha3_0,vlambda1,vlambda2,eps,max){
  P <- length(beta)

  beta <- beta_0
  alpha1 <- alpha1_0
  alpha2 <- alpha2_0
  alpha3 <- alpha3_0

  beta_t <- beta_0
  alpha1_t <- alpha1_0
  alpha2_t <- alpha2_0
  alpha3_t <- alpha3_0

  phi <- rep(0,5*P)

  convergence <- 1000
  i <- 1
  while (convergence>eps & i<max){
    new.theta <- update_theta(t,status,x,index,beta,alpha1,alpha2,alpha3,phi,beta_t,alpha1_t,alpha2_t,alpha3_t)
    new.beta <- new.theta$beta
    new.alpha1 <- new.theta$alpha1
    new.alpha2 <- new.theta$alpha2
    new.alpha3 <- new.theta$alpha3

    new.psi <- update_thetat(new.beta,new.alpha1,new.alpha2,new.alpha3,phi,vlambda1,vlambda2)
    new.beta_t <- new.psi$beta_t
    new.alpha1_t <- new.psi$alpha1_t
    new.alpha2_t <- new.psi$alpha2_t
    new.alpha3_t <- new.psi$alpha3_t

    new.phi <- update_phi(phi,new.beta,new.alpha1,new.alpha2,new.alpha3,new.beta_t,new.alpha1_t,new.alpha2_t,new.alpha3_t)

    convergence <- crossprod(beta-new.beta)+crossprod(alpha1-new.alpha1)+crossprod(alpha2-new.alpha2)+crossprod(alpha3-new.alpha3)+
      crossprod(beta_t-new.beta_t)+crossprod(alpha1_t-new.alpha1_t)+crossprod(alpha2_t-new.alpha2_t)+crossprod(alpha3_t-new.alpha3_t)+
      crossprod(phi-new.phi)

    beta <- new.beta
    alpha1 <- new.alpha1
    alpha2 <- new.alpha2
    alpha3 <- new.alpha3
    beta_t <- new.beta_t
    alpha1_t <- new.alpha1_t
    alpha2_t <- new.alpha2_t
    alpha3_t <- new.alpha3_t
    phi <- new.phi

    i <- i+1
    convergence
    beta_t
    alpha1_t
    alpha2_t
    alpha3_t
  }
  return(list(beta=beta_t,alpha1=alpha1_t,alpha2=alpha2_t,alpha3=alpha3_t))
}


