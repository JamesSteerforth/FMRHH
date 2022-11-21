#' Identify of homogeneous and heterogeneous variables in a mixture of linear regression model
#'
#' Fit a mixture of linear regression model using EM and ADMM algorithm
#' @param formula a formula object, with the response in the left of a ~, and the terms on the right. The response must be a survival object as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the formula.
#' @param K the number of components, K=2 or 3.
#' @param lambda specify the tuning parameters. The first number is the tuning for average effect beta, and the second number is the tuning for deviance effect alpha.
#' @param eps the convergence criteria.
#' @param max maximum number of iterations.
#'
#' @return coefficients estimates of beta and alpha
#' @import flexmix
#' @export
#' Internal function for admm_linear2
#' @description an internal function to update theta in admm through gradient descent
#' @keywords internal
#' ADMM algorithm for a 2-component mixture linear regression model
#'
#' Update the parameter estimates using ADMM algorithm in EM algorithm
#' @param y observed response variable.
#' @param x covariates.
#' @param w weights.
#' @param index component index.
#' @param beta_0 a vector of the average effect beta.
#' @param alpha1_0 a vector of the deviance of effects in the 1st component from the average effect.
#' @param alpha2_0 a vector of the deviance of effects in the 2nd component from the average effect.
#' @param vlambda1 a vector of tuning parameters for beta.
#' @param vlambda2 a vector of tuning parameters for alpha.
#' @param eps the convergence criteria.
#' @param max maximum number of iterations.
#'
#' @return coefficients estimates of beta and alpha
#' @export
gd_theta <- function(y,x,w,beta,alpha1,alpha2,phi,beta_t,alpha1_t,alpha2_t){
  P <- length(beta)
  f <- function(theta){
    temp1 <- (y-theta[1:P]%*%t(x)-theta[(P+1):(2*P)]%*%t(x))^2*w[1,]/2+(y-theta[1:P]%*%t(x)-theta[(2*P+1):(3*P)]%*%t(x))^2*w[2,]/2
    temp2 <- phi[1:P]*apply(rbind(theta[(P+1):(2*P)],theta[(2*P+1):(3*P)]),2,sum)+phi[(P+1):(2*P)]*theta[1:P]+phi[2*P+2*(1:P)-1]*theta[(P+1):(2*P)]+phi[2*P+2*(1:P)]*theta[(2*P+1):(3*P)]
    temp3 <- (apply(rbind(theta[(P+1):(2*P)],theta[(2*P+1):(3*P)]),2,sum))^2+(theta-theta_t)^2
    sum(temp1)+sum(temp2)+1/2*sum(temp3)
  }

  theta_t <- c(beta_t,alpha1_t,alpha2_t)
  out <- optim(c(beta,alpha1,alpha2),f,method="BFGS")$par
  return(list(beta=out[1:P],alpha1=out[(P+1):(2*P)],alpha2=out[(2*P+1):(3*P)]))
}

glasso <- function(beta,alpha1,alpha2,phi,vlambda1,vlambda2){
  P <- length(beta)
  beta_t <- rep(0,P)
  alpha1_t <- rep(0,P)
  alpha2_t <- rep(0,P)
  # update beta_t
  for(g in 1:P){
    f <- function(x){0.5*(x-beta[g]-phi[P+g])^2+vlambda1[g]*abs(x)}
    beta_t[g] <- round(optim(0,f,method="BFGS")$par,4)
  }

  # update alpha_t
  for(g in 1:P){
    f <- function(x){0.5*(x[1]-alpha1[g]-phi[2*P+2*g-1])^2+0.5*(x[2]-alpha2[g]-phi[2*P+2*g])^2+vlambda2[g]*sqrt(crossprod(c(x[1],x[2])))}
    temp <- optim(c(0,0),f,method="BFGS")$par
    alpha1_t[g] <- round(temp[1],4)
    alpha2_t[g] <- round(temp[2],4)
  }
  return(list(beta_t=beta_t,alpha1_t=alpha1_t,alpha2_t=alpha2_t))
}

gd_phi <- function(phi_0,beta,alpha1,alpha2,beta_t,alpha1_t,alpha2_t){
  phi <- phi_0
  P <- length(beta)
  phi[1:P] <- phi_0[1:P]+alpha1+alpha2
  phi[(P+1):(2*P)] <- phi_0[(P+1):(2*P)]+beta-beta_t
  phi[2*P-1+2*(1:P)] <- phi_0[2*P-1+2*(1:P)]+alpha1-alpha1_t
  phi[2*P+2*(1:P)] <- phi_0[2*P+2*(1:P)]+alpha2-alpha2_t
  return(phi)
}

admm <- function(y,x,w,beta_0,alpha1_0,alpha2_0,vlambda1,vlambda2,eps=1e-7,max=50){
  P <- length(beta)

  beta <- beta_0
  alpha1 <- alpha1_0
  alpha2 <- alpha2_0

  beta_t <- beta_0
  alpha1_t <- alpha1_0
  alpha2_t <- alpha2_0

  phi <- rep(0,4*P)

  convergence <- 1000
  i <- 1
  while (convergence>eps & i<max){
    new.theta <- gd_theta(y,x,w,beta,alpha1,alpha2,phi,beta_t,alpha1_t,alpha2_t)
    new.beta <- new.theta$beta
    new.alpha1 <- new.theta$alpha1
    new.alpha2 <- new.theta$alpha2

    new.psi <- glasso(new.beta,new.alpha1,new.alpha2,phi,vlambda1,vlambda2)
    new.beta_t <- new.psi$beta_t
    new.alpha1_t <- new.psi$alpha1_t
    new.alpha2_t <- new.psi$alpha2_t

    new.phi <- gd_phi(phi,new.beta,new.alpha1,new.alpha2,new.beta_t,new.alpha1_t,new.alpha2_t)

    convergence <- crossprod(beta-new.beta)+crossprod(alpha1-new.alpha1)+crossprod(alpha2-new.alpha2)+
      crossprod(beta_t-new.beta_t)+crossprod(alpha1_t-new.alpha1_t)+crossprod(alpha2_t-new.alpha2_t)+
      crossprod(phi-new.phi)

    beta <- new.beta
    alpha1 <- new.alpha1
    alpha2 <- new.alpha2
    beta_t <- new.beta_t
    alpha1_t <- new.alpha1_t
    alpha2_t <- new.alpha2_t
    phi <- new.phi

    i <- i+1
  }
  return(list(beta=beta_t,alpha1=alpha1_t,alpha2=alpha2_t))
}

linearHH <- function(formula,data,lambda,emmax=50,eps=1e-6){
  data <- na.omit(data) # remove subjects with missing covariates
  mf <- model.frame(formula, data)
  n <- dim(data)[1]
  y <- model.extract(mf, "response")
  x <- model.matrix(attr(mf, "terms"), mf)

  # Initial value from flexmix
  fullfit <- flexmix(formula,data,k=2)
  par <- parameters(fullfit)
  ini <- par[-nrow(par),]
  beta_0 <- apply(ini,1,mean)
  alpha1_0 <- ini[,1]-beta_0
  alpha2_0 <- ini[,2]-beta_0
  pi_0 <- prior(fullfit)

  # unpenalized MLEs are used as weights for adaptive LASSO
  vlambda1 <- lambda[1]/abs(beta_0)
  vlambda2 <- lambda[2]/sqrt(apply(rbind(alpha1_0,alpha2_0),2,crossprod))

  convergence <- 1000
  i <- 1
  beta <- beta_0
  alpha1 <- alpha1_0
  alpha2 <- alpha2_0
  pi <- pi_0

  temp1 <- pi[1]*dnorm(y,beta%*%t(x)+alpha1%*%t(x),1)
  temp2 <- pi[2]*dnorm(y,beta%*%t(x)+alpha2%*%t(x),1)
  w1 <- temp1/apply(rbind(temp1,temp2),2,sum)
  w2 <- temp2/apply(rbind(temp1,temp2),2,sum)
  w <- rbind(w1,w2)

  while (convergence>eps & i<emmax){

    #update theta = (beta,alpha1,alpha2)
    new.theta <- admm(y=y,x=x,w=w,beta_0=beta,alpha1_0=alpha1,alpha2_0=alpha2,vlambda1=vlambda1,vlambda2=vlambda2)
    new.beta <- new.theta$beta
    new.alpha1 <- new.theta$alpha1
    new.alpha2 <- new.theta$alpha2

    # update pi
    temp1 <- pi[1]*dnorm(y,new.beta%*%t(x)+new.alpha1%*%t(x),1)
    temp2 <- pi[2]*dnorm(y,new.beta%*%t(x)+new.alpha2%*%t(x),1)
    w1 <- temp1/apply(rbind(temp1,temp2),2,sum)
    w2 <- temp2/apply(rbind(temp1,temp2),2,sum)
    w <- rbind(w1,w2)
    new.pi <- apply(w,1,mean)

    convergence <- crossprod(beta-new.beta)+crossprod(alpha1-new.alpha1)+crossprod(alpha2-new.alpha2)+crossprod(pi-new.pi)
    beta <- new.beta
    alpha1 <- new.alpha1
    alpha2 <- new.alpha2
    pi <- new.pi

    i <- i+1
  }
  lik <- sum(log(pi[1]*dnorm(y,beta%*%t(x)+alpha1%*%t(x),1)+pi[2]*dnorm(y,beta%*%t(x)+alpha2%*%t(x),1)))
  BIC <- -2*lik+log(nrow(data))*(1+sum(beta!=0)+sum(alpha1!=0))
  return(list(beta=beta,alpha1=alpha1,alpha2=alpha2,BIC=BIC))
}

tuning <- function(vlambda1=1:5,vlambda2=1:5,formula,data){
  n1 <- length(vlambda1)
  n2 <- length(vlambda2)
  BIC <- matrix(0,n1,n2)
  for(i in 1:n1){
    for(j in 1:n2){
      lambda <- c(vlambda1[i],vlambda2[j])
      BIC[i,j] <- main_fun(formula,data,lambda)$BIC
      print(c(i,j))
    }
  }
  min.BIC <- min(BIC)
  loc <- which(BIC==min.BIC,arr.ind=T)
  optimal.lambda <- c(vlambda1[loc[1]],vlambda2[loc[2]])
  return(list(min.BIC=min.BIC,optimal.lambda=optimal.lambda))
}
