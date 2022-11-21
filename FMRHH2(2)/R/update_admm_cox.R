#' Internal function for admm_cox
#' #' Identify of homogeneous and heterogeneous variables in a mixture of Cox regression model
#'
#' Fit a mixture of  Cox regression model using EM and ADMM algorithm
#' @param formula a formula object, with the response in the left of a ~, and the terms on the right. The response must be a survival object as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the formula.
#' @param lambda specify the tuning parameters. The first number is the tuning for average effect beta, and the second number is the tuning for deviance effect alpha.
#' @param eps the convergence criteria.
#' @param max maximum number of iterations.
#'
#' @return coefficients estimates of beta and alpha
#' @import survival
#' @export
#'

update_theta_cox <- function(t,status,x,index,beta,alpha1,alpha2,alpha3,phi,beta_t,alpha1_t,alpha2_t,alpha3_t){
  P <- length(beta)
  n1 <- sum(index==1)
  n2 <- sum(index==2)
  n3 <- sum(index==3)
  f <- function(theta){
    event.index1 <- which(status[index==1]==1)
    event.index2 <- which(status[index==2]==1)
    event.index3 <- which(status[index==3]==1)
    temp1 <- 0
    for(i in 1:length(event.index1)){
      temp1 <- temp1+(x[index==1,])[event.index1[i], ]%*%(theta[1:P]+theta[(P+1):(2*P)])-log(sum(exp((x[index==1,])[event.index1[i]:n1,]%*%(theta[1:P]+theta[(P+1):(2*P)]))))
    }
    for(i in 1:length(event.index2)){
      temp1 <- temp1+(x[index==2,])[event.index2[i], ]%*%(theta[1:P]+theta[(2*P+1):(3*P)])-log(sum(exp((x[index==2,])[event.index2[i]:n2,]%*%(theta[1:P]+theta[(2*P+1):(3*P)]))))
    }
    for(i in 1:length(event.index3)){
      temp1 <- temp1+(x[index==3,])[event.index3[i], ]%*%(theta[1:P]+theta[(3*P+1):(4*P)])-log(sum(exp((x[index==3,])[event.index3[i]:n3,]%*%(theta[1:P]+theta[(3*P+1):(4*P)]))))
    }
    temp2 <- phi[1:P]*apply(rbind(theta[(P+1):(2*P)],theta[(2*P+1):(3*P)],theta[(3*P+1):(4*P)]),2,sum)
    +phi[(P+1):(2*P)]*theta[1:P]+phi[2*P+3*(1:P)-2]*theta[(P+1):(2*P)]+phi[2*P+3*(1:P)-1]*theta[(2*P+1):(3*P)]+phi[2*P+3*(1:P)]*theta[(3*P+1):(4*P)]
    temp3 <- (apply(rbind(theta[(P+1):(2*P)],theta[(2*P+1):(3*P)],theta[(3*P+1):(4*P)]),2,sum))^2+(theta-theta_t)^2
    -temp1+sum(temp2)+1/2*sum(temp3)
  }
  theta_t <- c(beta_t,alpha1_t,alpha2_t,alpha3_t)
  out <- optim(c(beta,alpha1,alpha2,alpha3),f,method="BFGS")$par
  return(list(beta=out[1:P],alpha1=out[(P+1):(2*P)],alpha2=out[(2*P+1):(3*P)],alpha3=out[(3*P+1):(4*P)]))
}


#' Internal function for admm_cox
#'
#' @param beta
#' @param alpha1
#' @param alpha2
#' @param alpha3
#' @param phi
#' @param vlambda1
#' @param vlambda2
#'
#' @description an internal function to update thetat in admm_cox
#' @keywords internal
#'
update_thetat_cox <- function(beta,alpha1,alpha2,alpha3,phi,vlambda1,vlambda2){
  P <- length(beta)
  beta_t <- rep(0,P)
  alpha1_t <- rep(0,P)
  alpha2_t <- rep(0,P)
  alpha3_t <- rep(0,P)
  # update beta_t
  for(g in 1:P){
    f <- function(x){0.5*(x-beta[g]-phi[P+g])^2+vlambda1[g]*abs(x)}
    beta_t[g] <- round(optim(0,f,method="BFGS")$par,4)
  }

  # update alpha_t
  for(g in 1:P){
    f <- function(x){0.5*(x[1]-alpha1[g]-phi[2*P+3*g-2])^2+0.5*(x[2]-alpha2[g]-phi[2*P+3*g-1])^2+0.5*(x[3]-alpha3[g]-phi[2*P+3*g])^2+vlambda2[g]*sqrt(crossprod(c(x[1],x[2],x[3])))}
    temp <- optim(c(0,0,0),f,method="BFGS")$par
    alpha1_t[g] <- round(temp[1],4)
    alpha2_t[g] <- round(temp[2],4)
    alpha3_t[g] <- round(temp[3],4)
  }
  return(list(beta_t=beta_t,alpha1_t=alpha1_t,alpha2_t=alpha2_t,alpha3_t=alpha3_t))
}


#' Internal function for admm_cox
#' @description an internal function to update thetat in admm_cox
#' @keywords internal
#'
update_phi_cox <- function(phi_0,beta,alpha1,alpha2,alpha3,beta_t,alpha1_t,alpha2_t,alpha3_t){
  phi <- phi_0
  P <- length(beta)
  phi[1:P] <- phi_0[1:P]+alpha1+alpha2+alpha3
  phi[(P+1):(2*P)] <- phi_0[(P+1):(2*P)]+beta-beta_t
  phi[2*P+3*(1:P)-2] <- phi_0[2*P+3*(1:P)-2]+alpha1-alpha1_t
  phi[2*P+3*(1:P)-1] <- phi_0[2*P+3*(1:P)-1]+alpha2-alpha2_t
  phi[2*P+3*(1:P)] <- phi_0[2*P+3*(1:P)]+alpha3-alpha3_t
  return(phi)
}
