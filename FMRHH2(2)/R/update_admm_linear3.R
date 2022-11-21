#' Internal function for admm_linear3
#' @description an internal function to update theta in admm through gradient descent
#' @keywords internal
#'
gd_theta3 <- function(y,x,w,beta,alpha1,alpha2,alpha3,phi,beta_t,alpha1_t,alpha2_t,alpha3_t){
  P <- length(beta)
  f <- function(theta){
    temp1 <- (y-theta[1:P]%*%t(x)-theta[(P+1):(2*P)]%*%t(x))^2*w[1,]/2
    +(y-theta[1:P]%*%t(x)-theta[(2*P+1):(3*P)]%*%t(x))^2*w[2,]/2
    +(y-theta[1:P]%*%t(x)-theta[(3*P+1):(4*P)]%*%t(x))^2*w[3,]/2
    temp2 <- phi[1:P]*apply(rbind(theta[(P+1):(2*P)],theta[(2*P+1):(3*P)],theta[(3*P+1):(4*P)]),2,sum)
    +phi[(P+1):(2*P)]*theta[1:P]+phi[2*P+3*(1:P)-2]*theta[(P+1):(2*P)]+phi[2*P+3*(1:P)-1]*theta[(2*P+1):(3*P)]+phi[2*P+3*(1:P)]*theta[(3*P+1):(4*P)]
    temp3 <- (apply(rbind(theta[(P+1):(2*P)],theta[(2*P+1):(3*P)],theta[(3*P+1):(4*P)]),2,sum))^2+(theta-theta_t)^2
    sum(temp1)+sum(temp2)+1/2*sum(temp3)
  }

  theta_t <- c(beta_t,alpha1_t,alpha2_t,alpha3_t)
  out <- optim(c(beta,alpha1,alpha2,alpha3),f,method="BFGS")$par
  return(list(beta=out[1:P],alpha1=out[(P+1):(2*P)],alpha2=out[(2*P+1):(3*P)],alpha3=out[(3*P+1):(4*P)]))
}

#' Internal function for admm
#' @description an internal function for group lasso in admm
#' @keywords internal
#'
glasso3 <- function(beta,alpha1,alpha2,alpha3,phi,vlambda1,vlambda2){
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

#' Internal function for admm
#' @description an internal function to update phi in admm through gradient descent
#' @keywords internal
#'
gd_phi3 <- function(phi_0,beta,alpha1,alpha2,alpha3,beta_t,alpha1_t,alpha2_t,alpha3_t){
  phi <- phi_0
  P <- length(beta)
  phi[1:P] <- phi_0[1:P]+alpha1+alpha2
  phi[(P+1):(2*P)] <- phi_0[(P+1):(2*P)]+beta-beta_t
  phi[2*P+3*(1:P)-2] <- phi_0[2*P+3*(1:P)-2]+alpha1-alpha1_t
  phi[2*P+3*(1:P)-1] <- phi_0[2*P+3*(1:P)-1]+alpha2-alpha2_t
  phi[2*P+3*(1:P)] <- phi_0[2*P+3*(1:P)]+alpha3-alpha3_t
  return(phi)
}
