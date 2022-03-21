library(rgl)
library(invgamma)
library(matlib)
library(MASS)

load("rain.rda")#Set working directory correctly
plot(rain)
plot(rain$day,rain$n.rain/rain$n.years)
sapply(rain,class)

#a
plot(rain$day,rain$n.rain)

print(rain$n.years[61])

#b

pi_func <- function(tau){
   1/(1+exp(-tau))
}

pi_inv <- function(pi){
  log(pi/(1-pi))
}

dnormal <- function(my,chol_sig,n) {
  d = length(my)
  #y <- matrix(nrow=n, ncol=d)
  x <- matrix(rnorm(d), d, n)
  A = chol_sig
  #my_matrix <- matrix(rep(my,n), nrow=d, ncol=n)
  #y <- my_matrix + A%*%x
  y <- my + A%*%x
  return(y)
}



mcmc <- function(N){
  t0 <-proc.time()[3]
  
  alpha <- 2
  beta <- 0.05
  
  T <- length(rain$day)
  
  tau <- matrix(nrow = N+1, ncol = T)
  tau[1,] <- rep(pi_inv(0.3), T)
  
  
  
  Q <- matrix(data = 0,nrow = T, ncol = T)
  diag(Q) <- rep(2,T)
  diag(Q[-nrow(Q),-1]) <- rep(-1,T-1)
  diag(Q[-1,-ncol(Q)]) <- rep(-1,T-1)
  Q[1,1] <- 1
  Q[T,T] <- 1
  
  
  sigma_squared <- c(1:N)*0
  
  
  for (i in c(1:N)){
    sigma_squared[i] <- rinvgamma(n = 1, shape = alpha + (T-1)/2, scale = beta + 1/2 * t(tau[i,]) %*% Q %*% tau[i,] )
    sigma <- sqrt(sigma_squared[i])
    
    
    for (d in c(1:T)){
      
      if (d == 1){
        new_tau_d <- rnorm(n=1, mean = tau[i,d+1], sd = sigma)
      }
      else if (d == T){
        new_tau_d <- rnorm(n=1, mean = tau[i,d-1], sd = sigma)
      }
      else{
        new_tau_d <- rnorm(n=1, mean = 1/2*(tau[i,d-1]+tau[i,d+1]), sd = sigma/sqrt(2))
      }
      
      old_tau_d <- tau[i, d]
      
      prob <- dbinom(rain$n.rain[d], size = rain$n.years[d], pi_func(new_tau_d))/dbinom(rain$n.rain[d], size = rain$n.years[d], pi_func(old_tau_d))
      
      u <- runif(1)
    
      if (u < min(1, prob)){
        tau[i, d] <- new_tau_d
      }
      else{
        #tau[i, d] <- old_tau_d
      }
    }
    
    tau[(i+1),] <- tau[i,]
  }
  
  t <- proc.time()[3]
  print(t-t0)
  
  apply(tau, 2, pi_func)
}




mcmc_block<- function(N,M){
  t0 <-proc.time()[3]
  
  alpha <- 2
  beta <- 0.05
  
  T <- length(rain$day)
  
  tau <- matrix(nrow = N+1, ncol = T)
  tau[1,] <- rep(pi_inv(0.3), T)
  
  
  
  Q <- matrix(data = 0,nrow = T, ncol = T)
  diag(Q) <- rep(2,T)
  diag(Q[-nrow(Q),-1]) <- rep(-1,T-1)
  diag(Q[-1,-ncol(Q)]) <- rep(-1,T-1)
  Q[1,1] <- 1
  Q[T,T] <- 1
  
  
  M_last <- T%%M
  
  Num_blocks <- T%/%M+1
  
  
  Q_first <- Q[1:M,1:M]
  QAB_first <- c(rep(0,M-1),-1)
  Q_inv_first <- inv(Q_first)
  chol_first <- t(chol(Q_inv_first))
  QQ_first <- -Q_inv_first%*%QAB_first
  
  Q_mid <- Q[(M+1):(2*M),(M+1):(2*M)]
  QAB_mid <- matrix(data = 0, nrow = M, ncol = 2)
  QAB_mid[1,1] <- -1
  QAB_mid[M,2] <- -1
  Q_inv_mid <- inv(Q_mid)
  chol_mid <- t(chol(Q_inv_mid))
  QQ_mid <- -Q_inv_mid%*%QAB_mid
  
  Q_last <- Q[(T-M_last+1):T,(T-M_last+1):T]
  QAB_last <- c(-1,rep(0,M_last-1))
  Q_inv_last <- inv(Q_last)
  chol_last <- t(chol(Q_inv_last))
  QQ_last <- -Q_inv_last%*%QAB_last
  
  
  
  sigma_squared <- c(1:N)*0
  
  
  for (i in c(1:N)){
    sigma_squared[i] <- rinvgamma(n = 1, shape = alpha + (T-1)/2, scale = beta + 1/2 * t(tau[i,]) %*% Q %*% tau[i,] )
    sigma <- sqrt(sigma_squared[i])
    
    
    for (d in c(1:Num_blocks)){
      
      if (d == 1){
        a <- 1
        b <- M
        #new_tau_d <- mvrnorm(n=1, mu = QQ_first%*%tau[i,b+1], Sigma = sigma_squared[i]*Q_inv_first)
        new_tau_d <- dnormal(my = QQ_first%*%tau[i,b+1], chol_sig = sigma*chol_first, n = 1)
        #new_tau_d <- rnorm(n=1, mean = tau[i,d+1], sd = sigma)
      }
      else if (d == Num_blocks){
        a <- T - M_last + 1
        b <- T
        #new_tau_d <- mvrnorm(n=1, mu = QQ_last%*%tau[i,a-1], Sigma = sigma_squared[i]*Q_inv_last)
        new_tau_d <- dnormal(my = QQ_last%*%tau[i,a-1], chol_sig = sigma*chol_last, n = 1)
        #new_tau_d <- rnorm(n=1, mean = tau[i,d-1], sd = sigma)
      }
      else{
        a <- (d-1)*M + 1
        b <- d*M
        #new_tau_d <- mvrnorm(n=1, mu = QQ_mid%*%c(tau[i,a-1],tau[i,b+1]), Sigma = sigma_squared[i]*Q_inv_mid)
        new_tau_d <- dnormal(my = QQ_mid%*%c(tau[i,a-1],tau[i,b+1]), chol_sig = sigma*chol_mid, n = 1)
        #new_tau_d <- rnorm(n=1, mean = 1/2*(tau[i,d-1]+tau[i,d+1]), sd = sigma/sqrt(2))
      }
      
      
      old_tau <- tau[i,]
      new_tau <- tau[i,]
      new_tau[a:b] <- new_tau_d
      
      logprob <- 0
      
      for (j in c(a:b)){
        logprob = logprob + log(dbinom(rain$n.rain[j], size = rain$n.years[j], pi_func(new_tau[j]))) - log(dbinom(rain$n.rain[j], size = rain$n.years[j], pi_func(old_tau[j])))
      }
      
      
      u <- runif(1)
      
      if (u < min(1, exp(logprob))){
        tau[i,] <- new_tau
      }
      else{
        #tau[i, d] <- old_tau_d
      }
    }
    
    tau[(i+1),] <- tau[i,]
  }
  
  t <- proc.time()[3]
  print(t-t0)
  
  apply(tau, 2, pi_func)
}

pipipi <- mcmc_block(1000,100)
plot(pipipi[1000,])

mcmcpi <- mcmc(1000)
plot(mcmcpi[1000,])

plot(rain$n.rain/rain$n.years)


#2##2##2##2##2##2##2##2##2##2##2##2#


library(INLA)








