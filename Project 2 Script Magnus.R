library(rgl)
library(invgamma)

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


mcmc <- function(N){
  t0 <-proc.time()[3]
  
  alpha <- 2
  beta <- 0.05
  
  T <- length(rain$day)
  
  #pi <- matrix(nrow = N+1, ncol = T)
  #tau <- matrix(nrow = N+1, ncol = T)
  #pi[1,] <- rep(0.3,T)
  #tau[1,] <- pi_inv(pi[1,])
  
  pi <- rep(0.3,T)
  tau <- pi_inv(pi)
  
  
  
  Q <- matrix(data = 0,nrow = T, ncol = T)
  diag(Q) <- rep(2,T)
  diag(Q[-nrow(Q),-1]) <- rep(-1,T-1)
  diag(Q[-1,-ncol(Q)]) <- rep(-1,T-1)
  Q[1,1] <- 1
  Q[T,T] <- 1
  
  
  sigma_squared <- 0
  
  
  
  
  
  for (i in c(1:N)){
    sigma_squared <- rinvgamma(n = 1, shape = alpha + (T-1)/2, scale = beta + 1/2 * t(tau) %*% Q %*% tau )
    sigma <- sqrt(sigma_squared)
    
    
    for (d in c(1:T)){
      
      if (d == 1){
        new_tau_d <- rnorm(n=1, mean = tau[d+1], sd = sigma)
      }
      else if (d == T){
        new_tau_d <- rnorm(n=1, mean = tau[d-1], sd = sigma)
      }
      else{
        new_tau_d <- rnorm(n=1, mean = 1/2*(tau[d-1]+tau[d+1]), sd = sigma/sqrt(2))
      }
      
      
      new_tau <- tau
      new_tau[d] <- new_tau_d
      new_pi <- pi_func(new_tau)
      
      u <- runif(1)
      
      #logprob <- 0
      
      #for (d in c(1:T)){
      #  logprob = logprob + log(dbinom(rain$n.rain[d], size = rain$n.years[d], new_pi[d])) - log(dbinom(rain$n.rain[d], size = rain$n.years[d], pi[d]))
      #}
      
      prob <- dbinom(rain$n.rain[d], size = rain$n.years[d], new_pi[d])/dbinom(rain$n.rain[d], size = rain$n.years[d], pi[d])
      
      if (u < min(1, prob)){
        tau <- new_tau
        pi <- new_pi
      }
      else{
        
      }
      
    }
    
  }
  t <- proc.time()[3]
  print(t-t0)
  return(pi)
  
}


mcmc2 <- function(N){
  t0 <-proc.time()[3]
  
  alpha <- 2
  beta <- 0.05
  
  T <- length(rain$day)
  
  pi <- matrix(nrow = N+1, ncol = T)
  tau <- matrix(nrow = N+1, ncol = T)
  pi[1,] <- rep(0.3,T)
  tau[1,] <- pi_inv(pi[1,])
  
  
  
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
      
      
      new_tau <- tau[i,]
      new_tau[d] <- new_tau_d
      new_pi <- pi_func(new_tau)
      
      u <- runif(1)
      
      #logprob <- 0
      
      #for (d in c(1:T)){
      #  logprob = logprob + log(dbinom(rain$n.rain[d], size = rain$n.years[d], new_pi[d])) - log(dbinom(rain$n.rain[d], size = rain$n.years[d], pi[d]))
      #}
      
      prob <- dbinom(rain$n.rain[d], size = rain$n.years[d], new_pi[d])/dbinom(rain$n.rain[d], size = rain$n.years[d], pi[i,d])
      
      if (u < min(1, prob)){
        tau[(i+1),] <- new_tau
        pi[(i+1),] <- new_pi
      }
      else{
        tau[(i+1),] <- tau[i,]
        pi[(i+1),] <- pi[i,]
      }
      
    }
    
  }
  t <- proc.time()[3]
  print(t-t0)
  return(pi)
  
}







mcmcpi <- mcmc2(ceiling(500000/366))
ceiling(500000/366)
plot(mcmcpi[1366,])
plot(rain$n.rain/rain$n.years)

