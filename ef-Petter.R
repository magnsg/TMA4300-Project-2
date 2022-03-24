library(rgl)
library(invgamma)
library(matlib)
library(MASS)

load("rain.rda")#Set working directory correctly
sapply(rain,class)

# Tau to pi
pi_func <- function(tau){
  1/(1+exp(-tau))
}

# Pi to tau
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
  accepted_count <- 0
  count <- 0
  
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
  
  
  sigma_squared <- 1:N
  
  
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
        accepted_count <- accepted_count + 1
      }
    }
    
    tau[(i+1),] <- tau[i,]
  }
  
  print('Processing time:')
  print(proc.time()[3] - t0)
  
  print('Acceptance probability:')
  print(accepted_count / (366*N))
  
  list(
    apply(tau, 2, pi_func),
    sigma_squared
  )
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
  
  list(
    apply(tau, 2, pi_func),
    sigma_squared
  )
}

# Equi-taled credible interval. 
# Finds an interval where we estimate a minimum of 95% chance of getting new values
cred_int <- function(vec){
  vec <- sort(vec)
  len <- length(vec)
  size <- ceiling(0.95 * len)
  start_index <- floor((len - size) / 2) + 1

  list(vec[start_index], vec(start_index + size))
}

# Simulation
N = 1000
mcmc_sim <- mcmc(N)
sim_pi <- matrix(unlist(mcmc_sim[1]), ncol = 366)
sim_sigsq <- mcmc_sim[2]

# Traceplots
par(mfrow = c(4, 1))
for (day in c(1, 201, 366)){
  plot(sim_pi[, day], main = c('Rain probability trace day ', day), xlab = 'Simulation #', ylab = 'Rain probability')
}
plot(1:N, sim_sigsq, main = 'Sigma trace', xlab = 'Simulation #', ylab = 'Sigma squared', ylim = c(0, 0.01))

# Histogram
par(mfrow = c(4, 1))
for (day in c(1, 201, 366)){
  hist(sim_pi[, day], main = c('Day ', day, ' rain probability histogram') freq = TRUE, xlab = 'Rain probability', ylab = 'Frequency')
}
hist(sim_sigsq, main = 'Sigma squared histogram' freq = TRUE, xlab = 'Sigma squared', ylab = 'Frequency', xlim = c(-0.1, 0.1), breaks = 0.001)

# Autocorrelation
par(mfrow = c(4, 1))
for (day in c(1, 201, 366)){
  acf_plot <- acf(sim_pi[, day], lag.max = 2000, plot = FALSE)
  plot(acf_plot, main = c('Autocorrelation for day ', day),  xlab = 'Autocorrelation interval', ylab = 'Correlation')
}
acf_plot <- acf(sim_sigsq, lag.max = 50, plot = FALSE)
plot(acf_plot, main = 'Autocorrelation sigma squared', xlab = 'Autocorrelation interval', ylab = 'Correlation', ylim = c(-0.03, 0.03))

# Means
print('Means')
for (day in c(1, 201, 366)){
  print(c('Day ', day, ': ', mean(sim_pi[, day])))
}
print('Sigma: ', mean(sim_sigsq))

# 95% credible intervals
print('95 % credible intervals')
for (day in c(1, 201, 366)){
  print(paste('Day ', day, ': '))
  print(paste(cred_int(sim_pi[, day])))
}
print(paste('Sigma: ', cred_int(sim_sigsq)[1], ' - ', cred_int(sim_sigsq)[2]))

# Plot to compare tau (red) and day rain probability (blue)
plot(rain$n.rain/rain$n.years, col = 'Blue')
points(sim_pi[N, ], col = 'Red')

