require(runjags)
require(rjags)
require(cubature)
require(mvtnorm)

# functions to calculate LSI from MCMC output of a latent space model

get_LSI_integral <- function(Z, alpha, Sigma){
  N <- nrow(Z) # number of individuals
  LSI <- matrix(0, nrow = N, ncol = N) # matrix to hold LSAI
  D = ncol(Z) #dimensions of latent space
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      z1 <- Z[i,] # location of i
      z2 <- Z[j,] # location of j
      # functions to integrate over
      fneither <- function(x){
        p1 <- plogis(alpha - sqrt(sum((x-z1)^2)))
        p2 <- plogis(alpha - sqrt(sum((x-z2)^2)))
        dmvnorm(x, mean = rep(0,D), sigma = Sigma)*(1-p1)*(1-p2)
      }
      fboth <- function(x){
        p1 <- plogis(alpha - sqrt(sum((x-z1)^2)))
        p2 <- plogis(alpha - sqrt(sum((x-z2)^2)))
        dmvnorm(x, mean = rep(0,D), sigma = Sigma)*p1*p2
      }
      # probability of i and j
      pboth <- hcubature(fboth, lowerLimit = rep(-Inf,D), upperLimit = rep(Inf,D), tol = 1e-4)$integral
      # probability of neither individual
      pneither <- hcubature(fneither, lowerLimit = rep(-Inf,D), upperLimit = rep(Inf,D), tol = 1e-4)$integral
      # latent space association index
      LSI[i,j] <- LSI[j,i] <- pboth/(1 - pneither)
    }
  }
  return(LSI)
}
get_LSI_MonteCarlo <- function(Z,alpha,Sigma,nsim){
  N <- nrow(Z) # number of individuals
  LSI <- matrix(0, nrow = N, ncol = N) # matrix to hold LSAI
  D = ncol(Z) #dimensions of latent space
  groups <- MASS::mvrnorm(nsim,rep(0,D),Sigma) # get simulated groups
  # get the probability that each individual is in each simulated group
  p_in <- apply(Z, 1, function(z){
    apply(groups, 1, function(w){
      plogis(alpha - sqrt(sum((z-w)^2)))
    })
  })
  pboth <- (t(p_in)%*%p_in)/nsim # probability of both
  pneither <- (t(1-p_in)%*%(1-p_in))/nsim #probability of either
  LSI <- pboth/(1-pneither) # LSI
  diag(LSI) <- 0
  return(LSI)
}

# function to fit the model and get the LSI matrices

fit_LSI <- function(gbi, D = 2, thin = 50, method = "MC", ...){
 
  G <- nrow(gbi)
  N <- ncol(gbi)
  
  modelstring <- "model{

  for(g in 1:G){
    for(i in 1:N){
      
      gbi[g,i] ~ dbern(p[g,i]) # observation process
      logit(p[g,i]) <- alpha - d[g,i] # get probability of joining group
      d[g,i] <- sqrt(sum((Z[i,1:D]-W[g,1:D])^2)) # get distance of i to g
      
    }
  }
  
  alpha ~ dnorm(0,0.01) #prior on intercept
  
  # Priors on individual latent positions
  for(i in 1:N){
    for(d in 1:D){
      Z[i,d] ~ dnorm(0,0.01)
    }
  }
  
  # Generating group latent positions from MV normal
  for(g in 1:G){
    W[g,1:D] ~ dmnorm(zero[1:D],iSw[1:D,1:D])
  }
  
  iSw[1:D,1:D] ~ dwish(R[1:D,1:D],D+1) # prior on group VCV matrix

 }"
  
 print("Estimating Model Parameters")
  
 est_model <- run.jags(modelstring,
                        data = list(
                          gbi = gbi,
                          N = N,
                          D = D,
                          R = diag(D),
                          G = nrow(gbi),
                          zero = rep(0,D)
                        ),
                        monitor = c("alpha","Z","W","iSw"),
                        ...)
  
  est_mcmc <- do.call(rbind,est_model$mcmc)
  
  print("Calculating LSI from model posterior")
  
  index <- seq(thin,nrow(est_mcmc),thin)
  LSI <- array(dim = c(length(index),N,N))
  
  pb <- txtProgressBar(min = 0, max = length(index), style = 3)
  
  for(i in 1:length(index)){
    setTxtProgressBar(pb,i)
    Zest <- matrix(est_mcmc[index[i],substr(colnames(est_mcmc),1,1) == "Z"],nrow=N,ncol=D)
    alpha <- est_mcmc[index[i],"alpha"]
    Sigma <- solve(matrix(est_mcmc[index[i],substr(colnames(est_mcmc),1,3) == "iSw"],nrow=D,ncol=D))
    if(method == "MC"){
      LSI[i,,] <- get_LSI_MonteCarlo(Zest,alpha,Sigma,nsim=10000)
    }else{
      LSI[i,,] <- get_LSI_integral(Zest,alpha,Sigma)
    }
  }
  
  return(LSI)
  
}
