require(runjags)
require(rjags)

get_LSAI <- function(Z, alpha, Sigma){
  N <- nrow(Z) # number of individuals
  LSAI <- matrix(0, nrow = N, ncol = N) # matrix to hold LSAI
  D = ncol(Z) #dimensions of latent space
  for(i in 1:(N-1)){
    print(i)
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
      LSAI[i,j] <- LSAI[j,i] <- pboth/(1 - pneither)
    }
  }
  return(LSAI)
}

fit_LSAI <- function(gbi, D = 2, nthin = 100, ...){
 
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
  
est_model <- run.jags(latent_space_GoG,
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
  
  LSAI <- array(dim = c(nthin,N,N))
  thin_samp <- sample(nrow(est_mcmc),nthin)
  
  for(i in 1:nthin){
    s <- thin_samp[i]
    Zest <- matrix(est_mcmc[s,substr(colnames(est_mcmc),1,1) == "Z"],nrow=N,ncol=D)
    alpha <- est_mcmc[s,"alpha"]
    Sigma <- solve(matrix(est_mcmc[s,substr(colnames(est_mcmc),1,3) == "iSw"],nrow=D,ncol=D))
    LSAI[i,,] <- get_LSAI(Zest,alpha,Sigma)
    plot(apply(LSAI,c(2,3),mean,na.rm=T) ~ true_SRI)
    abline(a = 0, b = 1)
    plot(est_SRI ~ true_SRI)
    abline(a = 0, b = 1)
  }
  
  return(LSAI)
  
}