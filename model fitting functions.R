require(rjags)
require(runjags)

get_LSI_MonteCarlo <- function(Z,alpha,sigma,nsim){
  N <- nrow(Z) # number of individuals
  LSI <- matrix(0, nrow = N, ncol = N) # matrix to hold LSAI
  D = ncol(Z) #dimensions of latent space
  groups <- matrix(data = rnorm(n = nsim*D, mean = 0, sd = sigma), ncol = D)
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

fit_LSI <- function(gbi, D = 2, nsamp = 1000, nsim = 10000, ...){
  
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
      Z[i,d] ~ dnorm(0,tau.ind)
    }
  }
  
  tau.ind ~ dgamma(1,0.1)
  
  # Generating group latent positions from MV normal
  for(g in 1:G){
    for(d in 1:D){
      W[g,d] ~ dnorm(0,tau.group)
    }
  }
  
  tau.group ~ dgamma(1,0.1)

 }"
  
  print("Estimating Model Parameters")
  
  est_model <- run.jags(modelstring,
                        data = list(
                          gbi = gbi,
                          N = N,
                          D = D,
                          G = nrow(gbi)
                        ),
                        monitor = c("alpha","Z","W","tau.group"),
                        ...)
  
  est_mcmc <- do.call(rbind,est_model$mcmc)
  
  print("Calculating LSI from model posterior")
  
  index <- seq(1,nrow(est_mcmc),length.out=nsamp)
  LSI <- array(dim = c(length(index),N,N))
  
  pb <- txtProgressBar(min = 0, max = length(index), style = 3)
  
  for(i in 1:length(index)){
    setTxtProgressBar(pb,i)
    Zest <- matrix(est_mcmc[index[i],substr(colnames(est_mcmc),1,1) == "Z"],nrow=N,ncol=D)
    alpha <- est_mcmc[index[i],"alpha"]
    tau.group <- est_mcmc[index[i],"tau.group"]
    sigma <- sqrt(1/tau.group)
    LSI[i,,] <- get_LSI_MonteCarlo(Zest,alpha,sigma,nsim=nsim)
  }
  
  return(LSI)
  
}

fit_MVB <- function(gbi, nsamp = 1000, nsim = 10000, ...){
  N <- ncol(gbi)
  G <- nrow(gbi)
  multivariate_bernoulli <- "model{

  for(g in 1:G){
    for(i in 1:N){
      gbi[g,i] ~ dbern(p[g,i])
      logit(p[g,i]) <- z[g,i]
    }
  }
  
  for(g in 1:G){
    z[g,1:N] ~ dmnorm(mu[1:N],iSigma[1:N,1:N])
  }
  
  for(i in 1:N){
    mu[i] ~ dnorm(0,0.001)
  }
  
  iSigma[1:N,1:N] ~ dwish(R[1:N,1:N],N+1)

  }"
  
  MVB_fit <- run.jags(model = multivariate_bernoulli,
                      data = list(
                        gbi = gbi,
                        G = G,
                        N = N,
                        R = diag(N)*(N+1)
                      ), monitor = c("mu","iSigma"),
                      ...)
  
  MVB_mcmc <- do.call(rbind, MVB_fit$mcmc)
  
  iSig <- array(data = MVB_mcmc[,substr(colnames(MVB_mcmc),1,2) == "iS"], dim = c(nrow(MVB_mcmc),N,N))
  Sigma <- iSig
  for(i in 1:nrow(MVB_mcmc)){
    Sigma[i,,] <- solve(iSig[i,,])
  }
  mu <- MVB_mcmc[,substr(colnames(MVB_mcmc),1,1) == "m"]
  MVB_index <- array(dim = c(nsamp,N,N))
  
  thin <- seq(from = 1, to = nrow(MVB_mcmc), length.out = nsamp)
  
  pb <- txtProgressBar(min = 0, max = nsamp, style = 3)
  
  for(i in 1:nsamp){
    setTxtProgressBar(pb,i)
    sim <- plogis(MASS::mvrnorm(n = nsim, mu = mu[thin[i],], Sigma = Sigma[thin[i],,]))
    pboth <- (t(sim) %*% sim)/nsim
    pneither <- (t(1-sim)%*%(1-sim))/nsim
    MVB_index[i,,] <- pboth/(1-pneither)
    diag(MVB_index[i,,]) <- 0
  }
  
  return(MVB_index)
  
}

fit_hub <- function(gbi, nsamp = 1000, ...){
  
  N <- ncol(gbi)
  G <- nrow(gbi)
  
  hub_model <- "model{

  for(g in 1:G){
    for(i in 1:N){
      gbi[g,i] ~ dbern(A[i,leader[g]])
    }
  }
  
  for(g in 1:G){
    leader[g] ~ dcat(p[1:N])
  }
  
  p[1:N] ~ ddirch(a[1:N])
  
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      A[j,i] <- A[i,j]
      A[i,j] ~ dunif(0,1)
    }
  }
  
  for(i in 1:N){
    A[i,i] <- 1
  }

  }"
  
  hub_fit <- run.jags(model = hub_model,
                      data = list(
                        gbi = gbi,
                        G = G,
                        N = N,
                        a = rep(1,N)
                      ), monitor = c("A","p"),
                      ...) 
  
  hub_mcmc <- do.call(rbind, hub_fit$mcmc)
  A <- array(data = hub_mcmc[,substr(colnames(hub_mcmc),1,1) == "A"], dim = c(nrow(hub_mcmc),N,N))
  p <- hub_mcmc[,substr(colnames(hub_mcmc),1,1) == "p"]
  
  hub_index <- array(0, dim = c(nsamp,N,N))
  
  thin <- seq(from = 1, to = nrow(hub_mcmc), length.out = nsamp)
  
  pb <- txtProgressBar(min = 0, max = nsamp, style = 3)
  
  for(k in 1:nsamp){
    setTxtProgressBar(pb,k)
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        both <- sum(p[k,]*A[k,i,]*A[k,j,])
        neither <- sum(p[k,]*(1-A[k,i,])*(1-A[k,j,]))
        hub_index[k,i,j] <- hub_index[k,j,i] <- both/(1-neither)
      }
    }
  }
  
  return(hub_index)
  
}