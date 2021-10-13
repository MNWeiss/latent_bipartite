require(runjags)
require(rjags)
require(asnipe)
require(aninet)

# define the latent space model for GoG data

latent_space_GoG <- "model{

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

# define parameters for simulation

G <- 100 # number of sampled groups
N <- 10 # number of individuals

# Make some simulated groups following Franks et al. 2010

# set up a preference matrix
B <- matrix(data = runif(N*N,-5,5), nrow = N, ncol = N)
diag(B) <- 0 
B[upper.tri(B)] <- t(B)[upper.tri(B)]

# function to sample groups
sample_groups <- function(B, M, mean_gs = 3){
  N <- nrow(B) # number of individuals
  groups <- list() # list to hold GBI for each sampling period
  for(i in 1:M){
    gs <- 0 # initialize group sizes
    for(j in 1:N){
      repeat{
        gs[j] <- rpois(1,mean_gs) # get a group size
        if(gs[j] != 0) break; # if it's a valid GS, break
      }
      if(sum(gs) >= N) break; # if we have group size totalling at least the pop size, break
    }
    if(sum(gs) > N){
      gs[length(gs)] <- gs[length(gs)] - (sum(gs)-N) #trim final group so total size of all groups matches pop size
    }
    gbi <- matrix(0, nrow = length(gs), ncol = N) # set up GBI
    while(any(rowSums(gbi) < gs)){
      s <- exp(gbi %*% B) # get preferences of individuals to groups
      s[rowSums(gbi) == gs,] <- 0 # full groups have 0 preference
      s.max <- apply(s,2,max) # get each individual's maximum preference score
      # choose an individual
      if(sum(colSums(gbi) == 0) == 1){
        a <- which(colSums(gbi) == 0)
      }else{
        a <- sample(which(colSums(gbi) == 0), 1, prob = s.max[colSums(gbi) == 0])
      }
      # allocate them to an empty group
      if(sum(rowSums(gbi) < gs) == 1){
        g <- which(rowSums(gbi) < gs)
      }else{
        g <- sample(which(rowSums(gbi) < gs), 1, prob = s[rowSums(gbi) < gs,a])
      }
      # record their group membership
      gbi[g,a] <- 1
    }
    groups[[i]] <- gbi # save the GBI
  }
  return(groups) # return the GBI
}

# generate all groups over a large number of sampling periods
all_groups <- sample_groups(B,10000,mean_gs = 3)

# full GBI
full_gbi <- do.call(rbind,all_groups)
# the "true" probabilities of co-occurrence in groups
true_SRI <- get_network(full_gbi)

# sampled/observed GBI
gbi <- do.call(rbind,
               lapply(all_groups[sample(length(all_groups),G)], function(z){
                 z[sample(nrow(z),1),,drop=F]
               }))

# empirical estimates of SRI from sampled GBI
est_SRI <- get_network(gbi)

# fit LS model
est_model <- run.jags(latent_space_GoG,
                      data = list(
                        gbi = gbi,
                        N = N,
                        D = 2,
                        R = diag(2),
                        G = nrow(gbi),
                        zero = rep(0,2)
                      ),
                      monitor = c("alpha","Z","W","iSw"),
                      adapt = 20000, sample = 50000)

# MCMC samples
est_mcmc <- do.call(rbind,est_model$mcmc)

require(cubature)
require(mvtnorm)

# Function to get LSAI
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
par(mfrow = c(1,2))
LSAI <- array(dim = c(100,N,N))
thin_samp <- sample(nrow(est_mcmc),100)
for(i in 1:length(thin_samp)){
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

SRI_b <- array(dim = c(1000,N,N))
for(i in 1:1000){
  gbi.b <- gbi[sample(nrow(gbi),rep=T),]
  SRI_b[i,,] <- get_network(gbi.b)
}

sri_x <- get_numerator(gbi,data_format="GBI",return="matrix")
sri_d <- get_denominator(gbi,data_format="GBI",return="matrix")

sri_a <- sri_x + 1
sri_b <- sri_d - sri_x + 1

SRI_bayes <- array(0,dim = c(1000,N,N))
for(i in 1:(N-1)){
  for(j in (i+1):N){
    SRI_bayes[,i,j] <- SRI_bayes[,j,i] <- rbeta(1000,sri_a[i,j],sri_b[i,j])
  }
}

str_b <- apply(SRI_b, 1, colSums)
str_bayes <- apply(SRI_bayes, 1, colSums)
str_LSAI <- apply(LSAI, 1, colSums)
str_true <- colSums(true_SRI)

str_b_ci <- apply(str_b, 1, quantile, probs = c(0.025,0.975))
str_bayes_ci <- apply(str_bayes, 1, quantile, probs = c(0.025,0.975))
str_LSAI_ci <- apply(str_LSAI, 1, quantile, probs = c(0.025,0.975),na.rm=T)

