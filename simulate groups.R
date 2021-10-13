generate_groups <- function(B, M, mean_gs){
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

generate_preferences <- function(N,delta){
  # set up a preference matrix
  B <- matrix(data = runif(N*N,-delta, delta), nrow = N, ncol = N)
  diag(B) <- 0 
  B[upper.tri(B)] <- t(B)[upper.tri(B)]
  return(B)
}

simulate_groups <- function(N,M,mean_gs,delta){
  B <- generate_preferences(N,delta)
  groups <- generate_groups(B,M,mean_gs=mean_gs)
  return(groups)
}

simulate_sampling <- function(groups, nsamp){
  do.call(rbind,
          lapply(groups[sample(length(all_groups),nsamp)], function(z){
            z[sample(nrow(z),1),,drop=F]
          }))
}

