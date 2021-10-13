source("LSAI fitting.R")
source("simulate groups.R")
source("bootstrap methods.R")

full_simulation <- function(N,groups,mean_gs,delta,
                            D, nthin,
                            prior, nboot){
  
  all_groups <- simulate_groups(N = N, M = 100000, mean_gs = mean_gs, delta = delta)
  gbi <- simulate_sampling(all_groups,groups)
  
  LSAI <- fit_LSAI(gbi, D = D, nthin = nthin)
  bayesSRI <- bayesian_bootstrap(gbi,nboot=nboot,prior=prior)
  bootSRI <- nonparametric_bootstrap(gbi,nboot=nboot)
  trueSRI <- asnipe::get_network(do.call(rbind,all_groups))
  
  return(list(
    trueSRI = trueSRI,
    bayesSRI = bayesSRI,
    bootSRI = bootSRI,
    LSAI = LSAI
  ))
  
}