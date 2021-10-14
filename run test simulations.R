source("simulate groups.R")
source("LSI fitting.R")
source("bootstrap methods.R")

N <- c(10,20,50)
ngroups <- c(50, 100, 200, 500)
mean_gs = c(2,5)
delta = c(1,5,10)

params <- expand.grid(N=N,sampled_groups=ngroups,mean_gs=mean_gs,delta=delta)

for(i in 1:nrow(params)){
  groups <- simulate_groups(N = params$N[i], M = 50000, mean_gs = params$mean_gs[i], delta = params$delta[i])
  full_gbi <- do.call(rbind,groups)
  gbi <- simulate_sampling(groups,nsamp=params$sampled_groups[i])
  true_SRI <- asnipe::get_network(full_gbi)
  emp_SRI <- asnipe::get_network(gbi)
  boot_SRI <- nonparametric_bootstrap(gbi,nboot = 1000)
  bayes_SRI <- bayesian_bootstrap(gbi,nboot = 1000)
  LSI <- fit_LSI(gbi,D=2,thin=100,adapt=50000,sample=20000,n.chains=4)
}