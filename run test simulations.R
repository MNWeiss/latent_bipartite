source("simulate groups.R")
source("LSI fitting.R")
source("bootstrap methods.R")

N <- c(10,20,50)
ngroups <- c(100, 200, 500)
mean_gs = c(2,5)
delta = c(1,5,10)

params <- expand.grid(N=N,sampled_groups=ngroups,mean_gs=mean_gs,delta=delta)

params$LSI_dyad <- params$bayes_dyad <- params$boot_dyad <- NA
params$LSI_str <- params$bayes_str <- params$boot_str <- NA
params$LSI_cl <- params$bayes_cl <- params$boot_cl <- NA
params$LSI_mean <- params$bayes_mean <- params$boot_mean <- NA
params$LSI_eff <- params$bayes_eff <- params$boot_eff <- NA

for(i in 1:nrow(params)){
  
  err_LSI_dyad <- err_bayes_dyad <- err_boot_dyad <- NA
  err_LSI_str <- err_bayes_str <- err_boot_str <- NA
  err_LSI_cl <- err_bayes_cl <- err_boot_cl <- NA
  err_LSI_mean <- err_bayes_mean <- err_boot_mean <- NA
  err_LSI_eff <- err_bayes_eff <- err_boot_eff <- NA
  
  for(s in 1:20){
    # Generate data
    
    groups <- simulate_groups(N = params$N[i], M = 100000, mean_gs = params$mean_gs[i], delta = params$delta[i])
    full_gbi <- do.call(rbind,groups)
    gbi <- simulate_sampling(groups,nsamp=params$sampled_groups[i])
    
    # Get some SRI measures
    
    true_SRI <- asnipe::get_network(full_gbi)
    emp_SRI <- asnipe::get_network(gbi)
    boot_SRI <- nonparametric_bootstrap(gbi,nboot = 1000)
    bayes_SRI <- bayesian_bootstrap(gbi,nboot = 1000)
    
    # Fit the LSI
    
    LSI <- fit_LSI(gbi,D=2,thin=400,adapt=1000,sample=50000,n.chains=4)
    
    # Get CIs and error for each dyad
    
    LSI_CI <- apply(LSI, c(2,3), quantile, probs = c(0.025, 0.975))
    bayes_CI <- apply(bayes_SRI, c(2,3), quantile, probs = c(0.025, 0.975))
    boot_CI <- apply(boot_SRI, c(2,3), quantile, probs = c(0.025, 0.975))
    
    err_LSI_dyad[s] <- mean(true_SRI[lower.tri(true_SRI)] < LSI_CI[1,,][lower.tri(true_SRI)] | true_SRI[lower.tri(true_SRI)] > LSI_CI[2,,][lower.tri(true_SRI)])
    err_bayes_dyad[s] <- mean(true_SRI[lower.tri(true_SRI)] < bayes_CI[1,,][lower.tri(true_SRI)] | true_SRI[lower.tri(true_SRI)] > bayes_CI[2,,][lower.tri(true_SRI)])
    err_boot_dyad[s] <- mean(true_SRI[lower.tri(true_SRI)] < boot_CI[1,,][lower.tri(true_SRI)] | true_SRI[lower.tri(true_SRI)] > boot_CI[2,,][lower.tri(true_SRI)])
    
    # Get CIs and error for strength
    
    LSI_str_CI <- apply(apply(LSI,c(1,2),sum),1,quantile,probs=c(0.025,0.975))
    bayes_str_CI <- apply(apply(bayes_SRI,c(1,2),sum),1,quantile,probs=c(0.025,0.975))
    boot_str_CI <- apply(apply(boot_SRI,c(1,2),sum),1,quantile,probs=c(0.025,0.975))
    true_str <- colSums(true_SRI)
    
    err_LSI_str[s] <- mean(true_str > LSI_str_CI[2,] | true_str < LSI_str_CI[1,])
    err_bayes_str[s] <- mean(true_str > bayes_str_CI[2,] | true_str < bayes_str_CI[1,])
    err_boot_str[s] <- mean(true_str > boot_str_CI[2,] | true_str < boot_str_CI[1,])
    
    # Get CIs and error for closeness
    
    true_g <- graph.adjacency(true_SRI, weighted = T, mode = "undirected")
    true_d <- igraph::distances(true_g,weights = 1/E(true_g)$weight)
    diag(true_d) <- NA
    
    LSI_cl_CI <- apply(apply(LSI, 1, function(z){
      g <- graph.adjacency(z, weighted = T, mode = "undirected")
      d <- igraph::distances(g, weights = 1/E(g)$weight)
      diag(d) <- NA
      e <- 1/d
      colMeans(e,na.rm = T)
    }),1,quantile,probs=c(0.025,0.975))
    bayes_cl_CI <- apply(apply(bayes_SRI, 1, function(z){
      g <- graph.adjacency(z, weighted = T, mode = "undirected")
      d <- igraph::distances(g, weights = 1/E(g)$weight)
      diag(d) <- NA
      e <- 1/d
      colMeans(e,na.rm = T)
    }),1,quantile,probs=c(0.025,0.975))
    boot_cl_CI <- apply(apply(boot_SRI, 1, function(z){
      g <- graph.adjacency(z, weighted = T, mode = "undirected")
      d <- igraph::distances(g, weights = 1/E(g)$weight)
      diag(d) <- NA
      e <- 1/d
      colMeans(e,na.rm = T)
    }),1,quantile,probs=c(0.025,0.975))
    true_cl <- rowMeans(1/true_d,na.rm=T)
    
    err_LSI_cl[s] <- mean(true_cl < LSI_cl_CI[1,] | true_cl > LSI_cl_CI[2,])
    err_bayes_cl[s] <- mean(true_cl < bayes_cl_CI[1,] | true_cl > bayes_cl_CI[2,])
    err_boot_cl[s] <- mean(true_cl < boot_cl_CI[1,] | true_cl > boot_cl_CI[2,])
    
    
    # Get CIs and error for mean association probability
    
    LSI_mean_CI <- quantile(apply(LSI, 1, function(z) mean(z[lower.tri(z)])  ),c(0.025,0.975))
    bayes_mean_CI <- quantile(apply(bayes_SRI, 1, function(z) mean(z[lower.tri(z)])  ),c(0.025,0.975))
    boot_mean_CI <- quantile(apply(boot_SRI, 1, function(z) mean(z[lower.tri(z)])  ),c(0.025,0.975))
    true_mean <- mean(true_SRI[lower.tri(true_SRI)])
    
    err_LSI_mean[s] <- mean(true_mean < LSI_mean_CI[1] | true_mean > LSI_mean_CI[2])
    err_bayes_mean[s] <- mean(true_mean < bayes_mean_CI[1] | true_mean > bayes_mean_CI[2])
    err_boot_mean[s] <- mean(true_mean < boot_mean_CI[1] | true_mean > boot_mean_CI[2])
    
    # Get CIs and error for network efficiency
    
    LSI_eff_CI <- quantile(apply(LSI, 1, function(z){
      g <- graph.adjacency(z, weighted = T, mode = "undirected")
      d <- igraph::distances(g, weights = 1/E(g)$weight)
      mean(1/d[lower.tri(d)])
    }), c(0.025,0.975))
    bayes_eff_CI <- quantile(apply(bayes_SRI, 1, function(z){
      g <- graph.adjacency(z, weighted = T, mode = "undirected")
      d <- igraph::distances(g, weights = 1/E(g)$weight)
      mean(1/d[lower.tri(d)])
    }), c(0.025,0.975))
    boot_eff_CI <- quantile(apply(boot_SRI, 1, function(z){
      g <- graph.adjacency(z, weighted = T, mode = "undirected")
      d <- igraph::distances(g, weights = 1/E(g)$weight)
      mean(1/d[lower.tri(d)])
    }), c(0.025,0.975))
    true_eff <- mean(1/true_d[lower.tri(true_d)])
    
    err_LSI_eff[s] <- mean(true_eff < LSI_eff_CI[1] | true_eff > LSI_eff_CI[2])
    err_bayes_eff[s] <- mean(true_eff < bayes_eff_CI[1] | true_eff > bayes_eff_CI[2])
    err_boot_eff[s] <- mean(true_eff < boot_eff_CI[1] | true_eff > boot_eff_CI[2])
  }
  
  params$LSI_dyad[i] <- mean(err_LSI_dyad)
  params$bayes_dyad[i] <- mean(err_boot_dyad)
  params$boot_dyad[i] <- mean(err_bayes_dyad)
  
  params$LSI_str[i] <- mean(err_LSI_str)
  params$bayes_str[i] <- mean(err_boot_str)
  params$boot_str[i] <- mean(err_bayes_str)
  
  params$LSI_cl[i] <- mean(err_LSI_cl)
  params$bayes_cl[i] <- mean(err_boot_cl)
  params$boot_cl[i] <- mean(err_bayes_cl)
  
  params$LSI_mean[i] <- mean(err_LSI_mean)
  params$bayes_mean[i] <- mean(err_boot_mean)
  params$boot_mean[i] <- mean(err_bayes_mean)
  
  params$LSI_eff[i] <- mean(err_LSI_eff)
  params$bayes_eff[i] <- mean(err_boot_eff)
  params$boot_eff[i] <- mean(err_bayes_eff)
  
}