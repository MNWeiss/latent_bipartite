nonparametric_bootstrap <- function(gbi,nboot){
  N <- ncol(gbi)
  SRI_b <- array(dim = c(nboot,N,N))
  for(i in 1:nboot){
    gbi.b <- gbi[sample(nrow(gbi),rep=T),]
    SRI_b[i,,] <- asnipe::get_network(gbi.b)
  }
  return(SRI_b)
}

bayesian_bootstrap <- function(gbi,nboot,prior = c(1,1)){
  N <- ncol(gbi)
  sri_x <- aninet::get_numerator(gbi,data_format="GBI",return="matrix")
  sri_d <- aninet::get_denominator(gbi,data_format="GBI",return="matrix")
  sri_a <- sri_x + prior[1]
  sri_b <- sri_d - sri_x + prior[2]
  SRI_bayes <- array(0,dim = c(1000,N,N))
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      SRI_bayes[,i,j] <- SRI_bayes[,j,i] <- rbeta(1000,sri_a[i,j],sri_b[i,j])
    }
  }
  return(SRI_bayes)
}