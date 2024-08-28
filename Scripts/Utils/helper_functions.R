library(Rcpp)

##Optim function to get the intercept parameter of psi1 for a fixed initial range filling
get_mean_psi1 <- function(ini.range.scn, beta.psi1=3, nsim=10000,
                          interval= c(-9,5)){
  
  nll <- function(mean.psi1){
    ini.occu <- vector("numeric", nsim)
    for (i in 1:nsim){
      psi0 <- plogis(mean.psi1+beta.psi1*(field))
      z <- rbinom(nsites, 1, psi0)
      ini.occu[i] <- length(which(z==1))/nsites
    }
    return(abs(mean(ini.occu)-ini.range.scn))
  }
  
  fm <- optimise(nll, interval=interval)
  return(fm$minimum)
}

##simulate autologistic spatial dynamic occupancy model
##must predefine adjacency list named index_neighs
simulate.occurrence <- function(mean.psi1, beta.psi1,
                                mean.gamma, beta.gamma.auto,
                                mean.phi, beta.phi.auto){  
  
  
  psi0 <- plogis(mean.psi1+beta.psi1*(field))
  z <- array(dim = c(nsites, nyears)) # Occupancy, occurrence
  z[,1] <- rbinom(nsites, 1, psi0)
  
  
  #### Ocuppancy at ti  ####
  
  for (t in 2:nyears){
    ##Connectivity
    sp_autocovariate <- sapply(index_neighs,
                               function(x) sum(z[x,t-1], na.rm=T)/length(x) ) 
    ##Col probability
    prob.col <- plogis(mean.gamma+beta.gamma.auto*sp_autocovariate)
    
    ##Surv probability
    prob.surv <- plogis(mean.phi+beta.phi.auto*sp_autocovariate)
    
    ##Update occu
    z[,t] <- rbinom(nsites, 1, z[,t-1]*prob.surv+(1-z[,t-1])*prob.col)
    
    
  }
  return(z)
  
  
}


##Simulate detection data for all the probability of detectio scenarios
sim.detection.data <- function(occu.data, .p.det.scenarios=p.det.scenarios, 
                               .nsurveys=nsurveys){
  ##array (bin 0-1): sites, surveys, years, p.det.scenarios, replicates
  y.array <- array(NA, dim=c(nsites, .nsurveys, nyears,
                             length(p.det.scenarios), num.replicates))
  for (replicate in 1:num.replicates){
    z <- occu.data[,,replicate]
    for (p.det.scn in 1:length(.p.det.scenarios)){
      p.det <- .p.det.scenarios[p.det.scn]
      for (s in 1:.nsurveys){
        y.array[,s,,p.det.scn,replicate] <- rbinom(nsites*nyears, 1, z*p.det)
      }
    }
  }
  return(y.array)
}

##Calculate change in total number of occupied sites t+1 - t
occu_year_diff <- function(occu.vector){
  l <- length(occu.vector)
  diff.vector <- occu.vector[2:l]-occu.vector[1:(l-1)]
  diff.vector
}

##Extracts model mean posterior estimates and calls number.occupied.sites.year.cpp
number.occupied.sites.year.bayesian.model <- function(model.fit, cpp=TRUE){
  
  ##Occu t0 
  mean.psi1 <- model.fit$mean[["alpha.lpsi1"]]
  beta.psi1 <- model.fit$mean[["beta.lpsi1"]]
  
  ##Gamma
  mean.gamma <- model.fit$mean[["alpha.lgamma"]]
  beta.gamma.auto <- model.fit$mean[["beta.lgamma"]]
  
  ##Phi
  mean.phi <- model.fit$mean[["alpha.lphi" ]]
  beta.phi.auto <- model.fit$mean[["beta.lphi"]]
  
  if (cpp){
    
    return(number.occupied.sites.year.cpp(mean.psi1, beta.psi1,
                                          mean.gamma, beta.gamma.auto,
                                          mean.phi, beta.phi.auto) )
  }
  else {
    print("no R function for this yet")
  }
  
}

##Predict total number of occupied sites per year of nsim (variable to be predifined) simulations
##must predefine "numN" vector (length=nsites) w nº neighs per site and a "neighID" matrix (nrow=nsites, ncol=max nº neighs) w neighs ID
number.occupied.sites.year.cpp <- function(mean.psi1, beta.psi1,
                                           mean.gamma, beta.gamma.auto,
                                           mean.phi, beta.phi.auto){  
  num_zs <- array(0L, dim=c(nyears, nsim))
  for (i in 1:nsim){
    num_zs[,i] <- apply(simulate_occurrence_cpp(numN, neighID-1, nyears, nsites, 
                                                as.numeric(field),
                                                c(mean.psi1, beta.psi1),
                                                c(mean.phi, beta.phi.auto), 
                                                c(mean.gamma, beta.gamma.auto) )
                        , 2,
                        sum)
    
  }
  
  
  return(num_zs)
  
}

##Install package containing simulate_occurrence_cpp
if(!require("RcppSampleZ")){
  
  Rcpp.package.skeleton(name= "RcppSampleZ",
                        cpp_files="./Utils/sample_z.cpp",
                        example_code=FALSE) 
  
  install.packages("./RcppSampleZ", repos = NULL, type="source")
}