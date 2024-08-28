#### Load libraries and define working space  ####
# library(raster) 
library(foreach)
library(doParallel)
library(tictoc)
library(jagsUI)

Dir.Base <- getwd() ##Should be repository main directory
setwd(Dir.Base)

source(file.path(Dir.Base, "Utils",
                 "helper_functions.R"))


output.path <- file.path(Dir.Base, "Rdata")

if (!file.exists(output.path)){
  dir.create(file.path(output.path))
}


#### Load simulated data ####
#### Load landscape (psi1 covariable)
load(file.path(output.path, "landscape_1.Rdata"))

#### Load ini.occu and colext scenarios
load(file.path(output.path, "scenarios.Rdata"))

####Load detection data
##y.list (list of lists): detection array for each ini.occu scenario and colext scenario
##detection array (bin 0-1): sites, surveys, years, p.det.scenarios, replicates
load(file.path(output.path, paste0("detection_data.Rdata") )) ##y.list

####Load sampling data
##surveyed.cells.array (bool): sites, sampling.scenarios, replicates
load(file.path(output.path, "surveyed_cells_array.Rdata" ) )

#### Define variables ####
# ini.scn <- 2
# colext.scn <- 2

num.replicates <- dim(y.list$ini.range.1$high.inc)[5]
## Surveys parameters
nsurveys <- dim(y.list$ini.range.1$high.inc)[2]
y.array <- y.list[[ini.scn]][[colext.scn]]

prop.surveyed.scenarios.jags.i <- c(1,2,4,6,9) ##only analyse spatial coverage 5%, 10%, 20%, 40% and 80%


####  Specify model in BUGS language  ####
jags_model_file <- file.path(output.path, "jags_model.txt")
if (!file.exists(jags_model_file)){
  cat(file = jags_model_file,"
  model {
  # Priors
  psi1 ~ dunif(0, 1) # Initial occupancy
  alpha.lpsi1 <- logit(psi1)
  beta.lpsi1 ~ dnorm(0, 0.01)
  phi.int ~ dunif(0, 1) # Persistence
  alpha.lphi <- logit(phi.int)
  beta.lphi ~ dnorm(0, 0.01)
  gamma.int ~ dunif(0, 1) # Colonization
  alpha.lgamma <- logit(gamma.int)
  beta.lgamma ~ dnorm(0, 0.01)
  p ~ dunif(0, 1) # Detection
  # Likelihood
  # Ecological submodel
  for (i in 1:nsites){
  z[i,1] ~ dbern( ilogit(alpha.lpsi1 + beta.lpsi1 * field[i]) )
  for (t in 2:nyears){
  z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1] + (1-z[i,t-1]) * gamma[i,t-1])
  # Compute autocovariate and specify its effects on phi and gamma
  autocov[i,t-1] <- sum(z[neighID[i,1:numN[i]], t-1]) / numN[i]
  logit(phi[i,t-1]) <- alpha.lphi + beta.lphi * autocov[i,t-1]
  logit(gamma[i,t-1]) <- alpha.lgamma + beta.lgamma * autocov[i,t-1]
  }
  }
  # Observation model
  for (i in 1:nsites){
  for (j in 1:nsurveys){
  for (t in 1:nyears){
  y[i,j,t] ~ dbern(z[i,t] * p)
  }
  }
  }
  }
  ")
}


#### Fit JAGS SpDynOcc models ####
for (ini.scn in 1:length(ini.range.filling.scenarios)){
  for (colext.scn in 1:length(colext.scenarios)){
    
    y.array <- y.list[[ini.scn]][[colext.scn]]
    tic()
    
    for (nyears in nyears.scenarios){
      
      if (nyears!=12){
        prop.surveyed.scenarios.jags.i <- c(1,2,4,6,9) ##only analyse spatial coverage 5%, 10%, 20%, 40% and 80%
      } else{
        prop.surveyed.scenarios.jags.i <- c(1,2,4,6) ##only analyse spatial coverage 5%, 10%, 20% and 40%
      }
      print(paste0("Fitting all scenarios: ini scn: ",ini.scn, 
                   ", colext scn: ", colext.scn, " nyears: ", nyears))
      ##Parallel iteration over detection scenarios
      tic()
      n.cores <- 25
      #create the cluster
      my.cluster <- parallel::makeCluster(
        n.cores, 
        type = "FORK"
      )
      #register it to be used by %dopar%
      doParallel::registerDoParallel(cl = my.cluster)
      
      writeLines(c(""), "log_jags.txt")
      
      jags_model_fits <- 
        foreach(
          p.det.scn= 1:length(p.det.scenarios),
          .packages = c("jagsUI")
        )  %:%
        foreach(prop.surveyed.scn= prop.surveyed.scenarios.jags.i)  %:%
        foreach(i = 1:(num.replicates/2)) %dopar% { ##Only fit first 25 replicates

          cat(paste0("Fitting scenario pdet: ", p.det.scn, 
                     " prop: ", prop.surveyed.scenarios[prop.surveyed.scn],
                     " nyears: ", nyears,
                     " rep: ", i,"\n"), 
              file="log_jags.txt", append=TRUE)

          
          not_surveyed_cells <- !surveyed.cells.array[,prop.surveyed.scn, i]
          
          y_jags <- y.array[, ,1:nyears, p.det.scn, i]
          y_jags[not_surveyed_cells,,]<- NA
          
          
          bdata <- list(y = y_jags, nsites = nsites, nsurveys = nsurveys, nyears = nyears,
                        neighID = neighID, numN = numN,
                        field=as.vector(field) ) 
          
          
          # Initial values
          zst <- array(1, dim = c(nsites, nyears))
          inits <- function(){ list(z = zst)}
          # Parameters monitored
          params <- c("alpha.lpsi1", "beta.lpsi1", "alpha.lgamma", "beta.lgamma",
                      "alpha.lphi", "beta.lphi", "p") # could also monitor "z"
          # MCMC settings
          na <- 1000 ; ni <- 2000 ; nt <- 2 ; nb <- 1000 ; nc <- 3
          # Call JAGS (ART 30 min), check convergence and summarize posteriors
          
          out1 <- jags(bdata, inits, params, jags_model_file, n.adapt = na,
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       parallel = FALSE, verbose = FALSE)
          
        }
      parallel::stopCluster(cl = my.cluster)
      toc()
      
      save(jags_model_fits, file=file.path(output.path,
                                           paste0("jags_fits_ini_", ini.scn, "colext_", colext.scn, 
                                                  "_nyears_", nyears, ".Rdata")))
      # nyears.i <- which(nyears.scenarios==nyears)
      # jags_fits[[nyears.i]] <- jags_model_fits
    } ##for each study duration scenario
    toc()
  } ##for each ini.occu scenario
} ##for each colext scenario



#### Update models that did not converge 3000 iterations  ####


# ini.scn <- 1
# colext.scn <- 1

for (ini.scn in 1:length(ini.range.filling.scenarios)){
  for (colext.scn in 1:length(colext.scenarios)){
    for (nyears.i in nyears.scenarios){

     print(paste0("Updating all scenarios: ini scn: ",ini.scn, 
                   ", colext scn: ", colext.scn, " nyears: ", nyears.i))
      
      
      load(file.path(output.path, 
                     paste0("jags_fits_ini_", ini.scn, "colext_", colext.scn,
                            "_nyears_", nyears.i,".Rdata")) ) ##load fitted models
      
      ##Vector with all models
      unlisted_jags_fits <- unlist(jags_model_fits, recursive = F)
      unlisted_jags_fits <- unlist(unlisted_jags_fits, recursive = F)
      
      ##only select models that did not converge
      rhats <- sapply(unlisted_jags_fits, function(x) max(unlist(x$Rhat[-8])) )
      jags_fits_to_update <- unlisted_jags_fits[rhats>1.1]
      
      ##scenarios and replicate of models to update
      jags.scns <- expand.grid(i=1:(num.replicates/2), 
                               prop.surveyed.scn=1:length(prop.surveyed.scenarios.jags.i),
                               p.det.scn= 1:length(p.det.scenarios))
      jags.scns <- jags.scns[rhats>1.1,]
      
      ##update models
      print(paste0("Updating ", length(jags_fits_to_update) ," models"))
      ##parallelise mcmc update
      n.cores <- 50
      #create the cluster
      my.cluster <- parallel::makeCluster(
        n.cores, 
        type = "PSOCK"
      )
      #register it to be used by %dopar%
      doParallel::registerDoParallel(cl = my.cluster)
      
      writeLines(c(""), "log_jags.txt")
      tic()
      
      jags_fits_updated <- 
        foreach(
          i= 1:length(jags_fits_to_update),
          .packages = c("jagsUI"),
          .noexport = c("jags_model_fits", "unlisted_jags_fits")
        )  %dopar% {
          
          cat(paste0("Fitting ", i,"/",length(jags_fits_to_update) ,"\n"), 
              file="log_jags.txt", append=TRUE)
          
          out1 <- update(jags_fits_to_update[[i]], n.iter=3000)
        }
      
      parallel::stopCluster(cl = my.cluster)
      toc()
      
      
      ##Save updated models
      for (up.i in 1:length(jags_fits_updated)){
        i <- jags.scns$i[up.i]
        prop.surv.i <- jags.scns$prop.surveyed.scn[up.i]
        p.det.i <- jags.scns$p.det.scn[up.i]
        jags_model_fits[[p.det.i]][[prop.surv.i]][[i]] <- jags_fits_updated[[up.i]]
      }
      
      save(jags_model_fits, file=file.path(output.path, 
                                           paste0("jags_fits_ini_", ini.scn, "colext_", colext.scn,
                                                  "_nyears_", nyears.i, ".Rdata")))
      
      # rhats <- sapply(jags_fits_updated, function(x) max(unlist(x$Rhat[-8])) )
      # print(length(which(rhats<1.1)))
      
    } ##for each study duration scenario
  } ##for each ini.occu scenario
} ##for each colext scenario

#### Predict models using mean posteriors ####
nsim <- 1000
nyears <- 15

for (ini.scn in 1:length(ini.range.filling.scenarios)){
  for (colext.scn in 1:length(colext.scenarios)){
    print(paste0("Simulating occu models, ini scn: ",ini.scn, 
                 ", colext scn: ", colext.scn ))
    
    
    if (!file.exists(file=file.path(output.path,
                                    paste0("jags_number_occupied_sites_ini_",
                                           ini.scn, "colext_", colext.scn,".Rdata"))) ){
      jags.occu.sites <- list()
      for (nyears.i in nyears.scenarios){
        
        load(file.path(output.path,
                       paste0("jags_fits_ini_", ini.scn, "colext_", colext.scn,
                              "_nyears_", nyears.i,".Rdata")) )
        unlisted.model.fits <- unlist(jags_model_fits, recursive = F)
        unlisted.model.fits <- unlist(unlisted.model.fits, recursive = F)
        
        writeLines(c(""), "samples_log.txt")
        
        tic()
        
        n.cores <- parallel::detectCores()/2
        #create the cluster
        my.cluster <- parallel::makeCluster(
          n.cores, 
          type = "PSOCK"
        )
        #register it to be used by %dopar%
        doParallel::registerDoParallel(cl = my.cluster)
        
        jags.number.occupied.sites <- foreach(
          pred.i = 1:length(unlisted.model.fits),
          .packages = c("RcppSampleZ"),
          .noexport = "simulate_occurrence_cpp"
        ) %dopar% {
          cat(paste0("Sampling ", pred.i,"/", length(unlisted.model.fits), "\n"), 
              file="samples_log.txt", append=TRUE)
          pred <- unlisted.model.fits[[pred.i]]
          return(number.occupied.sites.year.bayesian.model(pred))
        }
        
        
        parallel::stopCluster(cl = my.cluster)
        toc()
        
        jags.occu.sites <- c(jags.occu.sites, jags.number.occupied.sites)
      } ## for each study duration scenario
      
      jags.number.occupied.sites <- jags.occu.sites
      save(jags.number.occupied.sites, file=
             file.path(output.path, paste0("jags_number_occupied_sites_ini_",
                                           ini.scn, "colext_", colext.scn,".Rdata")) )
      
    }
    
  } ##for each ini.occu scenario
} ##for each colext scenario

#### Update models that obtained bad fits and did not converge 5000 iterations ####
##Load true occu data
load(file.path(output.path, paste0("occu_zs.Rdata") )) ##zs.list


for (ini.scn in 1:length(ini.range.filling.scenarios)){
  
  for (colext.scn in 1:length(colext.scenarios)){
    
    #### Assess goodness of fit (predicted diff occu < max(true diff occu)) ####
    
    ##True diff occu among replicates
    true.occu <- zs.list[[ini.scn]][[colext.scn]]
    array.true.number.occupied.sites <- apply(true.occu, c(2,3), sum)
    array.true.diff.occu <- apply(array.true.number.occupied.sites, 2,
                                  occu_year_diff)
    mean.true.diff.occu <- apply(array.true.diff.occu, 2, mean)
    range.true <- max(mean.true.diff.occu) - min(mean.true.diff.occu)
    
    ##Load model predictions
    load(file.path(output.path, "jags_update",
                   paste0("jags_number_occupied_sites_ini_",
                          ini.scn, "colext_", colext.scn,".Rdata")))
    ##Scenarios analysed
    prop.surveyed.scenarios.jags <- c(1,2,4,6,9)
    jags.scns <- expand.grid(i=1:(num.replicates/2), 
                             prop.surveyed.scn=1:length(prop.surveyed.scenarios.jags),
                             p.det.scn= 1:length(p.det.scenarios),
                             nyears= c(3,6,12))
    jags.scns <- jags.scns[-which(jags.scns$prop.surveyed.scn==5 & jags.scns$nyears==12),] ##12 years scenarios do not include prop surveyed 9 (80%)
    
    ##Test goodness of fit each replicate
    for (jags.fit.i in 1:nrow(jags.scns)){
      i <- jags.scns$i[jags.fit.i]
      
      pred <- apply(jags.number.occupied.sites[[jags.fit.i]], 2, function(x) 
        occu_year_diff(x) )
      
      true <- array.true.diff.occu[,i]
      diff <- apply(pred, 2, function(x) x-true)
      
      ##bool: good fit & same tendency as true data (increasing or decreasing)
      jags.scns$correct.diff.occu[jags.fit.i] <- mean(pred)*mean(true) > 0 &
        (abs(mean(diff)) < range.true)
    }
    
    
    #### Assess model convergence ####
    jags.scns$rhats <- NA
    
    
    for (nyears.i in nyears.scenarios){
      load(file.path(output.path,
                     paste0("jags_fits_ini_", ini.scn, "colext_", colext.scn,
                            "_nyears_", nyears.i,".Rdata")) )
      rhats <- lapply(jags_model_fits, function(x) lapply(x, function(y)
        lapply(y, function(z) max(unlist(z$Rhat[-8])))))
      # print(median(unlist(rhats)))
      rhats_unlist <- unlist(rhats)
      # print(length(which(rhats_unlist<1.1)))
      jags.scns$rhats[jags.scns$nyears==nyears.i] <- rhats_unlist
    }
    
    
    jags.scns$colext <- colext.scn
    jags.scns$ini.occ <- ini.scn
    
    
    
    ####  Update models with bad predictions and no convergence ####
    for (nyears.i in nyears.scenarios){
      model.status <- jags.scns[jags.scns$nyears==nyears.i,] 
      
      load(file.path(output.path, 
                     paste0("jags_fits_ini_", ini.scn, "colext_", colext.scn,
                            "_nyears_", nyears.i,".Rdata")) )
      unlisted_jags_fits <- unlist(jags_model_fits, recursive = F)
      unlisted_jags_fits <- unlist(unlisted_jags_fits, recursive = F)
      
      
      jags_fits_to_update <- unlisted_jags_fits[
        which(!model.status$correct.diff.occu & model.status$rhats > 1.1)]
      
      model.status <- model.status[
        which(!model.status$correct.diff.occu & model.status$rhats > 1.1),]
      
      print(paste0("Updating all scenarios: ini scn: ",ini.scn, 
                   ", colext scn: ", colext.scn, " nyears: ", nyears.i))
      print(paste0("Updating ", length(jags_fits_to_update) ," models"))
      
      ##parallelise mcmc upload
      n.cores <- 30
      #create the cluster
      my.cluster <- parallel::makeCluster(
        n.cores, 
        type = "PSOCK"
      )
      #register it to be used by %dopar%
      doParallel::registerDoParallel(cl = my.cluster)
      
      writeLines(c(""), "log_jags.txt")
      tic()
      
      jags_fits_updated <- 
        foreach(
          i= 1:length(jags_fits_to_update),
          .packages = c("jagsUI"),
          .noexport = c("jags_model_fits", "unlisted_jags_fits")
        )  %dopar% {
          
          cat(paste0("Fitting ", i,"/",length(jags_fits_to_update) ,"\n"), 
              file="log_jags.txt", append=TRUE)
          
          
          
          out1 <- update(jags_fits_to_update[[i]], n.iter=10000)
          
        }
      
      parallel::stopCluster(cl = my.cluster)
      toc()
      
      ##Save updated models
      for (up.i in 1:length(jags_fits_updated)){
        i <- model.status$i[up.i]
        prop.surv.i <- model.status$prop.surveyed.scn[up.i]
        p.det.i <- model.status$p.det.scn[up.i]
        jags_model_fits[[p.det.i]][[prop.surv.i]][[i]] <- jags_fits_updated[[up.i]]
      }
      
      save(jags_model_fits, file=file.path(output.path, "jags_update",
                                           paste0("jags_fits_ini_", ini.scn, "colext_", colext.scn,
                                                  "_nyears_", nyears.i, ".Rdata")))
    } ##for each study duration scenario
  } ##for each ini.occu scenario
} ##for each colext scenario
