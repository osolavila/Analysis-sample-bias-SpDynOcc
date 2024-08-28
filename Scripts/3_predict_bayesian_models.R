#### Load libraries and define working space  ####
# library(raster) 
library(foreach)
library(doParallel)
library(tictoc)
# library(viridis)
# library(unmarked)
library(RcppSampleZ)
library(jagsUI)

Dir.Base <- getwd() ##Should be repository main directory
setwd(Dir.Base)

source(file.path(Dir.Base, "Utils",
                 "helper_functions.R"))


output.path <- file.path(Dir.Base, "Rdata")

if (!file.exists(output.path)){
  dir.create(file.path(output.path))
}

##Load scenarios definition
load(file.path(output.path, "scenarios.Rdata"))
## Load landscape
load(file.path(output.path, "landscape_1.Rdata"))

#### Predict occu using mean posteriors JAGS models ####
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
