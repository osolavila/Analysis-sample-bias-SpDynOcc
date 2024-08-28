library(raster)
library(gstat) #for variogram

Dir.Base <- getwd() ##Should be repository main directory
setwd(Dir.Base)

source(file.path(Dir.Base, "Utils",
                 "helper_functions.R"))


output.path <- file.path(Dir.Base, "Rdata")

if (!file.exists(output.path)){
  dir.create(file.path(output.path))
}

####Simulate occurrence scenarios

####  1. Landscape design  ####
if (file.exists(file.path(output.path, "landscape_1.Rdata"))){
  load(file.path(output.path, "landscape_1.Rdata"))
} else{
  
  side <- 50
  nsites <- side^2
 
  xcoord <- 1:side
  ycoord <- 1:side
  grid <- as.matrix(expand.grid(x=xcoord, y=ycoord))
  
  
  ##Raster landscape grids
  landscape_raster <- rasterFromXYZ(cbind(grid, 1:nrow(grid)))
  landscape_raster[] <- 1:length(landscape_raster)
 
  
  #### Get neighs information ####
  coord <- as.data.frame(raster::coordinates(landscape_raster)*1000) ##coordinates in m (we assume 1km grid)
  disp.par <- sqrt(2) + 0.1 ##neighborhood radius
  
  ##Adjacency list with neighbours indexes
  index_neighs <- list()
  
  for (i in 1:nrow(coord)){
    
    ##Neighbours distances
    points_distances <- pointDistance(coord[i ,],
                                      coord[,],
                                      lonlat=F)
    points_distances <- points_distances/1000
    ##Eliminate same cell distances
    points_distances[i] <- NA
    
    
    neigh <- which(points_distances<disp.par)
    
    
    index_neighs[[i]] <- neigh
    
  }
  
  ##Adjacency matrix for JAGS
  # Number of neighbors for each cell
  numN <- sapply(index_neighs, length)
  # Put the neighbor IDs into a matrix
  neighID <- array(NA, dim = c(nsites, 8))
  for(i in 1:nsites){
    neighID[i, 1:numN[i]] <- index_neighs[[i]]
  }
  
  ####  Simulate Gaussian random field (psi0 covariate)  ####

  ## define the gstat object (spatial model)
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025,model="Exp",range=2), nmax=20)
  
  ## make four simulations based on the stat object
  set.seed(21)
  gaussian_fields <- predict(g.dummy, newdata=coord/1000, nsim=4)
  
  ##Choose simulation
  raster::plot(rasterFromXYZ(gaussian_fields), col=topo.colors(20), box=FALSE) 
  field <- scale(gaussian_fields$sim1)
  
  #### Save landscape information ####
  save(side, nsites,  
       index_neighs, numN, neighID, field, 
       file=file.path(output.path, "landscape_1.Rdata"))
}

#### 2. Define and save scenarios of ini.occu, colext and p.det ####
ini.range.filling.scenarios <- c(0.01, 0.1, 0.9)
colext.scenarios <- c("high.inc",  "mod.inc",   "high.dec",  "mod.dec")
nyears.scenarios <- c(3,6,12)
prop.surveyed.scenarios <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
p.det.scenarios <- c(0.553,0.293,0.106) ##Single visit p.det: effective probability for 2 visits c(0.8, 0.5, 0.2)

save(ini.range.filling.scenarios,
     colext.scenarios,
     nyears.scenarios,
     prop.surveyed.scenarios,
     p.det.scenarios,
     file= file.path(output.path, "scenarios.Rdata") )

##Set and store model parameters values for each colext scenario
if (!(file.exists(file.path(output.path, paste0("occu_params.Rdata") )) )) {
  
  num.replicates <- 50
  ##Params
  nyears <- 15
  
  ##Calculate psi1 intercept so mean initial number of occupied sites corresponds to the ini.occu scenarios
  mean.psi1.vector <- c()
  beta.psi1 <- 3 ##Fixed beta.psi1 for all scenarios
  for (ini.range.scn in ini.range.filling.scenarios){
    mean.psi1 <- get_mean_psi1(ini.range.scn)
    mean.psi1.vector <- c(mean.psi1.vector, mean.psi1)
    
    ##Plot initial range filling distribution
    nsim <- 10000
    ini.occu <- vector("numeric", nsim)
    for (i in 1:nsim){
      psi0 <- plogis(mean.psi1+beta.psi1*(field))
      z <- rbinom(nsites, 1, psi0)
      ini.occu[i] <- length(which(z==1))/nsites
    }
    
    hist(ini.occu)
    print(mean(ini.occu))
  }
  
  ##Store colext and ini.occu params in list
  occu.params <- list()
  occu.params[["ini.range.1"]] <- list(
    mean.psi1=mean.psi1.vector[1],
    beta.psi1=beta.psi1,
    high.inc=list(mean.gamma= -8,
                  beta.gamma.auto= 24,
                  mean.phi=8,
                  beta.phi.auto=0),
    mod.inc=list(mean.gamma= -8,
                 beta.gamma.auto= 16,
                 mean.phi=8,
                 beta.phi.auto=0),
    high.dec=list(mean.gamma= -8,
                  beta.gamma.auto= 0,
                  mean.phi=0,
                  beta.phi.auto=10),
    mod.dec=list(mean.gamma= -8,
                 beta.gamma.auto= 0,
                 mean.phi=1.1,
                 beta.phi.auto=10)
  )
  
  occu.params[["ini.range.2"]] <- list(
    mean.psi1=mean.psi1.vector[2],
    beta.psi1=beta.psi1,
    high.inc=list(mean.gamma= -8,
                  beta.gamma.auto= 16,
                  mean.phi=8,
                  beta.phi.auto=0),
    mod.inc=list(mean.gamma= -8,
                 beta.gamma.auto= 10,
                 mean.phi=8,
                 beta.phi.auto=0),
    high.dec=list(mean.gamma= -8,
                  beta.gamma.auto= 0,
                  mean.phi=0,
                  beta.phi.auto=6),
    mod.dec=list(mean.gamma= -8,
                 beta.gamma.auto= 0,
                 mean.phi=2,
                 beta.phi.auto=6)
  )
  
  occu.params[["ini.range.3"]] <- list(
    mean.psi1=mean.psi1.vector[3],
    beta.psi1=beta.psi1,
    high.inc=list(mean.gamma= -8,
                  beta.gamma.auto= 10,
                  mean.phi=8,
                  beta.phi.auto=0),
    mod.inc=list(mean.gamma= -8,
                 beta.gamma.auto= 7,
                 mean.phi=8,
                 beta.phi.auto=0),
    high.dec=list(mean.gamma= -8,
                  beta.gamma.auto= 0,
                  mean.phi=0,
                  beta.phi.auto=4),
    mod.dec=list(mean.gamma= -8,
                 beta.gamma.auto= 0,
                 mean.phi=2,
                 beta.phi.auto=4)
  )
  
  ##Params
  save(
    occu.params,
    nyears,
    num.replicates,
    file= file.path(output.path, paste0("occu_params.Rdata") ))
  
} else{
  load(file.path(output.path, paste0("occu_params.Rdata") ))
}


####  3. Simulate Occurrence Scn  ####
##Simulate yearly occupancy data over 15 years for each scenario
if (!(file.exists(file.path(output.path, paste0("occu_zs.Rdata") ) ) ) ){

  zs.list <- list()
  set.seed(3)
  for (ini.range.scn in names(occu.params)){
    mean.psi1 <- occu.params[[ini.range.scn]]$mean.psi1
    beta.psi1 <- occu.params[[ini.range.scn]]$beta.psi1
    
    colext.sims <- list()
    for (colext.scn in colext.scenarios){
      
      mean.gamma <- occu.params[[ini.range.scn]][[colext.scn]]$mean.gamma
      beta.gamma.auto <- occu.params[[ini.range.scn]][[colext.scn]]$beta.gamma.auto
      mean.phi <- occu.params[[ini.range.scn]][[colext.scn]]$mean.phi
      beta.phi.auto <- occu.params[[ini.range.scn]][[colext.scn]]$beta.phi.auto
      
      zs <- array(numeric(), dim=c(nsites,nyears,num.replicates))
      for (i in 1:num.replicates){
        zs[,,i] <- simulate.occurrence(mean.psi1, beta.psi1,
                                       mean.gamma, beta.gamma.auto,
                                       mean.phi, beta.phi.auto)
      }
      colext.sims[[colext.scn]] <- zs
    }
    
    zs.list[[ini.range.scn]] <- colext.sims
    
  }
  
  # ##Plot trajectories (num occu/ years)
  # for (ini.range.scn in names(occu.params)){
  #   for (colext.scn in colext.scenarios){
  #     zs <- zs.list[[ini.range.scn]][[colext.scn]]
  #     plot(apply(zs, c(2,3), sum)[,1], type="l", 
  #          main=paste0(ini.range.scn, " ",colext.scn))
  #     for (i in 2:num.replicates){
  #       lines(apply(zs, c(2,3), sum)[,i])
  #     }
  #   }
  # }
  
  save(zs.list, 
       file=file.path(output.path, paste0("occu_zs.Rdata") ) )
} else{
  load(file.path(output.path, paste0("occu_zs.Rdata") )) ##zs.list
}

#### 4. Simulate Detection data  ####
if (!(file.exists(file.path(output.path, "detection_data.Rdata" ) ) ) ){
  nsurveys <- 2  
  set.seed(29)
  
  ##For each scenario get detection data:
  ##array (bin 0-1): sites, surveys, years, p.det.scenarios, replicates  
  y.list <- rapply(zs.list, sim.detection.data, how="replace")  
  
    
  ##Check  
  # a <- y.list$ini.range.3$high.inc
  # b <- a[,1,1,,]
  # apply(b, c(2,3), sum)
  
  save(y.list, 
         file=file.path(output.path, paste0("detection_data.Rdata") ) )
} else{
    load(file.path(output.path, paste0("detection_data.Rdata") )) ##y.list
}
  
####  5. Simulate Sampling design (spatial coverage)  ####
if (!(file.exists(file.path(output.path, "surveyed_cells_array.Rdata" )))){
  set.seed(14)
  ##array: sites,prop.surveyed.scenarios, replicates
  surveyed.cells.array <- array(NA, 
                                dim=c(nsites,length(prop.surveyed.scenarios),
                                      num.replicates) )
  
  for (replicate in 1:num.replicates){
    for (prop.surveyed.scn in 1:length(prop.surveyed.scenarios) ){
      prop.surveyed <- prop.surveyed.scenarios[prop.surveyed.scn]
      # sample_cells <- function(prop.surveyed)
      surveyed.cells <- sample(1:nsites, ceiling(nsites*prop.surveyed))
      surveyed.cells <- 1:nsites %in% surveyed.cells
      
      surveyed.cells.array[,prop.surveyed.scn,replicate] <- surveyed.cells
      
    }
  }
  save(surveyed.cells.array, 
       file=file.path(output.path, "surveyed_cells_array.Rdata" ) )
} else {
  load(file.path(output.path, "surveyed_cells_array.Rdata" ))
}

