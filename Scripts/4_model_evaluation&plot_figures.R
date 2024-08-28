#### Load libraries and define working space  ####
# library(raster) 

library(tictoc)
library(viridis)
library(dplyr)
library(cowplot)
library(gridExtra)
library(grid)
library(lme4)
# library(jagsUI)

Dir.Base <- getwd() ##Should be repository main directory
setwd(Dir.Base)

source(file.path(Dir.Base, "Utils",
                 "helper_functions.R"))


output.path <- file.path(Dir.Base, "Rdata")

if (!file.exists(output.path)){
  dir.create(file.path(output.path))
}

plots.path <- file.path(Dir.Base, "Plots")
if (!file.exists(plots.path)){
  dir.create(file.path(plots.path))
}

##Load scenarios definition
load(file.path(output.path, "scenarios.Rdata"))
##Load real occu parameters values including nyears and num.replicates
load(file.path(output.path, "occu_params.Rdata" ))

##Load true occu data
load(file.path(output.path, paste0("occu_zs.Rdata") )) ##zs.list
load(file.path(output.path, "landscape_1.Rdata"))

nyears <- 15
ini.scn <- 2
colext.scn <- 1

true.occu <- zs.list[[ini.scn]][[colext.scn]]

#### Pred vs Obs occu change 10% ####

plots.list <- list()
# ini.scn <- 2
# colext.scn <- 1
for (ini.scn in 1:length(ini.range.filling.scenarios)){
  plots.sublist <- list()
  for (colext.scn in 1:length(colext.scenarios)){
    ##Observed occupancy change
    true.occu <- zs.list[[ini.scn]][[colext.scn]]
    array.true.number.occupied.sites <- apply(true.occu, c(2,3), sum)
    array.true.diff.occu <- apply(array.true.number.occupied.sites, 2,
                                  occu_year_diff)
    mean.true.diff <- apply(array.true.diff.occu, 2, mean)
    
    ##Predicted occupancy change
    load(file.path(output.path, 
                   paste0("jags_number_occupied_sites_ini_",
                          ini.scn, "colext_", colext.scn,".Rdata")))
    prop.surveyed.scenarios.jags <- c(1,2,4,6,9)
    jags.scns <- expand.grid(i=1:(num.replicates/2), 
                             prop.surveyed.scn=prop.surveyed.scenarios.jags,
                             p.det.scn= 1:length(p.det.scenarios),
                             nyears= c(3,6,12))
    jags.scns <- jags.scns[-which(jags.scns$prop.surveyed.scn==9 & jags.scns$nyears==12),] ##12 years scenarios do not include prop surveyed 9 (80%)
    
    jags.scns$diff.occu <- NA
    
    for (jags.fit.i in 1:nrow(jags.scns)){
      i <- jags.scns$i[jags.fit.i]
      
      
      pred <- apply(jags.number.occupied.sites[[jags.fit.i]], 2, function(x) 
        occu_year_diff(x) )
      
      mean.pred <- mean(apply(pred, 2, mean))
      
      jags.scns$diff.occu[jags.fit.i] <- mean.pred
    }
    
    
    ##Scenarios as factors for plotting
    jags.scns[,1:4] <- lapply(jags.scns[,1:4], factor)
    levels(jags.scns$nyears) <- paste0(nyears.scenarios, "y")
    levels(jags.scns$prop.surveyed.scn) <- paste0(prop.surveyed.scenarios[
      prop.surveyed.scenarios.jags]*100, "%")
    levels(jags.scns$p.det.scn) <- c(0.8, 0.5, 0.2)
    jags.scns$p.det.scn <- factor(jags.scns$p.det.scn, levels=c(0.2,0.5,0.8))
    
    ##Select spatial coverage 10% only
    jags.scns <- filter(jags.scns, prop.surveyed.scn=="10%")
    
    ##Plot predicted occupancy change
    p<-ggplot(jags.scns, aes(x=p.det.scn, y=diff.occu,
                             fill=nyears)) +
      geom_boxplot(position="dodge2", outlier.size = 0.7) +
      # facet_wrap(~prop.surveyed.scn, nrow= 1) +
      # labs(y = value.name) + 
      scale_fill_manual(values=viridis(length(nyears.scenarios))) +
      theme_bw() +
      theme(plot.margin = unit(c(0,0.25,0.25,0), "cm"),
            legend.position = "none",
            axis.title.x = element_blank(),  
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 7),  
            axis.text.y = element_text(size = 7))   
    
    
    ##Plot observed max and min occupancy change
    p <- p + geom_hline(yintercept=max(mean.true.diff)) +
      geom_hline(yintercept=min(mean.true.diff)) 
    
    
    plots.sublist[[colext.scn]] <- p 
    
  } ##for each colext scenario
  
  plots.list[[ini.scn]] <- plots.sublist
} ##for each ini.occu scenario


#### Figure 1 - Plot boxplots occu change 10%  ####
plots <- plot_grid(
  plots.list[[1]][[1]], plots.list[[2]][[1]], plots.list[[3]][[1]],
  plots.list[[1]][[2]],plots.list[[2]][[2]], plots.list[[3]][[2]],
  plots.list[[1]][[3]], plots.list[[2]][[3]], plots.list[[3]][[3]],
  plots.list[[1]][[4]], plots.list[[2]][[4]], plots.list[[3]][[4]],
  ncol = 3 
  # labels = c("A1", "", "A3","A4","B1", "B2", "B3","B4","C1", "C2", "C3", "C4"),
  # hjust=0.08, vjust=1
) + # Add some space around the edges
  theme(plot.margin = unit(c(0,0,-0.3,0), "cm"))



##Add columns and row titles
col1 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "ini occu 1%") + 
  theme_void()
col2 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "ini occu 10%") + 
  theme_void()
col3 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "ini occu 90%") + 
  theme_void()

col.titles <- plot_grid(
  col1,col2,col3,
  ncol = 3 
) + # Add some space around the edges
  theme(plot.margin = unit(c(0,0,0,0), "cm"))



row1 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "HI") + 
  theme_void()
row2 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "MI") + 
  theme_void()
row3 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "HD") + 
  theme_void()
row4 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "MD") + 
  theme_void()

row.titles <- plot_grid(
  row1,row2,row3, row4,
  ncol = 1 
) + # Add some space around the edges
  theme(plot.margin = unit(c(0,0,0,-0.1), "cm"))

void.plot <- ggplot() + 
  theme_void()

plots.titles <- plot_grid(col.titles, void.plot, plots, row.titles,
                          rel_heights = c(1,15), rel_widths = c(18,1),
                          ncol=2)



y.grob <- textGrob("Mean Yearly Change Number of Occupied Cells (predicted vs true)", 
                   gp=gpar(  fontsize=11), rot=90)

x.grob <- textGrob("Probability of detection", 
                   gp=gpar(fontsize=11))

plots <- grid.arrange(arrangeGrob(plots.titles, left = y.grob, bottom = x.grob))

## Add legend
p<-ggplot(jags.scns, aes(x=p.det.scn, y=diff.occu,
                         fill=nyears)) +
  geom_boxplot(position="dodge2") +   
  scale_fill_manual(values=viridis(length(nyears.scenarios)), 
                    name = "study duration") +
  theme(plot.margin = unit(c(0,0.25,0.25,0), "cm"),
        legend.position = "bottom", legend.box = "horizontal",
        axis.title.x = element_text(size = 10),  # Adjust the size as needed
        axis.title.y = element_blank(),
        legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=10),
        legend.text = element_text(size=8))   # Adjust the size as needed)



legend <- get_legend(p)

full.plot <- plot_grid(plots, legend, ncol=1,
                       rel_heights = c(0.92,0.08)) 

ggsave(file.path(plots.path, "figure1_boxplots_occu_change.jpeg"), plot=full.plot, 
       width= 18, height=22, units="cm",  bg = "transparent")
#### Pred vs Obs occu change ####
plots.list <- list()
# ini.scn <- 2
# colext.scn <- 1
for (ini.scn in 1:length(ini.range.filling.scenarios)){
  plots.sublist <- list()
  for (colext.scn in 1:length(colext.scenarios)){
    
    true.occu <- zs.list[[ini.scn]][[colext.scn]]
    array.true.number.occupied.sites <- apply(true.occu, c(2,3), sum)
    array.true.diff.occu <- apply(array.true.number.occupied.sites, 2,
                                  occu_year_diff)
    mean.true.diff <- apply(array.true.diff.occu, 2, mean)
    
    ##Estimated diff occu
    load(file.path(output.path, 
                   paste0("jags_number_occupied_sites_ini_",
                          ini.scn, "colext_", colext.scn,".Rdata")))
    prop.surveyed.scenarios.jags <- c(1,2,4,6,9)
    jags.scns <- expand.grid(i=1:(num.replicates/2), 
                             prop.surveyed.scn=prop.surveyed.scenarios.jags,
                             p.det.scn= 1:length(p.det.scenarios),
                             nyears= c(3,6,12))
    jags.scns <- jags.scns[-which(jags.scns$prop.surveyed.scn==9 & jags.scns$nyears==12),] ##12 years scenarios do not include prop surveyed 9 (80%)
    
    jags.scns$diff.occu <- NA
    
    for (jags.fit.i in 1:nrow(jags.scns)){
      i <- jags.scns$i[jags.fit.i]
      
      
      pred <- apply(jags.number.occupied.sites[[jags.fit.i]], 2, function(x) 
        occu_year_diff(x) )
      
      mean.pred <- mean(apply(pred, 2, mean))
      
      jags.scns$diff.occu[jags.fit.i] <- mean.pred
    }
    
    
    
    jags.scns[,1:4] <- lapply(jags.scns[,1:4], factor)
    levels(jags.scns$nyears) <- paste0(nyears.scenarios, "y")
    levels(jags.scns$prop.surveyed.scn) <- paste0("spatial coverage ",
                                                  prop.surveyed.scenarios[
                                                    prop.surveyed.scenarios.jags]*100, "%")
    levels(jags.scns$p.det.scn) <- c(0.8, 0.5, 0.2)
    jags.scns$p.det.scn <- factor(jags.scns$p.det.scn, levels=c(0.2,0.5,0.8))
    
    
    
    p<-ggplot(jags.scns, aes(x=p.det.scn, y=diff.occu,
                             fill=nyears)) +
      geom_boxplot(position="dodge2", outlier.size = 0.7) +
      facet_wrap(~prop.surveyed.scn, nrow= 1) +
      # labs(y = value.name) + 
      scale_fill_manual(values=viridis(length(nyears.scenarios))) +
      theme_bw() +
      theme(plot.margin = unit(c(0,0.25,0.25,0), "cm"),
            legend.position = "none",
            axis.title.x = element_blank(),  
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 7),  
            axis.text.y = element_text(size = 7))   
    
    
    
    p <- p + geom_hline(yintercept=max(mean.true.diff)) +
      geom_hline(yintercept=min(mean.true.diff)) +
      labs(title=paste0("ini occu ", c("1%", "10%", "90%")[ini.scn] ))
    
    # print(p)
    
    plots.sublist[[colext.scn]] <- p 
    
  }
  
  plots.list[[ini.scn]] <- plots.sublist
}


#### Supplementary Figure - Plot boxplots occu change   ####

#specify path to save PDF to
destination <- file.path(plots.path, "figureS4_boxplots_occu_change.pdf")
#open PDF

pdf(file=destination, width = 8.3, height = 7)

scn <- 1
titles <- paste0("Colext scenario: ", c("high occupancy increase",
                                        "moderate occupancy increase",
                                        "high occupancy decrease",
                                        "moderate occupancy decrease"))

for (scn in 1:length(colext.scenarios)){
  plots <- plot_grid(
    plots.list[[1]][[scn]], plots.list[[2]][[scn]], 
    plots.list[[3]][[scn]],
    ncol = 1
  ) + # Add some space around the edges  
    theme(plot.margin = unit(c(0,0,-0.10,0.15), "cm")) 
  
  
  y.grob <- textGrob("Mean Yearly Change Number of Occupied Cells (predicted vs true)", 
                     gp=gpar(  fontsize=11), rot=90)
  
  x.grob <- textGrob("Probability of detection", 
                     gp=gpar(fontsize=11))
  
  p <- arrangeGrob(plots, left = y.grob, bottom = x.grob)
  
  title <- ggdraw() + draw_label(titles[scn],
                                 fontface='bold')
  p <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
  
  
  
  p.leg <-ggplot(jags.scns, aes(x=p.det.scn, y=diff.occu,
                                fill=nyears)) +
    geom_boxplot(position="dodge2") +   
    scale_fill_manual(values=viridis(length(nyears.scenarios)), 
                      name = "study duration") +
    theme_bw() +
    theme(plot.margin = unit(c(0,0.25,0.25,0), "cm"),
          # legend.position = "bottom", legend.box = "vertical",
          axis.title.x = element_text(size = 10),  # Adjust the size as needed
          axis.title.y = element_blank(),
          legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8),
          legend.text = element_text(size=7))   # Adjust the size as needed)
  
  
  
  legend <- get_legend(p.leg)
  
  full.plot <- plot_grid(p, legend, ncol=2,
                         rel_widths = c(0.90,0.10)) 
  
  print(full.plot)
}

dev.off()

#### Pred vs Obs occu 1y ####
plots.list <- list()
ini.scn <- 2
colext.scn <- 1
for (ini.scn in 1:length(ini.range.filling.scenarios)){
  plots.sublist <- list()
  for (colext.scn in 1:length(colext.scenarios)){
    
    true.occu <- zs.list[[ini.scn]][[colext.scn]]
    array.true.number.occupied.sites <- apply(true.occu, c(2,3), sum)
    true.number.occupied.sites.y1 <- array.true.number.occupied.sites[1,]
    
    ##Estimated diff occu
    load(file.path(output.path, 
                   paste0("jags_number_occupied_sites_ini_",
                          ini.scn, "colext_", colext.scn,".Rdata")))
    prop.surveyed.scenarios.jags <- c(1,2,4,6,9)
    jags.scns <- expand.grid(i=1:(num.replicates/2), 
                             prop.surveyed.scn=prop.surveyed.scenarios.jags,
                             p.det.scn= 1:length(p.det.scenarios),
                             nyears= c(3,6,12))
    jags.scns <- jags.scns[-which(jags.scns$prop.surveyed.scn==9 & jags.scns$nyears==12),] ##12 years scenarios do not include prop surveyed 9 (80%)
    
    jags.scns$num.occu <- NA
    
    for (jags.fit.i in 1:nrow(jags.scns)){
      i <- jags.scns$i[jags.fit.i]
      
      
      pred <- jags.number.occupied.sites[[jags.fit.i]][1,]
      
      mean.pred <- mean(pred)
      
      jags.scns$num.occu[jags.fit.i] <- mean.pred
    }
    
    
    
    jags.scns[,1:4] <- lapply(jags.scns[,1:4], factor)
    levels(jags.scns$nyears) <- paste0(nyears.scenarios, "y")
    levels(jags.scns$prop.surveyed.scn) <- paste0("spatial coverage ",
                                                  prop.surveyed.scenarios[
                                                    prop.surveyed.scenarios.jags]*100, "%")
    levels(jags.scns$p.det.scn) <- c(0.8, 0.5, 0.2)
    jags.scns$p.det.scn <- factor(jags.scns$p.det.scn, levels=c(0.2,0.5,0.8))
    
    
    
    p<-ggplot(jags.scns, aes(x=p.det.scn, y=num.occu,
                             fill=nyears)) +
      geom_boxplot(position="dodge2", outlier.size = 0.7) +
      facet_wrap(~prop.surveyed.scn, nrow= 1) +
      # labs(y = value.name) + 
      scale_fill_manual(values=viridis(length(nyears.scenarios))) +
      theme_bw() +
      theme(plot.margin = unit(c(0,0.25,0.25,0), "cm"),
            legend.position = "none",
            axis.title.x = element_blank(),  
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 7),  
            axis.text.y = element_text(size = 7))   
    
    
    
    p <- p + geom_hline(yintercept=max(true.number.occupied.sites.y1)) +
      geom_hline(yintercept=min(true.number.occupied.sites.y1)) +
      labs(title=paste0("ini occu ", c("1%", "10%", "90%")[ini.scn] ))
    
    # print(p)
    
    plots.sublist[[colext.scn]] <- p 
    
  }
  
  plots.list[[ini.scn]] <- plots.sublist
}


#### Supplementary Figure - Plot boxplots occu 1y   ####

#specify path to save PDF to
destination <- file.path(plots.path, "figureS5_boxplots_occu_y1.pdf")
#open PDF

pdf(file=destination, width = 8.3, height = 7)

scn <- 1
titles <- paste0("Colext scenario: ", c("high occupancy increase",
                                        "moderate occupancy increase",
                                        "high occupancy decrease",
                                        "moderate occupancy decrease"))

for (scn in 1:length(colext.scenarios)){
  plots <- plot_grid(
    plots.list[[1]][[scn]], plots.list[[2]][[scn]], 
    plots.list[[3]][[scn]],
    ncol = 1
  ) + # Add some space around the edges  
    theme(plot.margin = unit(c(0,0,-0.10,0.15), "cm")) 
  
  
  y.grob <- textGrob("Initial Number of Occupied Cells (predicted vs true)", 
                     gp=gpar(  fontsize=11), rot=90)
  
  x.grob <- textGrob("Probability of detection", 
                     gp=gpar(fontsize=11))
  
  p <- arrangeGrob(plots, left = y.grob, bottom = x.grob)
  
  title <- ggdraw() + draw_label(titles[scn],
                                 fontface='bold')
  p <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
  
  
  
  p.leg <-ggplot(jags.scns, aes(x=p.det.scn, y=num.occu,
                                fill=nyears)) +
    geom_boxplot(position="dodge2") +   
    scale_fill_manual(values=viridis(length(nyears.scenarios)), 
                      name = "study duration") +
    theme_bw() +
    theme(plot.margin = unit(c(0,0.25,0.25,0), "cm"),
          # legend.position = "bottom", legend.box = "vertical",
          axis.title.x = element_text(size = 10),  # Adjust the size as needed
          axis.title.y = element_blank(),
          legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size=8),
          legend.text = element_text(size=7))   # Adjust the size as needed)
  
  
  
  legend <- get_legend(p.leg)
  
  full.plot <- plot_grid(p, legend, ncol=2,
                         rel_widths = c(0.90,0.10)) 
  
  print(full.plot)
}

dev.off()

#### Calculate number of good fits each scenario (in terms of occu change) ####
# ini.scn <- 2
# colext.scn <- 1
plots.list <- list()
data.list <- list()

for (ini.scn in 1:length(ini.range.filling.scenarios)){
  plots.sublist <- list()
  data.sublist <- list()
  for (colext.scn in 1:length(colext.scenarios)){
    
    ##True diff occu
    true.occu <- zs.list[[ini.scn]][[colext.scn]]
    array.true.number.occupied.sites <- apply(true.occu, c(2,3), sum)
    array.true.diff.occu <- apply(array.true.number.occupied.sites, 2,
                                  occu_year_diff)
    mean.true.diff.occu <- apply(array.true.diff.occu, 2, mean)
    range.true <- max(mean.true.diff.occu) - min(mean.true.diff.occu)
    
    ##Estimated diff occu
    load(file.path(output.path, "jags_update_5000",
                   paste0("jags_number_occupied_sites_ini_",
                          ini.scn, "colext_", colext.scn,".Rdata")))
    prop.surveyed.scenarios.jags <- c(1,2,4,6,9)
    jags.scns <- expand.grid(i=1:(num.replicates/2), 
                             prop.surveyed.scn=prop.surveyed.scenarios.jags,
                             p.det.scn= 1:length(p.det.scenarios),
                             nyears= c(3,6,12))
    jags.scns <- jags.scns[-which(jags.scns$prop.surveyed.scn==9 & jags.scns$nyears==12),] ##12 years scenarios do not include prop surveyed 9 (80%)
    
    
    # jags.scns$diff.occu <- NA
    
    for (jags.fit.i in 1:nrow(jags.scns)){
      i <- jags.scns$i[jags.fit.i]
      
      
      pred <- apply(jags.number.occupied.sites[[jags.fit.i]], 2, function(x) 
        occu_year_diff(x) )
      
      true <- array.true.diff.occu[,i]
      
      diff <- apply(pred, 2, function(x) x-true)
      
      # jags.scns$diff.occu[jags.fit.i] <- mean(diff)
      
      ##bool: good fit & same tendency as true data (increasing or decreasing)
      jags.scns$correct.diff.occu[jags.fit.i] <- mean(pred)*mean(true) > 0 &
        (abs(mean(diff)) < range.true)
    }
    
    jags.scns[,1:4] <- lapply(jags.scns[,1:4], factor)
    levels(jags.scns$nyears) <- paste0(nyears.scenarios, "y")
    levels(jags.scns$prop.surveyed.scn) <- paste0(prop.surveyed.scenarios[
      prop.surveyed.scenarios.jags]*100, "%")
    levels(jags.scns$p.det.scn) <- paste0("p.det ",c(0.8, 0.5, 0.2))
    # jags.scns$p.det.scn <- factor(jags.scns$p.det.scn, levels=c(0.8,0.5,0.2))
    
    jags.scns$ini.occu <- ini.scn
    jags.scns$colext <- colext.scn
    
    ##Calculate total number of good fits among replicates of each scenario
    jags.mean.scns <- jags.scns %>%
      group_by(nyears, prop.surveyed.scn, p.det.scn) %>%
      summarise(correct.diff.occu = sum(correct.diff.occu))
    
    
    p <- ggplot(jags.mean.scns, aes(prop.surveyed.scn, nyears, fill= correct.diff.occu)) + 
      geom_tile() +
      facet_wrap(~ p.det.scn, nrow= 3) +
      scale_fill_gradient(low="white", high="green") +
      labs(x= "% surveyed cells",
           y = "num. years") +
      geom_text(aes(label = correct.diff.occu), color = "black", size = 2) +
      theme(plot.margin = unit(c(0,0.25,0.25,0), "cm"),
            strip.text = element_text(size = 8, margin = margin(b=1.5)),
            panel.spacing = unit(0.5, "lines"),
            legend.position = "none",
            axis.title.x = element_blank(),  
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 7),  
            axis.text.y = element_text(size = 7))   
    
    plots.sublist[[colext.scn]] <- p 
    data.sublist[[colext.scn]] <- jags.scns 
  }
  plots.list[[ini.scn]] <- plots.sublist
  data.list[[ini.scn]] <- data.sublist
}
# dev.off()


#### Predict model perfornace ####
for (ini.scn in 1:length(ini.range.filling.scenarios)){
  data.list[[ini.scn]] <- do.call(rbind, data.list[[ini.scn]])
}
jags.good.fits <- do.call(rbind, data.list)


jags.good.fits$ini.occu <- factor(jags.good.fits$ini.occu)
jags.good.fits$colext <- factor(jags.good.fits$colext)

##Predict probability of obtaining a good fit
##Study duration, spatial coverage and p.det fixed effects
##ini.occu scenario interaction with colext scenario as random intercept
# mylogit <- glm(correct.diff.occu ~ 0 + nyears + p.det.scn + prop.surveyed.scn, data = jags.good.fits, family = "binomial")
mylogit2 <- glmer(correct.diff.occu ~ 0 + nyears + p.det.scn + prop.surveyed.scn + (1| ini.occu:colext), data = jags.good.fits, family = "binomial")

##predict only fixed effects
jags.scns <- expand.grid(prop.surveyed.scn=prop.surveyed.scenarios.jags,
                         p.det.scn= 1:length(p.det.scenarios),
                         nyears= c(3,6,12))
jags.scns[,] <- lapply(jags.scns[,], factor)
levels(jags.scns$nyears) <- paste0(nyears.scenarios, "y")
levels(jags.scns$prop.surveyed.scn) <- paste0(prop.surveyed.scenarios[
  prop.surveyed.scenarios.jags]*100, "%")
levels(jags.scns$p.det.scn) <- paste0("p.det ",c(0.8, 0.5, 0.2))

jags.scns$pred <- predict(mylogit2, newdata=jags.scns, type = "response", re.form=NA)
#### Figure 2 - predicted model  performance  ####

##Plot fixed effects probability of obtaining good fits
p1 <- ggplot(jags.scns, aes(prop.surveyed.scn, nyears, fill= pred)) + 
  geom_tile() +
  facet_wrap(~ p.det.scn, nrow= 3) +
  scale_fill_gradient(low="white", high="green", limits = c(0, 1), 
                      breaks = seq(0,1, by=0.25)) +
  labs(x= "% surveyed cells",
       y = "num. years") +
  geom_text(aes(label = round(pred,2)), color = "black", size = 3) +
  xlab("survey spatial coverage") + ylab("study duration") +
  # scale_colour_gradient(limits = c(0, 1), 
  #                       breaks = seq(0,1, by=0.25)) +
  theme(plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position = "bottom", legend.box = "horizontal", 
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        # legend.position = "none",
        # axis.title.x = element_blank(),  
        # axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7))   


##Plot random effects probability of obtainig good fits (logit scale)
occu.scns <- expand.grid(colext.scn=colext.scenarios,
                         ini.occu.scn= factor(ini.range.filling.scenarios))

occu.scns$r.effect <- mylogit2@u  

levels(occu.scns$colext.scn) <- c("HI", "MI", "HD", "MD")
levels(occu.scns$ini.occu.scn) <- paste0(c(1,10,90), "%")


p2 <- ggplot(occu.scns, aes(ini.occu.scn, colext.scn, fill= r.effect)) + 
  geom_tile() +
  scale_fill_gradient() +
  labs(x= "ini occu scn" ,
       y = "colext scn") +
  geom_text(aes(label = round(r.effect,2)), color = "black", size = 3) +
  theme(legend.position = "bottom", legend.box = "horizontal", 
        legend.title = element_blank(),
        legend.text = element_text(size=8)) +
  scale_y_discrete(limits=rev, position = "right")


plots <- plot_grid(p1, p2, rel_widths = c(3,2), ncol = 2)

ggsave(file.path(plots.path, "fig_2_pred_model_performance.jpeg"), plot=plots, 
       width= 18, height=12, units="cm")
#### Figure 3 - observed number good fits per scenario ####

plots <- plot_grid(
  plots.list[[1]][[1]], plots.list[[2]][[1]], plots.list[[3]][[1]],
  plots.list[[1]][[2]],plots.list[[2]][[2]], plots.list[[3]][[2]],
  plots.list[[1]][[3]], plots.list[[2]][[3]], plots.list[[3]][[3]],
  plots.list[[1]][[4]], plots.list[[2]][[4]], plots.list[[3]][[4]],
  ncol = 3 
  # labels = c("A1", "", "A3","A4","B1", "B2", "B3","B4","C1", "C2", "C3", "C4"),
  # hjust=0.08, vjust=1
) + # Add some space around the edges
  theme(plot.margin = unit(c(0,0,-0.3,0), "cm"))



##Add columns and row titles
col1 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "ini occu 1%") + 
  theme_void()
col2 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "ini occu 10%") + 
  theme_void()
col3 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "ini occu 90%") + 
  theme_void()

col.titles <- plot_grid(
  col1,col2,col3,
  ncol = 3 
) + # Add some space around the edges
  theme(plot.margin = unit(c(0,0,0,0), "cm"))



row1 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "HI") + 
  theme_void()
row2 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "MI") + 
  theme_void()
row3 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "HD") + 
  theme_void()
row4 <- ggplot() + 
  annotate("text", x = 4, y = 10, size=4, label = "MD") + 
  theme_void()

row.titles <- plot_grid(
  row1,row2,row3, row4,
  ncol = 1 
) + # Add some space around the edges
  theme(plot.margin = unit(c(0,0,0,-0.1), "cm"))

void.plot <- ggplot() + 
  theme_void()

plots.titles <- plot_grid(col.titles, void.plot, plots, row.titles,
                          rel_heights = c(1,15), rel_widths = c(18,1),
                          ncol=2)


##Add x an y axis labels

y.grob <- textGrob("study duration", 
                   gp=gpar(  fontsize=11), rot=90)

x.grob <- textGrob("survey spatial coverage", 
                   gp=gpar(fontsize=11))

plots <- grid.arrange(arrangeGrob(plots.titles, left = y.grob, bottom = x.grob))

ggsave(file.path(plots.path, "num_good_fits.jpeg"), plot=plots, 
       width= 18, height=22, units="cm")






