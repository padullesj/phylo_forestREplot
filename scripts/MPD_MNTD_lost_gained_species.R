
#############
# Calculate MPD and MNTD for lost and gained species
#############

# Load libraries:
library(dplyr)
library(tidyr)
library(picante)
library(PhyloMeasures)

# Load tree
tree<-read.tree("data/phylo.tree.tre")

# Load vegetation data:
veg_data<-read.table("data/veg_data_forestreplot.csv")

# Subset understory layer:
veg_data_und<- unique(subset(veg_data, layer == "H"))

# Get unique list of plots:
plots<-unique(veg_data_und$plotID)

# Prepare table for output:
out<-data.frame(plotID=plots, 
                herb_lost_mpd_obs=rep(NA, length(plots)), #Lost species MPD (observed)
                herb_lost_mntd_obs=rep(NA, length(plots)), #Lost species MNTD (observed)
                herb_lost_mpd_ses=rep(NA, length(plots)), #Lost species MPD (SES)
                herb_lost_mntd_ses=rep(NA, length(plots)), #Lost species MNTD (SES)
                herb_gain_mpd_obs=rep(NA, length(plots)), #Gained species MPD (observed)
                herb_gain_mntd_obs=rep(NA, length(plots)), #Gained species MNTD (observed)
                herb_gain_mpd_ses=rep(NA, length(plots)), #Gained species MPD (SES)
                herb_gain_mntd_ses=rep(NA, length(plots))) #Gained species MNTD (SES)

# Loop for each plot:
for (i in 1:length(plots)){ 
  print(i/length(plots)) #print progress
  
  # Subset surveys:
  sub_veg<-sort(unique((subset(veg_data_und, plotID==plots[i]))$sample)) #filter plot
  sub_veg<-veg_data_und[veg_data_und$sample %in% sub_veg, ] #subset plots from veg_data
  
  # Prepare table to get diversity values:
  veg_com <- unique(subset(sub_veg, select = c(sample, scientificName))) #subset columns
  veg_com$value<-1 #add  column with 1 as value
  veg_com<-reshape2::dcast(veg_com, sample ~ scientificName, value.var="value") #get long-format table
  rownames(veg_com)<-veg_com$sample #assign plot names to rownames
  veg_com$sample<-NULL #remove column with plot names
  veg_com <- veg_com  %>% mutate_all(funs(replace_na(.,0))) #replace NA with 0
  veg_com <- veg_com[, colSums(veg_com != 0) > 0] #remove species not present in the plot
  colnames(veg_com)<-gsub(" ", "_", colnames(veg_com)) #replace " " with "_" to match species names in the tree
  
  # Add rows with lost and gained species:
  veg_com[nrow(veg_com) + 1,] <- veg_com[1, ] - veg_com[2, ]
  veg_com[veg_com < 0] <- 0 #make sure it is only zeros and ones
  veg_com[nrow(veg_com) + 1,] <- veg_com[2, ] - veg_com[1, ]
  veg_com[veg_com < 0] <- 0 #make sure it is only zeros and ones
  veg_com<-veg_com[-c(1, 2), ] #remove original rows
  veg_com<-veg_com %>% select_if(colSums(.) != 0) #remove species that remained
  
  # Subset tree:
  treef<-keep.tip(tree, colnames(veg_com))
  
  # Get observed metrics:
  mpd<-mpd.query(treef, veg_com, standardize = F) #calculates MPD
  mntd<-mntd.query(treef, veg_com, standardize = F) #calculates MNTD
  
  # Save results:
  out$herb_lost_mpd_obs[i]<-mpd[1]
  out$herb_lost_mntd_obs[i]<-mntd[1]
  out$herb_gain_mpd_obs[i]<-mpd[2]
  out$herb_gain_mntd_obs[i]<-mntd[2]
  
  # Produce random values and get SES:
  rand_lost_mpd<-list()
  rand_lost_mntd<-list()
  rand_gain_mpd<-list()
  rand_gain_mntd<-list()
  
  for (j in 1:999){
    #print(j)
    set.seed(123+j)
    treef2<-tree #create copy of the tree
    treef2$tip.label<-sample(treef2$tip.label) #shuffle tips
    treef2<-keep.tip(treef2, colnames(veg_com)) #keep tips
    
    # Get metrics:
    mpd<-mpd.query(treef2, veg_com, standardize = F) #calculates MPD
    mntd<-mntd.query(treef2, veg_com, standardize = F) #calculates MNTD
    
    # Store results
    rand_lost_mpd[[j]]<-mpd[1]
    rand_lost_mntd[[j]]<-mntd[1]
    rand_gain_mpd[[j]]<-mpd[2]
    rand_gain_mntd[[j]]<-mntd[2]
  }
  
  # Calculate and save SES values:
  out$herb_lost_mpd_ses[i]<-(out$herb_lost_mpd_obs[i]-mean(unlist(rand_lost_mpd)))/sd(unlist(rand_lost_mpd))
  out$herb_lost_mntd_ses[i]<-(out$herb_lost_mntd_obs[i]-mean(unlist(rand_lost_mntd)))/sd(unlist(rand_lost_mntd))
  out$herb_gain_mpd_ses[i]<-(out$herb_gain_mpd_obs[i]-mean(unlist(rand_gain_mpd)))/sd(unlist(rand_gain_mpd))
  out$herb_gain_mntd_ses[i]<-(out$herb_gain_mntd_obs[i]-mean(unlist(rand_gain_mntd)))/sd(unlist(rand_gain_mntd))
}

# Save result:
write.table(out, "data/plot_data_MPD_MNTD_metrics.csv")

# Clean-up environment:
rm(list = ls())