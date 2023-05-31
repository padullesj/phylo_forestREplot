
#############
# Calculate Dpw and Dnn focusing on lost species
#############

# Load libraries:
library(dplyr)
library(tidyr)
library(picante)

# Load tree
tree<-read.tree("data/phylo.tree.tre")

# Calculate patristic distances between species in the tree
tree_dist<-as.data.frame(as.matrix(adephylo::distTips(tree, method = "patristic"))) 

# Load vegetation data:
veg_data<-read.table("data/veg_data_forestreplot.csv")

# Subset understory layer:
veg_data_und<- unique(subset(veg_data, layer == "H"))

# Get unique list of plots:
plots<-unique(veg_data_und$plotID)

# Prepare table for output:
out<-data.frame(plotID=plots, 
                herb_dnn_obs=rep(NA, length(plots)), #Dnn (observed)
                herb_dnn_ses=rep(NA, length(plots)), #Dnn (SES)
                herb_dpw_obs=rep(NA, length(plots)), #Dpw (observed)
                herb_dpw_ses=rep(NA, length(plots))) #Dpw (SES)

# Loop for each plot:
for (i in 1:length(plots)){ 
  print(i/length(plots)) #print progress
  
  # Keep plots from baseline and last year
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

  # Subset lost and gained species:
  veg_com2<-as.data.frame(t(veg_com))
  spp_lost<-rownames(subset(veg_com2, `3`==1))
  spp_gain<-rownames(subset(veg_com2, `4`==1))
  
  if(length(spp_lost)<1 | length(spp_gain)<1) { #assign NA if there were no lost and gained species
    out$herb_dnn_obs[i]<-NA
    out$herb_dnn_ses[i]<-NA
  } else { #otherwise
    
    # Subset columns with lost species:
    res<-subset(tree_dist, select = spp_lost)
    
    # Subset rows with gained species
    res<-subset(res, rownames(res) %in% spp_gain)
    
    # Get the observed average of all minimum (Dnn) and mean (Dpw) values across columns:
    out$herb_dnn_obs[i]<-mean(apply(res,2,min))
    out$herb_dpw_obs[i]<-mean(apply(res,2,mean))
    
    # Produce random values and get SES:
    rand_lost_mntd<-list()
    rand_lost_mpd<-list()
    
    for (j in 1:999){
      #print(j)
      set.seed(123+j)
      tree_dist2<-tree_dist #create copy of the tree
      rownames(tree_dist2)<-sample(rownames(tree_dist2)) #shuffle tips
      colnames(tree_dist2)<-rownames(tree_dist2)
      
      # Subset columns with lost species:
      res<-subset(tree_dist2, select = spp_lost)
      
      # Subset rows with gained species:
      res<-subset(res, rownames(res) %in% spp_gain)
      
      # Get the random average of all minimum (Dnn) and mean (Dpw) values across columns:
      rand_lost_mntd[[j]]<-mean(apply(res,2,min))
      rand_lost_mpd[[j]]<-mean(apply(res,2,mean))
    }
    out$herb_dnn_ses[i]<-(out$herb_dnn_obs[i]-mean(unlist(rand_lost_mntd)))/sd(unlist(rand_lost_mntd))
    out$herb_dpw_ses[i]<-(out$herb_dpw_obs[i]-mean(unlist(rand_lost_mpd)))/sd(unlist(rand_lost_mpd))
  }
}

# Save result:
write.table(out, "data/plot_data_Dpw_Dnn_metrics.csv")

# Clean-up environment:
rm(list = ls())
