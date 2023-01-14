
#############
# Calculate delta PD, lost PD, and gained PD in the plots 
#############

# Load libraries:
library(dplyr)
library(tidyr)
library(picante)

# Load phylogenetic tree:
tree<-read.tree("processed_data/phylo.tree.tre")

# Load vegetation data:
veg_data<-read.table("processed_data/veg_data_forestreplot.csv")

# Subset understory layer:
veg_data_und<- unique(subset(veg_data, layer == "H"))

# Get unique list of plots:
plots<-unique(veg_data_und$plotID)

# Prepare table for output:
out<-data.frame(plotID=plots, 
                herb_dif_sr=rep(NA, length(plots)), #delta SR
                herb_sr_lost=rep(NA, length(plots)), #lost SR
                herb_sr_gain=rep(NA, length(plots)), #gained SR
                herb_dif_pd=rep(NA, length(plots)),  #delta PD
                herb_pd_lost=rep(NA, length(plots)), #lost PD
                herb_pd_gain=rep(NA, length(plots))) #gained PD

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

  # Add row to matrix summing previous ones:
  veg_com[nrow(veg_com) + 1,] <- veg_com[1, ] + veg_com[2, ]
  veg_com[veg_com > 0] <- 1 #make sure it is only zeros and ones
  
  # Get SR and PD:
  res_pd<-pd(veg_com, tree)
  
  # Move results to the final table:
  out$herb_dif_sr[i]<-res_pd$SR[2]-res_pd$SR[1]
  out$herb_sr_gain[i]<-res_pd$SR[3]-res_pd$SR[1]
  out$herb_sr_lost[i]<-res_pd$SR[3]-res_pd$SR[2]
  
  out$herb_dif_pd[i]<-res_pd$PD[2]-res_pd$PD[1]
  out$herb_pd_gain[i]<-res_pd$PD[3]-res_pd$PD[1]
  out$herb_pd_lost[i]<-res_pd$PD[3]-res_pd$PD[2]
}

# Integrate with plot_data:
write.table(out, "processed_data/plot_data_PD_metrics.csv")

# Clean-up environment:
rm(list = ls())