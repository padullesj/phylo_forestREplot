
#############
# Calculate phylogenetic diversity (PD) relatedness (PR) metrics
#############

# Load libraries:
library(dplyr)
library(tidyr)
library(PhyloMeasures)
library(ape)

# Load phylogenetic tree:
tree<-read.tree("data/phylo.tree.tre")

# Load plot data:
plot_time<-read.table("data/plot_data_forestreplot.csv")

# Load vegetation data:
veg_data<-read.table("data/veg_data_forestreplot.csv")

# Subset understory layer:
veg_data_und<- unique(subset(veg_data, layer == "H"))

# Get unique list of plots:
plots<-unique(veg_data_und$plotID)

# Remove plots with no herb layer in both surveys:
rem<-NULL
for (i in 1:length(plots)){
  sub_veg<-sort(unique((subset(veg_data_und, plotID==plots[i]))$sample)) #filter plot
  if(length(sub_veg)==1) {
    rem<-append(rem, plots[i])
  }
}
veg_data_und<-veg_data_und[!veg_data_und$plotID %in% rem, ]

# Get new list of plots:
plots<-unique(veg_data_und$plotID)

# Prepare table for output:
out<-data.frame(plotID=plots,
                herb_delta_sr=rep(NA, length(plots)), 
                herb_lost_sr=rep(NA, length(plots)),
                herb_gain_sr=rep(NA, length(plots)),
                herb_delta_pd_obs=rep(NA, length(plots)), herb_delta_pd_ses=rep(NA, length(plots)),
                herb_lost_pd_obs=rep(NA, length(plots)), herb_lost_pd_ses=rep(NA, length(plots)),
                herb_gain_pd_obs=rep(NA, length(plots)), herb_gain_pd_ses=rep(NA, length(plots)),
                herb_lost_mpd_obs=rep(NA, length(plots)), herb_lost_mpd_ses=rep(NA, length(plots)),
                herb_lost_mntd_obs=rep(NA, length(plots)), herb_lost_mntd_ses=rep(NA, length(plots)),
                herb_gain_mpd_obs=rep(NA, length(plots)), herb_gain_mpd_ses=rep(NA, length(plots)),
                herb_gain_mntd_obs=rep(NA, length(plots)), herb_gain_mntd_ses=rep(NA, length(plots)), 
                herb_lost_dpw_obs=rep(NA, length(plots)), herb_lost_dpw_ses=rep(NA, length(plots)),
                herb_lost_dnn_obs=rep(NA, length(plots)), herb_lost_dnn_ses=rep(NA, length(plots)),
                herb_gain_dpw_obs=rep(NA, length(plots)), herb_gain_dpw_ses=rep(NA, length(plots)),
                herb_gain_dnn_obs=rep(NA, length(plots)), herb_gain_dnn_ses=rep(NA, length(plots))
)

# Loop for each plot:
for (i in 1:length(plots)){ 
  print(i/length(plots)) #print progress
  
  # Get the number of years between surveys:
  yr<-subset(plot_time, plotID==plots[i])
  yr<-yr$year_final-yr$year_baseline_survey
  
  # Keep plots from the baseline survey and the resurvey:
  sub_veg<-sort(unique((subset(veg_data_und, plotID==plots[i]))$sample)) #filter plot
  sub_veg<-c(head(sub_veg, 1), tail(sub_veg, 1)) #select the baseline and final plots
  sub_veg<-veg_data_und[veg_data_und$sample %in% sub_veg, ] #subset plots from veg_data

  # Prepare table to get diversity values:
  veg_com <- subset(sub_veg, angiosperm="angiosperm") #subset angiosperms
  veg_com <- unique(subset(veg_com, select = c(sample, species_id))) #subset columns
  veg_com$value<-1 #add  column with 1 as value
  veg_com<-reshape2::dcast(veg_com, sample ~ species_id, value.var="value") #get long-format table
  rownames(veg_com)<-veg_com$sample #assign plot names to row names
  veg_com$sample<-NULL #remove column with plot names
  veg_com <- veg_com  %>% mutate_all(funs(replace_na(.,0))) #replace NA with 0
  veg_com <- veg_com[, colSums(veg_com != 0) > 0] #remove species not present in the plot
  colnames(veg_com)<-gsub(" ", "_", colnames(veg_com)) #replace " " with "_" to match species names in the tree
  
  # Create copy of community data:
  veg_com2<-veg_com
  
  # Create table with lost and gained species:
  veg_com2[nrow(veg_com2) + 1,] <- veg_com2[1, ] - veg_com2[2, ]
  veg_com2[veg_com2 < 0] <- 0 #make sure it is only zeros and ones
  veg_com2[nrow(veg_com2) + 1,] <- veg_com2[2, ] - veg_com2[1, ]
  veg_com2[veg_com2 < 0] <- 0 #make sure it is only zeros and ones
  veg_com2<-veg_com2[-c(1, 2), ] #remove original rows
  veg_com2<-veg_com2 %>% select_if(colSums(.) != 0) #remove persisting species
  
  # Create table with lost vs. persisting species:
  veg_com3<-veg_com
  veg_com3[nrow(veg_com3) + 1,] <- veg_com3[1, ] - veg_com3[2, ] #row with lost species
  veg_com3[veg_com3 < 0] <- 0 #make sure it is only zeros and ones
  veg_com3[nrow(veg_com3) + 1,] <- veg_com3[1, ] + veg_com3[2, ] #row for persisting species
  veg_com3[c(4),][veg_com3[c(4),]==1] <- 0 #make sure it is only zeros and ones
  veg_com3[c(4),][veg_com3[c(4),]==2] <- 1 #make sure it is only zeros and ones
  veg_com3<-veg_com3[-c(1, 2), ] #remove original rows
  veg_com3<-veg_com3 %>% select_if(colSums(.) != 0) #remove persisting species that remained
  
  # Create table with gained vs. persisting species:
  veg_com4<-veg_com
  veg_com4[nrow(veg_com4) + 1,] <- veg_com4[2, ] - veg_com4[1, ] #row with gained species
  veg_com4[veg_com4 < 0] <- 0 #make sure it is only zeros and ones
  veg_com4[nrow(veg_com4) + 1,] <- veg_com4[1, ] + veg_com4[2, ] #row for persisting species
  veg_com4[c(4),][veg_com4[c(4),]==1] <- 0 #make sure it is only zeros and ones
  veg_com4[c(4),][veg_com4[c(4),]==2] <- 1 #make sure it is only zeros and ones
  veg_com4<-veg_com4[-c(1, 2), ] #remove original rows
  veg_com4<-veg_com4 %>% select_if(colSums(.) != 0) #remove species that remained
  
  # Subset tree:
  tree1<-keep.tip(tree, colnames(veg_com))
  tree2<-keep.tip(tree, colnames(veg_com2))
  tree3<-keep.tip(tree, colnames(veg_com3))
  tree4<-keep.tip(tree, colnames(veg_com4))
  
  # Get PD for all species:
  pd1<-pd.query(tree1, veg_com, standardize = FALSE)
  
  # Get observed metrics for SR, PD, and PR:
  if(min(rowSums(veg_com2))==0) {sr<-c(0, 0)
  } else {sr<-rowSums(veg_com2)} #SR
  
  if(min(rowSums(veg_com2))==0) {pd2<-c(0, 0)
  } else {pd2<-pd.query(tree2, veg_com2, standardize = FALSE)} #PD
  
  if(min(rowSums(veg_com2))==0) {mpd<-c(NA, NA)
  } else {mpd<-mpd.query(tree2, veg_com2, standardize = FALSE)} #MPD
  
  if(min(rowSums(veg_com2))==0) {mntd<-c(NA, NA)
  } else {mntd<-mntd.query(tree2, veg_com2, standardize = FALSE)} #MNTD
  
  # Get observed metrics for Dpw:
  if(min(rowSums(veg_com3))==0) {dpw1<-NA
  } else {dpw1<-cd.query(tree3, veg_com3, standardize = FALSE)[1,2]} #Dpw
  
  if(min(rowSums(veg_com4))==0) {dpw2<-NA
  } else {dpw2<-cd.query(tree4, veg_com4, standardize = FALSE)[1,2]} #Dpw
  
  # Get observed metrics for Dnn:
  if(min(rowSums(veg_com3))==0) {dnn1<-NA
  } else {dnn1<-cdnt.query(tree3, veg_com3)[1,2]} #Dnn
  
  if(min(rowSums(veg_com4))==0) {dnn2<-NA
  } else {dnn2<-cdnt.query(tree4, veg_com4)[1,2]} #Dnn
  
  # Save results:
  out$herb_delta_sr[i]<-(log(rowSums(veg_com)[2]/rowSums(veg_com)[1]))/yr
  out$herb_delta_pd_obs[i]<-(log(pd1[2]/pd1[1]))/yr
  out$herb_lost_sr[i]<-sr[1]
  out$herb_gain_sr[i]<-sr[2]
  out$herb_lost_pd_obs[i]<-pd2[1]
  out$herb_gain_pd_obs[i]<-pd2[2]
  out$herb_lost_mpd_obs[i]<-mpd[1]
  out$herb_gain_mpd_obs[i]<-mpd[2]
  out$herb_lost_mntd_obs[i]<-mntd[1]
  out$herb_gain_mntd_obs[i]<-mntd[2]
  out$herb_lost_dpw_obs[i]<-dpw1[1]
  out$herb_gain_dpw_obs[i]<-dpw2[1]
  out$herb_lost_dnn_obs[i]<-dnn1[1]
  out$herb_gain_dnn_obs[i]<-dnn2[1]
  
  # Produce random values and get SES:
  rand_delta_pd<-list()
  rand_lost_pd<-list()
  rand_gain_pd<-list()
  rand_lost_mpd<-list()
  rand_gain_mpd<-list()
  rand_lost_mntd<-list()
  rand_gain_mntd<-list()
  rand_lost_dpw<-list()
  rand_gain_dpw<-list()
  rand_lost_dnn<-list()
  rand_gain_dnn<-list()
  
  for (j in 1:999){ #repeat 999 times
    #print(j)
    set.seed(12+j, kind="Mersenne-Twister", normal.kind = "Inversion")
    treer<-tree #create copy of the tree
    treer$tip.label<-sample(treer$tip.label) #shuffle tips
    
    #Create trees:
    tree1<-keep.tip(treer, colnames(veg_com)) #keep tips
    tree2<-keep.tip(treer, colnames(veg_com2)) #keep tips
    tree3<-keep.tip(treer, colnames(veg_com3)) #keep tips
    tree4<-keep.tip(treer, colnames(veg_com4)) #keep tips
    
    # Get PD for all species:
    pd1<-pd.query(tree1, veg_com, standardize = FALSE)
    
    # Get metrics for PD and PR:
    if(min(rowSums(veg_com2))==0) {pd2<-c(0, 0)
    } else {pd2<-pd.query(tree2, veg_com2, standardize = FALSE)} #PD
    
    if(min(rowSums(veg_com2))==0) {mpd<-c(NA, NA)
    } else {mpd<-mpd.query(tree2, veg_com2, standardize = FALSE)} #MPD  
    
    if(min(rowSums(veg_com2))==0) {mntd<-c(NA, NA)
    } else {mntd<-mntd.query(tree2, veg_com2, standardize = FALSE)} #MNTD
    
    # Get observed metrics for Dpw:
    if(min(rowSums(veg_com3))==0) {dpw1<-NA
    } else {dpw1<-cd.query(tree3, veg_com3, standardize = FALSE)[1,2]} #Dpw
    
    if(min(rowSums(veg_com4))==0) {dpw2<-NA
    } else {dpw2<-cd.query(tree4, veg_com4, standardize = FALSE)[1,2]} #Dpw
    
    # Get observed metrics for Dnn:
    if(min(rowSums(veg_com3))==0) {dnn1<-NA
    } else {dnn1<-cdnt.query(tree3, veg_com3)[1,2]} #Dnn
    
    if(min(rowSums(veg_com4))==0) {dnn2<-NA
    } else {dnn2<-cdnt.query(tree4, veg_com4)[1,2]} #Dnn
    
    # Store results:
    rand_delta_pd[[j]]<-(log(pd1[2]/pd1[1]))/yr
    rand_lost_pd[[j]]<-pd2[1]
    rand_gain_pd[[j]]<-pd2[2]
    rand_lost_mpd[[j]]<-mpd[1]
    rand_gain_mpd[[j]]<-mpd[2]
    rand_lost_mntd[[j]]<-mntd[1]
    rand_gain_mntd[[j]]<-mntd[2]
    rand_lost_dpw[[j]]<-dpw1[1]
    rand_gain_dpw[[j]]<-dpw2[1]
    rand_lost_dnn[[j]]<-dnn1[1]
    rand_gain_dnn[[j]]<-dnn2[1]
    
  }
  
  # Calculate and save SES values:
  out$herb_delta_pd_ses[i]<-(out$herb_delta_pd_obs[i]-mean(unlist(rand_delta_pd)))/sd(unlist(rand_delta_pd))
  out$herb_lost_pd_ses[i]<-(out$herb_lost_pd_obs[i]-mean(unlist(rand_lost_pd)))/sd(unlist(rand_lost_pd))
  out$herb_gain_pd_ses[i]<-(out$herb_gain_pd_obs[i]-mean(unlist(rand_gain_pd)))/sd(unlist(rand_gain_pd))
  out$herb_lost_mpd_ses[i]<-(out$herb_lost_mpd_obs[i]-mean(unlist(rand_lost_mpd)))/sd(unlist(rand_lost_mpd))
  out$herb_gain_mpd_ses[i]<-(out$herb_gain_mpd_obs[i]-mean(unlist(rand_gain_mpd)))/sd(unlist(rand_gain_mpd))
  out$herb_lost_mntd_ses[i]<-(out$herb_lost_mntd_obs[i]-mean(unlist(rand_lost_mntd)))/sd(unlist(rand_lost_mntd))
  out$herb_gain_mntd_ses[i]<-(out$herb_gain_mntd_obs[i]-mean(unlist(rand_gain_mntd)))/sd(unlist(rand_gain_mntd))
  out$herb_lost_dpw_ses[i]<-(out$herb_lost_dpw_obs[i]-mean(unlist(rand_lost_dpw)))/sd(unlist(rand_lost_dpw))
  out$herb_gain_dpw_ses[i]<-(out$herb_gain_dpw_obs[i]-mean(unlist(rand_gain_dpw)))/sd(unlist(rand_gain_dpw))
  out$herb_lost_dnn_ses[i]<-(out$herb_lost_dnn_obs[i]-mean(unlist(rand_lost_dnn)))/sd(unlist(rand_lost_dnn))
  out$herb_gain_dnn_ses[i]<-(out$herb_gain_dnn_obs[i]-mean(unlist(rand_gain_dnn)))/sd(unlist(rand_gain_dnn))
  
}

# Save final table:
write.table(out, "data/plot_data_PD_PR_metrics.csv")

# Clean-up environment:
rm(list = ls())