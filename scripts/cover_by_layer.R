
#############
# Calculate covers of vegetation layers
#############

# Load packages
library(stringi)

# Load vegetation data:
veg_data<-read.table("data/veg_data_forestreplot.csv")

# Subset layer of interest (T=tree, S=shrub; H=herb):
veg_data_tree<- unique(subset(veg_data, layer == "T"))

# Get unique list of plots:
plots<-sort(unique(veg_data_tree$sample))

# Loop to get % cover of the layer:
out<-data.frame(sample=plots, tree_cover=rep(NA, length(plots))) #empty table to save results
for (i in 1:length(plots)){  
  sub_veg_data<-subset(veg_data_tree, sample==plots[i])
  out$tree_cover[i]<-1-exp(sum(log(1-(sub_veg_data$abundance/100)))) #apply Jennings-Fischer's formula
}

# Create new column with plotID:
out$plotID<-stri_replace_all_regex(out$sample,
                                   pattern=c('_B', '_R1', '_R2', '_R3', '_R4', '_R5'), #replace endings with sample time 
                                   replacement=c('', '', '', '', '', ''), #to get only ploIDs
                                   vectorize=FALSE)

# Calculate difference between baseline and resurvey:
pairs<-sort(unique(veg_data$plotID))

# Create empty data frame to save results:
lcv<-data.frame(plotID=pairs, 
                baseline_tree_cover=rep(NA, length(pairs)), 
                resurvey_tree_cover=rep(NA, length(pairs)), 
                dif_tree_cover=rep(NA, length(pairs)))

# Loop for each pair of surveys:
for (i in 1:length(pairs)){
  veg_com <- subset(out, plotID==pairs[i]) #subset the same plot at different times
  veg_com <- rbind(head(veg_com,1), tail(veg_com,1)) #merge first and last row
  lcv$dif_tree_cover[i]<-veg_com$tree_cover[2]-veg_com$tree_cover[1]
  lcv$baseline_tree_cover[i]<-veg_com$tree_cover[1]
  lcv$resurvey_tree_cover[i]<-veg_com$tree_cover[2]
}

# Save result:
write.table(lcv, "data/tree_covers.csv")

# Clean-up environment:
rm(list = ls())
