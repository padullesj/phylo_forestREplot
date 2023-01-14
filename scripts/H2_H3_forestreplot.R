
####
# Core analyses for H2 and H3
####

#####################
#
# Table of contents
#
# A: Load plot and vegetation data
#
# B: H2
# B.1: Calculate U values
# B.2: Load phylogenetic tree and U values
# B.3: Calculate phylogenetic signal
# B.4: Get nodes with higher or lower mean U values
# B.5: Plot phylogenetic tree
#
# C: H3
# C.1: Load U values, trait data, and phylogenetic tree
# C.2: Integrate data
# C.3: Run PGLS
#
# D: Create final figure
#
#####################

# Load packages
library(plyr)
library(phytools)
library(ggplot2)
library(ggtree)
library(dplyr)
library(plyr)
library(scico)
library(grid)
library(cowplot)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(caper)

##
# A: Load plot and vegetation data
##

# Load plot data:
plot_data<-read.table("processed_data/plot_data_forestreplot.csv")

#Subset plots with resurvey time >= 20 years for H1, H2 & H3:
plot_data2<-subset(plot_data, time_resurvey>=20)

#Subset plots re-surveyed after 2000 for H1, H2 & H3:
plot_data2<-subset(plot_data2, year_final>2000)
plot_data2<-plot_data2[order(row.names(plot_data2)), ]

# Load vegetation data:
veg_data<-read.table("processed_data/veg_data_forestreplot.csv")

# Retain plots from the new "plot_data":
veg_data<-veg_data[veg_data$plotID %in% unique(plot_data2$plotID), ]

# Subset the understory layer:
veg_data<- unique(subset(veg_data, layer == "H"))

# Simplify vegetation data:
veg_data <- unique(subset(veg_data, select = c(scientificName, sample, plotID))) #simplify table

# Add "study" to vegetation data:
veg_data$study<-substr(veg_data$plotID, 1, 6) #add field

##
# B: H2
##

##
# B.1: Calculate U values
##

# Get list of species:
spp <-sort(unique(veg_data$scientificName)) #subset species occurring in at least 10 plots and 3 sites

# Get unique list of plots:
plots<-unique(veg_data$plotID)

# Prepare output table:
out<-data.frame(species=spp, uvalue=rep(NA, length(spp)))

# Loop to calculate U values per species
for (i in 1:length(spp)){ 
  print(i/length(spp)) #print progress
  sub_veg<-unique(subset(veg_data, scientificName==spp[i])) #subset species
  st_n<-length(unique(sub_veg$study)) #get number of study sites where species occurs
  
  if(st_n<3) { #if species occurs in less than 3 sites, assign NA
    out$uvalue[i]<-NA
    
  } else { #otherwise
    
    # Get plot names:
    plots<-unique(sub_veg$plotID) #get unique list of plots where the species occurs

    # Determine if a species was gained (1), lost (-1), or remained present (0) over time
    out2<-NULL #create empty object

    for (j in 1:length(plots)){
      sub_veg2<-subset(sub_veg, plotID==plots[j])
      
      if(nrow(sub_veg2) == 2) {
        a<-0
      } else {
        if(endsWith(sub_veg2$sample, 'B')) {
          a<--1
        } else {
          a<-1
        }
      }
      out2<-c(out2, a)
    }
    
    # Merge vectors into table:
    res<-data.frame(plotID=plots, change=out2)
    res$study<-substr(res$plotID, 1, 6)
    
    # Get number of plots per site:
    b<-table(res$study)
    
    # Determine if 3 of the sites contain at least 10 plots:
    if(length(b)>2) {
      c<-max(combn(b, 3, FUN = sum))
    } else {
      c<-2 #random number under the threshold
    }
    
    # Calculate values per species:
    if(c<10) { #if less than 10 plots in at least 3 sites, assign NA
      out$uvalue[i]<-NA
      
    } else {
      val<-list()
      for (k in 1:1000){ #loop 1000 times subsetting 10 plots stratified within sites
        #print(k)
        set.seed(1234+k) #set seed for randomization
        
        # Repeat process until we have 3 sites with at least 10 plots:
        success <- FALSE
        while (!success) {
          # do something
          sites<-sort(sample(b, 3, replace=FALSE))
          nplot<-sum(sites)
          # check for success
          success <- nplot >= 10
        }
        
        # Subset table with only specific sites:
        res3 <- res[res$study %in% names(sites), ]
        
        # Site 1:
        if(sites[1]>3) { #subset X number of plots
          sel_plots1<-sample(res3$plotID[res3$study %in% names(sites)[1]], 3, replace=F)
        } else {
          sel_plots1<-res3$plotID[res3$study %in% names(sites)[1]]
        }
        
        # Define number of plots to sample in Site 2:
        if(length(sel_plots1)<3) { #subset X number of plots
          stt<-4 } else { stt<-3 }
        
        # Site 2:
        if(sites[2]>3) { #subset X number of plots
          sel_plots2<-sample(res3$plotID[res3$study %in% names(sites)[2]], stt, replace=F)
        } else {
          sel_plots2<-res3$plotID[res3$study %in% names(sites)[2]]
        }
        
        # Site 3:
        stt<-10-(length(sel_plots1)+length(sel_plots2))
        sel_plots3<-sample(res3$plotID[res3$study %in% names(sites)[3]], stt, replace=F)
        
        # Put plots together:
        res4<-res[res$plotID %in% c(sel_plots1, sel_plots2, sel_plots3), ]
        
        # Get weighted mean:
        val[[k]]<-mean(res4$change) #calculate the weighted mean
      }
      out$uvalue[i]<-mean(unlist(val)) #get the mean across all randomizations
    }
  }
}

# Save result:
write.table(out, "processed_data/species_U_values.csv")


##
# B.2: Load phylogenetic tree and U values
##

# Load U values:
out<-read.table("processed_data/species_U_values.csv")

# Assign families to species and add to table with the U values:
fam<-read.table("processed_data/veg_data_forestreplot.csv")
fam<-unique(subset(fam, select = c(scientificName, family))) #simplify table
names(fam)[1]<-"species"
out2<-merge(out, fam, by="species", all.x=T)

# Reshape table:
out$species<-gsub(" ", "_", out$species) #change names of species
out<-out[!is.na(out$uvalue),] #remove NA
rownames(out)<-out$species #assign row names
out$species<-NULL #delete unnecessary field

# Load tree:
tree<-read.tree("processed_data/phylo.tree.tre")
treef<-keep.tip(tree, rownames(out)) #subset species present in the U values list


##
# B.3: Calculate phylogenetic signal
##

out2<-as.numeric(out$uvalue)
names(out2)<-rownames(out)
phylosig(treef, out2, method="lambda", test=TRUE, nsim=999)


##
# B.4: Get nodes with higher or lower mean U values
##

# Run the node.mean function (load from "node.mean_function_forestreplot)
# to detect nodes with higher or lower tendencies than under random expectation:
nodef<-node.mean(treef, out, 999) #this can take a while
#write.table(nodef, "processed_data/node_forestreplot.csv") #save the result
#nodef<-read.table("processed_data/node_forestreplot.csv") #load the result

# Clean-up result to highlight nodes in the tree:
significant<-nodef #create a copy of the main result
significant$P_value[significant$P_value>0.05]<-NA #replace non-significant with NA
significant$P_value[significant$SR_mis<3]<-NA #assign NA to nodes with 2 or less species
significant$P_value[1]<-NA #set the first node (the root node) to NA
significant$P_value  <- with(significant, ifelse(Obs>significant$Mean_Exp & P_value<0.05, "pos.05", P_value)) #identify significantly higher at alpha < .05
significant$P_value  <- with(significant, ifelse(Obs<significant$Mean_Exp & P_value<0.05, "neg.05", P_value)) #identify significantly lower at alpha < .05
significant<-c(rep(NA, length(treef$tip.label)), significant$P_value) #merge "tip nodes" with "inner" tree nodes
significant<-as.factor(significant) #convert values into factors
significant<-factor(significant, levels = c("neg.05", "pos.05" )) #change order of factors


###
# B.5: Plot phylogenetic tree
###

# Get list of taxa by family:
fam<-read.table("processed_data/veg_data_forestreplot.csv")
fam<-unique(subset(fam, select = c(scientificName, family))) #simplify table
fam$scientificName<-gsub(" ", "_", fam$scientificName) #adapt species nomenclature
fam<-fam[which(fam$scientificName %in% rownames(out)),] #subset species included in wl

# Get table with ranked families based on their number of species:
famf<-as.data.frame(table(fam$family)) #create dataframe
famf<-famf[order(-famf$Freq),] #order it
row.names(famf) <- NULL #get rid of rownames

# Create vector with unique names of species:
list.r <- unique(fam$family)

# Get nodes for families:
list.nod<-NULL #create empty object to store the results
for (i in 1:length(list.r)) {
  #print(i)
  tips3 <- as.vector(subset(fam, family == list.r[i])$scientificName) #vector with all species belonging to the given family
  if (length(tips3)>1) {
    num <- phytools::findMRCA(treef, tips3) #get node that contain those species
    num <-data.frame(list.r[i], num) #assign family name to node
    list.nod<-rbind(list.nod, num) #save and merge with other families
  } 
}

# Clean-up the result and merge with previous table:
names(list.nod)[1]<-paste("Var1")
famf<-join(famf,list.nod)

#To colour branches based on mean values:
svl <- as.matrix(out)[,1]
fit <- phytools::fastAnc(treef, svl, vars=TRUE, CI=TRUE)
td <- data.frame(node = nodeid(treef, names(svl)),
                 trait = svl)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)
treef <- full_join(treef, d, by = 'node')

# Get vector with names of families containing more species:
toplot<-as.character(head(famf$Var1, n=50)) #select the top 50 families with 5 for more species

# Plot tree:
p <- 
  ggtree(treef, layout="circular", size=0.5)+ # build circular tree
  geom_tree(aes(color=trait), continuous = 'colour', show.legend = F) +
  scale_color_scico(palette = "vik", direction=-1, na.value = "gray48", limits = c(-1, 1)) +
  
  ggnewscale::new_scale("color") +
  
  geom_point(aes(color=as.factor(significant)), size=2, alpha=1, show.legend = F) + # highlight nodes
  scale_colour_manual(values=c("#8A6000", "#006FA4"), labels=c("Loss", "Gain"), na.translate=FALSE)+ # set aesthetics for highlighted nodes
  
  geom_cladelabel(node=subset(famf, Var1==toplot[1])$num, label=toplot[1], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[2])$num, label=toplot[2], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[3])$num, label=toplot[3], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 15) +
  geom_cladelabel(node=subset(famf, Var1==toplot[4])$num, label=toplot[4], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[5])$num, label=toplot[5], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[6])$num, label=toplot[6], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 339) +
  geom_cladelabel(node=subset(famf, Var1==toplot[7])$num, label=toplot[7], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 278) +
  geom_cladelabel(node=subset(famf, Var1==toplot[8])$num, label=toplot[8], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[9])$num, label=toplot[9], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 308) +
  geom_cladelabel(node=subset(famf, Var1==toplot[10])$num, label=toplot[10], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 40) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[11])$num, label=toplot[11], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 71) +
  geom_cladelabel(node=subset(famf, Var1==toplot[12])$num, label=toplot[12], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[13])$num, label=toplot[13], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[14])$num, label=toplot[14], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 53) +
  geom_cladelabel(node=subset(famf, Var1==toplot[15])$num, label=toplot[15], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 78) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[16])$num, label=toplot[16], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[17])$num, label=toplot[17], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 87) +
  geom_cladelabel(node=subset(famf, Var1==toplot[18])$num, label=toplot[18], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 31) +
  geom_cladelabel(node=subset(famf, Var1==toplot[19])$num, label=toplot[19], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 328) +
  geom_cladelabel(node=subset(famf, Var1==toplot[20])$num, label=toplot[20], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[21])$num, label=toplot[21], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[22])$num, label=toplot[22], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 355) +
  geom_cladelabel(node=subset(famf, Var1==toplot[23])$num, label=toplot[23], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 318) +
  geom_cladelabel(node=subset(famf, Var1==toplot[24])$num, label=toplot[24], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 350) +
  geom_cladelabel(node=subset(famf, Var1==toplot[25])$num, label=toplot[25], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 295) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[26])$num, label=toplot[26], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[27])$num, label=toplot[27], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[28])$num, label=toplot[28], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 288) +
  geom_cladelabel(node=subset(famf, Var1==toplot[29])$num, label=toplot[29], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 322) +
  geom_cladelabel(node=subset(famf, Var1==toplot[30])$num, label=toplot[30], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[31])$num, label=toplot[31], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[32])$num, label=toplot[32], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 60) +
  geom_cladelabel(node=subset(famf, Var1==toplot[33])$num, label=toplot[33], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 47) +
  geom_cladelabel(node=subset(famf, Var1==toplot[34])$num, label=toplot[34], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 49) +
  geom_cladelabel(node=subset(famf, Var1==toplot[35])$num, label=toplot[35], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 291) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[36])$num, label=toplot[36], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[37])$num, label=toplot[37], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[38])$num, label=toplot[38], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 298) +
  geom_cladelabel(node=subset(famf, Var1==toplot[39])$num, label=toplot[39], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 2) +
  geom_cladelabel(node=subset(famf, Var1==toplot[40])$num, label=toplot[40], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 63) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[41])$num, label=toplot[41], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[42])$num, label=toplot[42], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 314) +
  geom_cladelabel(node=subset(famf, Var1==toplot[43])$num, label=toplot[43], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 58) +
  geom_cladelabel(node=subset(famf, Var1==toplot[44])$num, label=toplot[44], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 44) +
  geom_cladelabel(node=subset(famf, Var1==toplot[45])$num, label=toplot[45], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 284) +
  
  geom_cladelabel(node=subset(famf, Var1==toplot[46])$num, label=toplot[46], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[47])$num, label=toplot[47], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 301) +
  geom_cladelabel(node=subset(famf, Var1==toplot[48])$num, label=toplot[48], offset=5, fontsize=2.8, barsize = 0.2, angle = "auto") +
  geom_cladelabel(node=subset(famf, Var1==toplot[49])$num, label=toplot[49], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 82) +
  geom_cladelabel(node=subset(famf, Var1==toplot[50])$num, label=toplot[50], offset=5, fontsize=2.8, barsize = 0.2, hjust=1, angle = 358) +
  
  theme(plot.title = element_text(size = 23, face = "bold", hjust=0.5),
        legend.title=element_blank(), 
        legend.text=element_text(size=16),
        legend.key.size = unit(1, "cm"),
        legend.position="none")

# Add heatmap with tendency values:
p1 <- gheatmap(p, out, offset=0.3, width=.03, colnames = F, color = NULL) +
  scale_fill_scico(palette = "vik", direction=-1, na.value = "gray48",
                   limits = c(-1, 1)) +
  theme(legend.position="bottom",
        legend.text=element_text(size=15),
        legend.title=element_text(colour="white"))


##
# C: H3
##

##
# C.1: Load U values, trait data, and phylogenetic tree
##

# Load U values:
out<-read.table("processed_data/species_U_values.csv")
out$species<-gsub(" ", "_", out$species) #change names of species
out<-out[!is.na(out$uvalue),] #remove NA

# Load trait data:
spp_traits<-read.table("processed_data/traits_forestreplot.csv")
spp_traits$species<-gsub(" ", "_", spp_traits$species) #create column with species names
rownames(spp_traits)<-spp_traits$species #assign species as row names

# Load tree:
tree<-read.tree("processed_data/phylo.tree.tre")


##
# C.2: Integrate data
##

# Merge traits and U values:
spp_traits<-merge(out, spp_traits, by="species", all.x=T)
spp_traits<- spp_traits[complete.cases(spp_traits), ] #get complete cases

# Subset tree:
treef<-keep.tip(tree, spp_traits$species) #subset species present in the plots
treef$node.label<-NULL


##
# C.3: Run PGLS
##

# Merge spp_traits with phylogenetic tree:
dat <- comparative.data(treef, spp_traits, species, vcv=TRUE, vcv.dim=3) #assemble data

# Run model:
m1 <- pgls(uvalue ~ 
           scale(log(sla)) + scale(log(leaf_area)) + scale(log(seed_mass)) + scale(log(plant_height)),
           dat, lambda='ML')

# Extract coefficients:
rat<-as.data.frame(summary(m1)$coefficients)
rat<-rat[-c(1), ] #remove Intercept
rat$var<-c("SLA", "LA", "SM", "H") #assign trait names
rat$var <- factor(rat$var, levels = c("SLA", "LA", "SM", "H")) #changer order of factors

# Set confidence interval:
rat$low<-rat$Estimate - (rat$`Std. Error`* 1.959)
rat$up<-rat$Estimate + (rat$`Std. Error`* 1.959)

# Define groups of variables:
rat$col <- ifelse((rat$low < 0) & (rat$up < 0),"green", ifelse((rat$low > 0) & (rat$up > 0),"orange", ifelse("black")))
rat$col[is.na(rat$col)] <- "black"

# Plot effects:
p_rat<-ggplot(rat, aes(x=var, y=Estimate, ymin=low, ymax=up, shape=col)) +  
  geom_errorbar(size=0.5, width = 0.2) + 
  geom_point(size=5) +
  scale_shape_manual(values=c(1, 16, 16))+
  scale_y_continuous(breaks = c(0, 0.1)) +
  #scale_color_manual(values=c("black", "#4BB446", "#AF46B4")) +
  geom_hline(yintercept=0, lty=2, size=0.25) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Estimates \n(± 95% confident interval)") + 
  #ylim(-0.7, 0.7)+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position="none",
        axis.title.x = element_text(size = 8),
        plot.title = element_text(color="black", hjust = 0.5, size=8, face="bold"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"))


##
# D: Create final figure
##

# Create layout to print:
h <- ggdraw(p1)
ger <- ggplotGrob(p_rat) #boxplot
p2<- h + draw_grob(ger, 0.005, 0.01, 0.24, 0.22)

# Save result:
png("results/Fig2.png",
    res=300,height=8,width=8,units="in"); 
p2
grid.text("Gain", x = unit(0.7, "npc"), y = unit(0.08, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#006FA4"))
grid.text("Loss", x = unit(0.35, "npc"), y = unit(0.08, "npc"), gp=gpar(fontsize=17, fontface="bold", col="#8A6000"))

grid.text("a)", x = unit(0.01, "npc"), y = unit(0.98, "npc"), hjust=0, gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("b)", x = unit(0.01, "npc"), y = unit(0.25, "npc"), hjust=0, gp=gpar(fontsize=16, fontface="bold", col="black"))

grid.text("Gain", x = unit(0.18, "npc"), y = unit(0.25, "npc"), gp=gpar(fontsize=11, fontface="bold"))
grid.text("Loss", x = unit(0.09, "npc"), y = unit(0.25, "npc"), gp=gpar(fontsize=11, fontface="bold"))

grid.lines(x = c(0.11, 0.2),  y = 0.235, gp = gpar(col = "black"), arrow = arrow(angle = 30, length = unit(0.07, "inches"),
                                                                                 ends = "last", type = "closed"))
grid.lines(x = c(0.11, 0.07),  y = 0.235, gp = gpar(col = "black"), arrow = arrow(angle = 30, length = unit(0.07, "inches"),
                                                                                  ends = "last", type = "closed"))

dev.off()

# Clean-up environment:
rm(list = ls())