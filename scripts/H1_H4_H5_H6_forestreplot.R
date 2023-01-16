
####
# Core analyses for H1, H4, H5 and H6
####

#####################
#
# Table of contents
#
# A: Load plot and vegetation data
#
# B: H1
#
# C: H4 and H5
#
# D: H6
#
#####################

# Load packages:
library(nlme)
library(ggplot2)
library(ggpubr)
library(grid)
library(scico)
library(lsr)
library(performance)
library(MuMIn)


##
# A: Load and clean-up data
##

# Load vegetation data:
veg_data<-read.table("data/veg_data_forestreplot.csv")

# Load plot data:
plot_data<-read.table("data/plot_data_forestreplot.csv")

# Calculate temporal rates:
plot_data$herb_dif_pd_rr<-plot_data$herb_dif_pd/plot_data$time_resurvey
plot_data$herb_lost_pd_rr<-plot_data$herb_lost_pd/plot_data$time_resurvey
plot_data$herb_gain_pd_rr<-plot_data$herb_gain_pd/plot_data$time_resurvey

#Subset plots with resurvey time >= 20 years for H1, H2 & H3:
plot_data2<-subset(plot_data, time_resurvey>=20)

#Subset plots re-surveyed after 2000 for H1, H2 & H3:
plot_data2<-subset(plot_data2, year_final>2000)
plot_data2<-plot_data2[order(row.names(plot_data2)), ]

#Subset rows with non NA:
plot_data2<-subset(plot_data2, !is.na(plot_data2$herb_dif_pd_rr))

##
# B: H1
##

# Calculate Cohen's D:
cohensD(plot_data2$herb_dif_pd_rr, mu=0) #test H1

#Plot histograms for Figure 1:
a<-ggplot(plot_data2, aes(x=herb_dif_pd_rr)) + 
  geom_histogram(aes(y=..density..), colour="gray80", fill="white")+
  geom_vline(xintercept = mean(plot_data2$herb_dif_pd_rr), color="red", linetype="dashed", size=1)+
  geom_density(alpha=.2, fill="blue") +
  xlim(-125, 125)+
  labs(x=expression("\u0394PD/year"), y="Density")+
  theme_classic()

b<-ggplot(plot_data2, aes(x=-herb_lost_pd_rr)) + 
  geom_histogram(aes(y=..density..), colour="gray80", fill="white")+
  geom_vline(xintercept = -mean(plot_data2$herb_lost_pd_rr), color="red", linetype="dashed", size=1)+
  geom_density(alpha=.2, fill="gold") +
  xlim(-125, NA)+
  ylim(0, 0.042)+
  labs(x=expression("Lost PD/year"), y="Density")+
  theme_classic()

c<-ggplot(plot_data2, aes(x=herb_gain_pd_rr)) + 
  geom_histogram(aes(y=..density..), colour="gray80", fill="white")+
  geom_vline(xintercept = mean(plot_data2$herb_gain_pd_rr), color="red", linetype="dashed", size=1)+
  geom_density(alpha=.2, fill="green") +
  xlim(NA, 125)+
  ylim(0, 0.042)+
  labs(x=expression("Gained PD/year"), y="Density")+
  theme_classic()

# Put together:
d<-ggarrange(a,                                                 # First row with scatter plot
             ggarrange(b, c, ncol = 2, labels = c("b)", "c)")), # Second row with box and dot plots
             nrow = 2, heights = c(2, 1.3),
             labels = "a)"                                        # Labels of the scatter plot
)

# Save result:
png("Results/Fig1.png",
    res=600, height=4.5,width=4,units="in"); 
d
grid.text("-2.97", x = unit(0.575, "npc"), y = unit(0.98, "npc"), hjust=0, gp=gpar(fontsize=11, fontface="bold", col="red"))
grid.text("-17.76", x = unit(0.3, "npc"), y = unit(0.37, "npc"), hjust=0, gp=gpar(fontsize=11, fontface="bold", col="red"))
grid.text("14.79", x = unit(0.71, "npc"), y = unit(0.37, "npc"), hjust=0, gp=gpar(fontsize=11, fontface="bold", col="red"))
dev.off()


####
# C: H4 and H5
####

# Violin plots for MNTD.ses:

# Aggregate data to plot:
vl<-data.frame(div=c(plot_data2$herb_lost_mntd_ses, 
                     plot_data2$herb_gain_mntd_ses),
               dir=c(rep("Lost\nspecies", nrow(plot_data2)), 
                     rep("Gained\nspecies", nrow(plot_data2))))

# Change order of factors:
vl$dir <- factor(vl$dir, levels = c("Lost\nspecies", "Gained\nspecies"))

# Plot:
b<-ggplot(vl, aes(x=dir, y=div)) + 
  geom_violin()+
  annotate("rect", xmin = 0.25, xmax = 2.75, ymin = -1.96, ymax = 1.96,
           alpha = .2) +
  geom_hline(yintercept=0, linetype="dashed") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1.96),
               geom="pointrange", color="red")+
  ylim(-6,3)+
  labs(y="MNTD.ses", x=" ") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        legend.title=element_blank())

# Violin plot for Dnn.ses:

# Aggregate data to plot:
vl<-data.frame(div=c(plot_data2$herb_dnn_ses),
               type=c(rep("Lost vs. gained\nspecies", nrow(plot_data2)))) #create empty data frame

# Plot:
c<-ggplot(vl, aes(y=div, x=type)) + 
  geom_violin()+
  annotate("rect", xmin = 0.25, xmax = 1.75, ymin = -1.96, ymax = 1.96,
           alpha = .2) +
  geom_hline(yintercept=0, linetype="dashed") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1.96),
               geom="pointrange", color="red")+
  labs(y=expression(D[nn]*".ses"), x=" ") +
  ylim(-6,3)+
  theme_classic()+
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        legend.title=element_blank())

# Save result:
png("Results/Fig3.png",
    res=600, height=3.2,width=4.3,units="in"); 
ggarrange(b, c, labels=c("a)", "b)"), widths = c(2, 1.5), common.legend = TRUE, ncol = 2, nrow = 1)
dev.off()


####
# D: H6
####

# Transform variables:
plot_data$sqrt_herb_lost_pd<-sqrt(plot_data$herb_lost_pd)
plot_data$sqrt_herb_gain_pd<-sqrt(plot_data$herb_gain_pd)

# Get list of response variables:
resp<-c("sqrt_herb_lost_pd", "sqrt_herb_gain_pd", "herb_lost_mpd_ses", "herb_lost_mntd_ses",
        "herb_gain_mpd_ses", "herb_gain_mntd_ses", "herb_dpw_ses", "herb_dnn_ses")

# Create object to save results:
res <- list()
r2s <- list()

# Run loop for each response variable:
for (i in 1:length(resp)){ 
  
  # Select variable:
  plot_data$a<-plot_data[,which(colnames(plot_data)==resp[i])]
  
  # Run Linear Mixed-Effect Model:
  m1<-lme(a ~ 
            scale(dif_tmx_avg) + scale(dif_tmn_avg) + scale(dif_prec_avg) + scale(n_dep) +
            scale(log(plot_size)) + scale(sr_all_bl) + scale(herb_cover) + scale(dif_tree_cover) + 
            scale(time_resurvey) + scale(bl_tmx_avg) + scale(bl_tmn_avg) + scale(bl_prec_avg) + scale(soil_ph), 
          random = ~1|study, method="REML", data=plot_data, na.action = na.omit)

  # Print AIC:
  print(AIC(m1))
  
  # Save Marginal and Conditional R2s:
  r2s[[i]]<-r.squaredGLMM(m1)
  
  # Save model summary:
  m1<-summary(m1)
  
  # Save output:
  res[[i]]<-m1$tTable
}

# Aggregate results:
res<-do.call(cbind, res)
r2s<-round(do.call(rbind, r2s), 2)

# Save results:
write.table(res, "data/models_lme.csv")
write.table(r2s, "data/r2s_lme.csv")

# Clean-up environment:
rm(list = ls())