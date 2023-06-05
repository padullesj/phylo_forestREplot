
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

# Load plot data:
plot_data<-read.table("data/plot_data_forestreplot.csv")

# Add "study" to vegetation data:
plot_data$study<-substr(plot_data$plotID, 1, 6) #add field

# Calculate rates of change of environmental variables:
plot_data$rate_dif_tmx_avg<-(plot_data$dif_tmx_avg/plot_data$time_resurvey)
plot_data$rate_dif_tmn_avg<-(plot_data$dif_tmn_avg/plot_data$time_resurvey)
plot_data$rate_dif_prec_avg<-(plot_data$dif_prec_avg/plot_data$time_resurvey)
plot_data$rate_n_dep<-(plot_data$n_dep/plot_data$time_resurvey)

# Subset plots with resurvey time >= 20 years for H1, H2 & H3:
plot_data2<-subset(plot_data, time_resurvey>=20)

# Subset plots re-surveyed after 2000 for H1, H2 & H3:
plot_data2<-subset(plot_data2, year_final>2000)
plot_data2<-plot_data2[order(row.names(plot_data2)), ]

# Subset rows with non NA:
plot_data2<-subset(plot_data2, !is.infinite(plot_data2$herb_delta_pd_obs))
plot_data2<-subset(plot_data2, !is.na(plot_data2$herb_delta_pd_obs))

##
# B: H1
##

# Calculate Cohen's D:
cohensD(plot_data2$herb_delta_pd_obs*100, mu=0) #test H1
cohensD(plot_data2$herb_delta_pd_ses, mu=0) #test H1

# Prepare table for histogram of Figure 1:
vl<-data.frame(div=c(plot_data2$herb_delta_pd_obs*100, 
                     plot_data2$herb_delta_pd_ses),
               dir=c(rep("RR_PD x100", nrow(plot_data2)), rep("RR_PD.ses", nrow(plot_data2)))) #create empty data frame
vl$dir <- factor(vl$dir, levels = c("RR_PD x100", "RR_PD.ses")) #set order of factors

# Get mean by group:
vl2 <- vl %>%
  group_by(dir) %>%
  dplyr::summarise(mean.pd = mean(div, na.rm=T))

# Plot histograms for Figure 1:
d<-ggplot(vl, aes(x = div, y = ..count.. + 1, fill = ..x..)) +
  facet_grid(. ~ dir) +
  geom_histogram(binwidth = .2, pad = TRUE) +
  scale_fill_gradient2(name = "val", low = "darkred", mid="gray99", high = "darkblue")+
  geom_vline(data = vl2, mapping = aes(xintercept = mean.pd), color="red", linetype="dashed", size=.5) +
  labs(x=expression(" "), y="Frequency")+
  theme_classic()+
  ylim(0, 210)+
  theme(axis.title.y= element_blank(),
        strip.background = element_blank(),
        legend.position = "none")


####
# C: H4 and H5
####

#Tests for differences from zero:
cohensD(plot_data2$herb_lost_pd_ses, mu=0)
cohensD(plot_data2$herb_lost_mntd_ses, mu=0)
cohensD(plot_data2$herb_lost_dnn_ses, mu=0)

cohensD(plot_data2$herb_gain_pd_ses, mu=0)
cohensD(plot_data2$herb_gain_mntd_ses, mu=0)
cohensD(plot_data2$herb_gain_dnn_ses, mu=0)

# Violin plots for PR metrics:

# Aggregate data to plot:
vl<-data.frame(div=c(plot_data2$herb_lost_mntd_ses,
                     plot_data2$herb_gain_mntd_ses,
                     plot_data2$herb_lost_dnn_ses,
                     plot_data2$herb_gain_dnn_ses),
               dir=c(rep("MNTD.ses", nrow(plot_data2)), 
                     rep("MNTD.ses", nrow(plot_data2)),
                     rep("Dnn.ses", nrow(plot_data2)),
                     rep("Dnn.ses", nrow(plot_data2))),
               grp=c(rep("Lost species", nrow(plot_data2)), 
                     rep("Gained species", nrow(plot_data2)),
                     rep("Lost species", nrow(plot_data2)),
                     rep("Gained species", nrow(plot_data2)))) #create empty data frame
vl$dir <- factor(vl$dir, levels = c("MNTD.ses", "Dnn.ses")) #changer order of factors
vl$grp <- factor(vl$grp, levels = c("Lost species", "Gained species")) #changer order of factors
vl$tot <- "value" #add column

# To add the right labels into the plots:
pred_names <- c(
  `MNTD.ses` = "MNTD.ses",
  `Dnn.ses` = "D[nn]*.ses"
)

# Plot:
b<-ggplot(vl, aes(x=tot, y=div)) + 
  ggh4x::facet_nested(. ~ grp + dir, labeller = labeller(
    dir = as_labeller(pred_names, label_parsed)
  )) +
  geom_violin()+
  geom_hline(yintercept=0, linetype="dashed") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1.96), 
               geom="pointrange", color="red")+
  labs(y=" ", x=" ") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank())


####
# D: H6
####

# Models for RR_PD

#Subset rows with non NA
plot_data3<-subset(plot_data, !is.infinite(plot_data$herb_delta_pd_obs))
plot_data3<-subset(plot_data3, !is.na(plot_data3$herb_delta_pd_obs))
plot_data3$herb_delta_pd_obs2<-plot_data3$herb_delta_pd_obs*100

# Get names of response variables:
resp<-c("herb_delta_pd_obs2", "herb_delta_pd_ses")

# Create object to save results and run loop:
res <- list()
r2s <- list()

# Run loop
for (i in 1:length(resp)){ 
  
  # Select column:
  plot_data3$a<-plot_data3[,which(colnames(plot_data3)==resp[i])]
  
  # Run full model
  m1<-lme(a ~ 
            scale(rate_dif_tmx_avg) + scale(rate_dif_tmn_avg) + scale(rate_dif_prec_avg) + scale(rate_n_dep) +
            scale(log(plot_size)) + scale(sr_all_bl) + scale(herb_cover) + 
            scale(bl_tmx_avg) + scale(bl_tmn_avg) + scale(bl_prec_avg) + 
            scale(dif_tree_cover) + scale(soil_ph), 
          random = ~1|study, method="REML", data=plot_data3, na.action = na.omit)

  # Print AIC:
  print(AIC(m1))
  
  # Save R-squared:
  r2s[[i]]<-r.squaredGLMM(m1)
  
  # Prepare table with coefficients:
  m1<-as.data.frame(summary(m1)$tTable)
  m1$type<-resp[i]
  m1$var<-rownames(m1)
  
  # Set confidence interval:
  m1$low<-m1$Value - (m1$Std.Error* 1.959)
  m1$up<-m1$Value + (m1$Std.Error* 1.959)
  
  # Define groups of variables based on the sign of the effect:
  m1$col <- ifelse((m1$low < 0) & (m1$up < 0),"Negative", ifelse((m1$low > 0) & (m1$up > 0),"Positive", ifelse("Not significant")))
  m1$col[is.na(m1$col)] <- "Not significant"
  
  # Save output:
  res[[i]]<-m1
}

# Aggregate results:
res<-do.call(rbind, res)
r2s<-round(do.call(rbind, r2s), 2)

# Set order of the different factors:
res<-res[res$var != "(Intercept)", ] #remove intercept term
res$var <- plyr::revalue(res$var, c("scale(rate_dif_tmx_avg)" = "\u394 Max. summer temperature/year", 
                                    "scale(rate_dif_tmn_avg)" = "\u394 Min. winter temperature/year",
                                    "scale(rate_dif_prec_avg)" = "\u394 Annual precipitation/year",
                                    "scale(rate_n_dep)" = "N deposition/year",
                                    "scale(log(plot_size))" = "Plot area (log-transformed)",
                                    "scale(sr_all_bl)" = "Species richness (bl.)",
                                    "scale(herb_cover)" = "% herb-layer cover (bl.)",
                                    "scale(bl_tmx_avg)" = "Max. summer temperature (bl.)",
                                    "scale(bl_tmn_avg)" = "Min. winter temperature (bl.)",
                                    "scale(bl_prec_avg)" = "Annual precipitation (bl.)",
                                    "scale(dif_tree_cover)" = "\u394 % canopy cover",
                                    "scale(soil_ph)" = "Soil pH (2021)"))

# Rename variables:
res$type <- plyr::revalue(res$type, c("herb_delta_pd_obs2" = "RR_PD", 
                                      "herb_delta_pd_ses" = "RR_PD.ses"))

# Reorder factor levels:
res$type <- factor(res$type, levels = c("RR_PD", "RR_PD.ses")) #changer order of factors
res$var <- factor(res$var, levels = rev(c("\u394 Max. summer temperature/year", "\u394 Min. winter temperature/year","\u394 Annual precipitation/year",
                                          "N deposition/year","Plot area (log-transformed)", "Species richness (bl.)", "% herb-layer cover (bl.)",
                                          "Max. summer temperature (bl.)", "Min. winter temperature (bl.)", "Annual precipitation (bl.)",
                                          "\u394 % canopy cover", "Soil pH (2021)"))) #changer order of factors
res$col <- factor(res$col, levels = c("Not significant", "Negative", "Positive")) #changer order of factors

#Plot:
p1<-ggplot(res, aes(x=var, y=Value, ymin=low, ymax=up, shape=col, fill=col, col=col)) +  
  facet_grid(. ~ type) +
  geom_errorbar(size=0.3, width = 0.2) + 
  geom_point(size=3.5) +
  geom_vline(xintercept=8.5, linetype="dotted") +
  scale_shape_manual(values=c(21, 21, 21))+
  scale_color_manual(values=c("black", "darkorange", "blue")) +
  scale_fill_manual(values=c("white", "darkorange", "blue")) +
  geom_hline(yintercept=0, lty=2, size=0.1) +
  coord_flip() +
  xlab("") + ylab("Estimates (± 95% confidence interval)") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(angle=45, vjust=.5),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(color="black", hjust = 0.5, size=8, face="bold"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"))

# Put together plots from H1 and H6:
figure<-ggarrange(d, p1,
                  heights = c(0.6, 1),
                  ncol = 1, nrow = 2, align = "v")

# Save results:
png("Fig1.png",
    res=600, height=4.7,width=6,units="in"); 
figure

grid.text("Frequency", x = unit(0.28, "npc"), y = unit(0.825, "npc"), rot=90,
          gp=gpar(fontsize=10))

grid.text("a", x = unit(0.33, "npc"), y = unit(0.98, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))
grid.text("b", x = unit(0.68, "npc"), y = unit(0.98, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))
grid.text("c", x = unit(0.33, "npc"), y = unit(0.65, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))
grid.text("d", x = unit(0.68, "npc"), y = unit(0.65, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))

grid.text("-0.19", x = unit(0.445, "npc"), y = unit(0.92, "npc"), hjust=0, gp=gpar(fontsize=10, fontface="bold", col="red"))
grid.text("0.26", x = unit(0.84, "npc"), y = unit(0.92, "npc"), hjust=0, gp=gpar(fontsize=10, fontface="bold", col="red"))

grid.text(expression(italic(d)*"=0.13"), x = unit(0.61, "npc"), y = unit(0.73, "npc"), gp=gpar(fontsize=10))
grid.text(expression(italic(d)*"=0.23"), x = unit(0.93, "npc"), y = unit(0.73, "npc"), gp=gpar(fontsize=10))

dev.off()


# Models for PR metrics:

# Get names of response variables:
resp<-c("herb_lost_sr", "herb_gain_sr",
        "herb_lost_pd_ses", "herb_gain_pd_ses", 
        "herb_lost_mpd_ses", "herb_lost_mntd_ses",
        "herb_gain_mpd_ses", "herb_gain_mntd_ses", 
        "herb_lost_dpw_ses", "herb_lost_dnn_ses",
        "herb_gain_dpw_ses", "herb_gain_dnn_ses")

# Split between lost and gained species:
grp<-c("Lost species", "Gained species", 
       "Lost species", "Gained species", 
       "Lost species", "Lost species", 
       "Gained species", "Gained species", 
       "Lost species", "Lost species",
       "Gained species", "Gained species")

# Create object to save results and run loop:
res <- list()
r2s <- list()
for (i in 1:length(resp)){ 
  
  # Select column:
  plot_data$a<-plot_data[,which(colnames(plot_data)==resp[i])]
  
  #Run model:
  m1<-lme(a ~ 
            scale(dif_tmx_avg) + scale(dif_tmn_avg) + scale(dif_prec_avg) + scale(n_dep) +
            scale(log(plot_size)) + scale(sr_all_bl) + scale(herb_cover) + 
            scale(bl_tmx_avg) + scale(bl_tmn_avg) + scale(bl_prec_avg) + 
            scale(dif_tree_cover) + scale(soil_ph) + scale(time_resurvey), 
          random = ~1|study, method="REML", data=plot_data, na.action = na.omit)
  
  # Print AIC:
  print(AIC(m1))
  
  # Save R-squared:
  r2s[[i]]<-r.squaredGLMM(m1)
  
  # Prepare table with coefficients:
  m1<-as.data.frame(summary(m1)$tTable)
  m1$type<-resp[i]
  m1$grp<-grp[i]
  m1$var<-rownames(m1)
  
  # Set confidence interval:
  m1$low<-m1$Value - (m1$Std.Error* 1.959)
  m1$up<-m1$Value + (m1$Std.Error* 1.959)
  
  # Define groups of variables based on the sign of the effect:
  m1$col <- ifelse((m1$low < 0) & (m1$up < 0),"Negative", ifelse((m1$low > 0) & (m1$up > 0),"Positive", ifelse("Not significant")))
  m1$col[is.na(m1$col)] <- "Not significant"
  
  # Save output:
  res[[i]]<-m1
}

# Aggregate results:
res<-do.call(rbind, res)
r2s<-round(do.call(rbind, r2s), 2)

# Set order of the different factors:
res<-res[res$var != "(Intercept)", ]
res$var <- plyr::revalue(res$var, c("scale(dif_tmx_avg)" = "\u394 Max. summer temperature", 
                                    "scale(dif_tmn_avg)" = "\u394 Min. winter temperature",
                                    "scale(dif_prec_avg)" = "\u394 Annual precipitation",
                                    "scale(n_dep)" = "N deposition",
                                    "scale(log(plot_size))" = "Plot area (log-transformed)",
                                    "scale(sr_all_bl)" = "Species richness (bl.)",
                                    "scale(herb_cover)" = "% herb-layer cover (bl.)",
                                    "scale(bl_tmx_avg)" = "Max. summer temperature (bl.)",
                                    "scale(bl_tmn_avg)" = "Min. winter temperature (bl.)",
                                    "scale(bl_prec_avg)" = "Annual precipitation (bl.)",
                                    "scale(dif_tree_cover)" = "\u394 % canopy cover",
                                    "scale(soil_ph)" = "Soil pH (2021)",
                                    "scale(time_resurvey)" = "Time between surveys"))

# Subset variables of interest:
res1<-res #create a copy
res1<-res1[res1$type %in% c("herb_lost_mntd_ses", "herb_gain_mntd_ses",
                            "herb_lost_dnn_ses", "herb_gain_dnn_ses"), ]

# Rename factors:
res1$type <- plyr::revalue(res1$type, c("herb_lost_mntd_ses" = "MNTD.ses",
                                        "herb_gain_mntd_ses" = "MNTD.ses",
                                        "herb_lost_dnn_ses" = "Dnn.ses",
                                        "herb_gain_dnn_ses" = "Dnn.ses"))

# Reorder factor levels:
res1$type <- factor(res1$type, levels = c("MNTD.ses", "Dnn.ses")) #changer order of factors
res1$var <- factor(res1$var, levels = rev(c("\u394 Max. summer temperature", "\u394 Min. winter temperature","\u394 Annual precipitation",
                                            "N deposition","Plot area (log-transformed)", "Species richness (bl.)", "% herb-layer cover (bl.)",
                                            "Max. summer temperature (bl.)", "Min. winter temperature (bl.)", "Annual precipitation (bl.)",
                                            "\u394 % canopy cover", "Soil pH (2021)", "Time between surveys"))) #changer order of factors
res1$col <- factor(res1$col, levels = c("Not significant", "Negative", "Positive")) #changer order of factors
res1$grp <- factor(res1$grp, levels = c("Lost species", "Gained species")) #changer order of factors

# Plot:
p2<-ggplot(res1, aes(x=var, y=Value, ymin=low, ymax=up, shape=col, fill=col, col=col)) +  
  ggh4x::facet_nested(. ~ grp + type) +
  geom_errorbar(size=0.3, width = 0.2) + 
  geom_point(size=3.5) +
  geom_vline(xintercept=9.5, linetype="dotted") +
  scale_shape_manual(values=c(21, 21, 21))+
  scale_color_manual(values=c("black", "darkorange", "blue")) +
  scale_fill_manual(values=c("white", "darkorange", "blue")) +
  geom_hline(yintercept=0, lty=2, size=0.1) +
  coord_flip() +
  xlab("") + ylab("Estimates (± 95% confidence interval)") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(angle=45, vjust=.5),
        #panel.spacing.x = unit(4, "mm"),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.title = element_text(color="black", hjust = 0.5, size=8, face="bold"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"))

# Put together plots from H4, H5, and H6:
figure <- ggarrange(b, p2,
                    heights = c(0.6, 1),
                    ncol = 1, nrow = 2, align = "v")

# Save final plot:
png("Fig3.png",
    res=600, height=6,width=7,units="in"); 
figure

grid.text("Standardized effect sizes", x = unit(0.25, "npc"), y = unit(0.79, "npc"), rot=90,
          gp=gpar(fontsize=10))

grid.text("a", x = unit(0.27, "npc"), y = unit(0.98, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))
grid.text("b", x = unit(0.64, "npc"), y = unit(0.98, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))
grid.text("c", x = unit(0.27, "npc"), y = unit(0.63, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))
grid.text("d", x = unit(0.64, "npc"), y = unit(0.63, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))

grid.lines(x=c(0.30, .63), y = c(0.9,0.9), default.units='npc' )
grid.lines(x=c(0.65, .97), y = c(0.9,0.9), default.units='npc' )

grid.lines(x=c(0.30, .63), y = c(0.94,0.94), default.units='npc' )
grid.lines(x=c(0.65, .97), y = c(0.94,0.94), default.units='npc' )

grid.text(expression(italic(d)*"=0.31"), x = unit(0.41, "npc"), y = unit(0.88, "npc"), gp=gpar(fontsize=10))
grid.text(expression(italic(d)*"=0.07"), x = unit(0.59, "npc"), y = unit(0.88, "npc"), gp=gpar(fontsize=10))

grid.text(expression(italic(d)*"=0.06"), x = unit(0.755, "npc"), y = unit(0.88, "npc"), gp=gpar(fontsize=10))
grid.text(expression(italic(d)*"=0.18"), x = unit(0.93, "npc"), y = unit(0.88, "npc"), gp=gpar(fontsize=10))

dev.off()


# Clean-up environment:
rm(list = ls())