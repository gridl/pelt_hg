#############################################################################
# Otter fur THg to organ THg otter sampling spot optimization
# Written in R Version 3.5.0
#############################################################################
# Load data
compiled = read.csv("pelts_compiled.csv")

# Library
library(ggplot2)
library(Metrics)
library(rgdal)
library(ggpubr)
library(dplyr)
library(tidyr)
library(rgeos)
library(smoothr)
library(sp)
library(spdep)
library(gridExtra)
library(adespatial)

#Eccles et al. (2017) fur THg to organ THg predictive model eqns
#Fur to brain: y = 0.15x + 0
#Fur to kidney: y = 0.62x + 0
#Fur to liver: y = 0.70x + 0
#Fur to muscle: y = 0.46x + 0

#######################################
# Import Shapefiles
# Composite Pelt 
comp_shp = readOGR(dsn = ".", layer = "Comp_Outline_2")
# Convert shapefile to a ggplot ready data frame
comp_shp_df <- fortify(comp_shp, region = "id")
comp_shp$id <- rownames(comp_shp@data)
comp_shp_df <- left_join(comp_shp_df, comp_shp@data,mby = "id")

#########################################################################
#### Compiled TC ####
# Brain
compiled$TC_est_brain= 0.15*compiled$THg_TC
compiled$TC_resid_brain=compiled$TC_est_brain-compiled$Brain
compiled$TC_perc_brain= abs(((compiled$TC_est_brain-compiled$Brain)/(compiled$Brain))*100)

# Liver
compiled$TC_est_liver= 0.70*compiled$THg_TC
compiled$TC_resid_liver=compiled$TC_est_liver-compiled$Liver
compiled$TC_perc_liver= abs(((compiled$TC_est_liver-compiled$Liver)/(compiled$Liver))*100)

# Kidney
compiled$TC_est_kidney= 0.46*compiled$THg_TC
compiled$TC_resid_kidney=compiled$TC_est_kidney-compiled$Kidney
compiled$TC_perc_kidney= abs(((compiled$TC_est_kidney-compiled$Kidney)/(compiled$Kidney))*100)

# Muscle
compiled$TC_est_muscle = 0.62*compiled$THg_TC
compiled$TC_resid_muscle = compiled$TC_est_muscle-compiled$Muscle
compiled$TC_perc_muscle = abs(((compiled$TC_est_muscle-compiled$Muscle)/(compiled$Muscle))*100)
#################################################################
# Average error TC
compiled$TC_comp_average_error = (compiled$TC_perc_brain/100*.4 + compiled$TC_perc_liver/100*.2 + 
                                 compiled$TC_perc_kidney/100*.2 + compiled$TC_perc_muscle/100*.2)*100

adj_x_compiled = aggregate(compiled$X_adj, list(compiled$factor), mean)
adj_y_compiled = aggregate(compiled$Y_adj, list(compiled$factor), mean)
TC_compiled = aggregate(compiled$TC_comp_average_error, list(compiled$factor), mean)
compiled_pelt_TC = as.data.frame(cbind(adj_x_compiled[,2], adj_y_compiled[,2],TC_compiled[,2]))
colnames(compiled_pelt_TC) <- c("x", "y", "TC")

#compiled_pelt_TC= subset(compiled_pelt_TC, TC < 40)

TC_error = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = compiled_pelt_TC, aes(x = x, y = y, color = TC), size=5)+
  scale_colour_gradient(low = "white", high = "red", limits =c(0,115))+
  labs(colour="Percent Error")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

#Set up spatial neighbourhood
comp_sppnt = SpatialPoints(cbind(compiled_pelt_TC$x, compiled_pelt_TC$y))
comp_nb <- chooseCN(coordinates(comp_sppnt), type = 6, k = 12, plot.nb = FALSE) 
comp_lw <- nb2listw(comp_nb, style = 'W', zero.policy = TRUE)

# Gi*
compiled_pelt_TC$Gi_TC = as.numeric(localG(compiled_pelt_TC$TC, listw = comp_lw))

# Gi* p-value
#Z score values for levels of statistical significance:
# 90% significant: >= 1.645
# 95% significant: >= 1.960
# 99% significant: >= 2.576
# 99.9% significant: >= 3.291

# Bin gi* scores based on whether it is past the critical cut off
compiled_pelt_TC$Gi_pval_TC = as.factor(compiled_pelt_TC$Gi_TC >= 1.96 |  compiled_pelt_TC$Gi_TC <= -1.96) 

#TC THg Gi*
# Plot outline and points
hs_comp_TC = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = compiled_pelt_TC, aes(x = x, y = y, color = compiled_pelt_TC$Gi_TC), size=5, shape=19) +
  geom_point(data = compiled_pelt_TC, aes(x = x, y = y, shape = factor(compiled_pelt_TC$Gi_pval_TC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi*score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-6, 6),
                        guide = guide_colourbar(title = "Gi*score")) +
  coord_quickmap()+
  theme_void()
########################################################################
#### Compiled UC ####
# Brain
compiled$UC_est_brain= 0.15*compiled$THg_UC
compiled$UC_resid_brain=compiled$UC_est_brain-compiled$Brain
compiled$UC_perc_brain= abs(((compiled$UC_est_brain-compiled$Brain)/(compiled$Brain))*100)

# Liver
compiled$UC_est_liver= 0.70*compiled$THg_UC
compiled$UC_resid_liver=compiled$UC_est_liver-compiled$Liver
compiled$UC_perc_liver= abs(((compiled$UC_est_liver-compiled$Liver)/(compiled$Liver))*100)

# Kidney
compiled$UC_est_kidney= 0.46*compiled$THg_UC
compiled$UC_resid_kidney=compiled$UC_est_kidney-compiled$Kidney
compiled$UC_perc_kidney= abs(((compiled$UC_est_kidney-compiled$Kidney)/(compiled$Kidney))*100)

# Muscle
compiled$UC_est_muscle = 0.62*compiled$THg_UC
compiled$UC_resid_muscle = compiled$UC_est_muscle-compiled$Muscle
compiled$UC_perc_muscle = abs(((compiled$UC_est_muscle-compiled$Muscle)/(compiled$Muscle))*100)

# Average error UC
compiled$UC_comp_average_error = (compiled$UC_perc_brain/100*.4 + compiled$UC_perc_liver/100*.2 + 
                                    compiled$UC_perc_kidney/100*.2 + compiled$UC_perc_muscle/100*.2)*100

adj_x_compiled = aggregate(compiled$X_adj, list(compiled$factor), mean)
adj_y_compiled = aggregate(compiled$Y_adj, list(compiled$factor), mean)
UC_compiled = aggregate(compiled$UC_comp_average_error, list(compiled$factor), mean)
compiled_pelt_UC = as.data.frame(cbind(adj_x_compiled[,2], adj_y_compiled[,2],UC_compiled[,2]))
colnames(compiled_pelt_UC) <- c("x", "y", "UC")

UC_error = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = compiled_pelt_UC, aes(x = x, y = y, color = UC), size=5)+
  scale_colour_gradient(low = "white", high = "red", limits =c(0,115))+
  labs(colour="Percent Error")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

#Set up spatial neighbourhood
comp_sppnt = SpatialPoints(cbind(compiled_pelt_UC$x, compiled_pelt_UC$y))
comp_nb <- chooseCN(coordinates(comp_sppnt), type = 6, k = 12, plot.nb = FALSE) 
comp_lw <- nb2listw(comp_nb, style = 'W', zero.policy = TRUE)

# Gi*
compiled_pelt_UC$Gi_UC = as.numeric(localG(compiled_pelt_UC$UC, listw = comp_lw))

# Gi* p-value
#Z score values for levels of statistical significance:
# 90% significant: >= 1.645
# 95% significant: >= 1.960
# 99% significant: >= 2.576
# 99.9% significant: >= 3.291

# Bin gi* scores based on whether it is past the critical cut off
compiled_pelt_UC$Gi_pval_UC = as.factor(compiled_pelt_UC$Gi_UC >= 1.96 |  compiled_pelt_UC$Gi_UC <= -1.96) 

#All THg Gi*
# Plot outline and points
hs_comp_UC = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = compiled_pelt_UC, aes(x = x, y = y, color = compiled_pelt_UC$Gi_UC), size=5, shape=19) +
  geom_point(data = compiled_pelt_UC, aes(x = x, y = y, shape = factor(compiled_pelt_UC$Gi_pval_UC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi*score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-6, 6),
                        guide = guide_colourbar(title = "Gi*score")) +
  coord_quickmap()+
  theme_void()

##############################################################
p1 = ggarrange(hs_comp_TC, hs_comp_UC,
          labels = c("TC Hotspot", "UC Hotspot"),
          vjust = 1,
          hjust = -.5,
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = "right")

#Plot figures with dpi=300
save_plot("composite_error_hs.tif", p1, width = 15, height = 15, dpi = 300)

p2 = ggarrange(TC_error, UC_error,
          labels = c("Average Error TC", "Average Error UC", "TC Hotspot", "UC Hotspot"),
          vjust = 1,
          hjust = -.5,
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = "right")

#Plot figures with dpi=300
save_plot("composite_error.tif", p2, width = 15, height = 15, dpi = 300)

#############################################################
#### Average TC and UC ####
compiled$all_comp_average_error = ((compiled$UC_perc_brain/100*.2 + compiled$UC_perc_liver/100*.1 + 
                                      compiled$UC_perc_kidney/100 *.1+ compiled$UC_perc_muscle/100*.1+ 
                                      compiled$TC_perc_brain/100 *.2 + compiled$TC_perc_liver/100*.1 + 
                                      compiled$TC_perc_kidney/100*.1 + compiled$TC_perc_muscle/100*.1))*100

adj_x_compiled = aggregate(compiled$X_adj, list(compiled$factor), mean)
adj_y_compiled = aggregate(compiled$Y_adj, list(compiled$factor), mean)
all_compiled = aggregate(compiled$all_comp_average_error, list(compiled$factor), mean)
compiled_pelt_all = as.data.frame(cbind(adj_x_compiled[,2], adj_y_compiled[,2],all_compiled[,2]))
colnames(compiled_pelt_all) <- c("x", "y", "all")

all = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = compiled_pelt_all, aes(x = x, y = y, color = all), size=4)+
  scale_colour_gradient(low = "white", high = "red")+
  labs(colour="Percent Error")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

### Compiled Pelt Hot Spot Analysis ####
#Set up spatial neighbourhood
comp_sppnt = SpatialPoints(cbind(compiled_pelt_all$x, compiled_pelt_all$y))
comp_nb <- chooseCN(coordinates(comp_sppnt), type = 6, k = 12, plot.nb = FALSE) 
comp_lw <- nb2listw(comp_nb, style = 'W', zero.policy = TRUE)

# Gi*
compiled_pelt_all$Gi_all = as.numeric(localG(compiled_pelt_all$all, listw = comp_lw))

# Gi* p-value
#Z score values for levels of statistical significance:
# 90% significant: >= 1.645
# 95% significant: >= 1.960
# 99% significant: >= 2.576
# 99.9% significant: >= 3.291

# Bin gi* scores based on whether it is past the critical cut off
compiled_pelt_all$Gi_pval_all = as.factor(compiled_pelt_all$Gi_all >= 1.96 |  compiled_pelt_all$Gi_all <= -1.96) 

#All THg Gi*
# Plot outline and points
hs_comp_all = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = compiled_pelt, aes(x = x, y = y, color = compiled_pelt$Gi_TC), size=5, shape=19) +
  geom_point(data = compiled_pelt, aes(x = x, y = y, shape = factor(compiled_pelt$Gi_pval_TC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi*score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-6, 6),
                        guide = guide_colourbar(title = "Gi*score")) +
  coord_quickmap()+
  theme_void()
##############################################################
plot1=ggarrange(all, hs_comp_all,
          labels = c("          Average Error", "Average Error Hot Spots"),
          vjust = 1,
          hjust = -0.15,
          ncol = 2, nrow = 1,
          common.legend = FALSE,
          legend = "right")

#Plot figures with dpi=300
save_plot("composite_error_hs.tif", plot1, width = 20, height = 15, dpi = 300)
