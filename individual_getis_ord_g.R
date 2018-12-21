#############################################################################
# Cluster analysis of pelt Hg using Getis and Ord's Gi*
# Written in R Version 3.5.0
#############################################################################
# Load data
pelt1 = read.csv("PELT1.csv")
pelt2 = read.csv("PELT2.csv")
pelt3 = read.csv("PELT3.csv")
pelt4 = read.csv("PELT4.csv")

# Load Libraries
library(spdep)
library(rgdal)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(adespatial)

##############################
### Pelt1 ###
#Set up spatial neighbourhood
p1_sppnt = SpatialPoints(cbind(pelt1$X_coord, pelt1$Y_coord))
p1_nb <- chooseCN(coordinates(p1_sppnt), type = 6, k = 12, plot.nb = FALSE) 
p1_lw <- nb2listw(p1_nb, style = 'W', zero.policy = TRUE)

# Gi*
pelt1$Gi_TC = as.numeric(localG(pelt1$THg_TC, listw = p1_lw))
pelt1$Gi_UC = as.numeric(localG(pelt1$THg_UC, listw = p1_lw))

# Gi* p-value
#Z score values for levels of statistical significance:
# 90% significant: >= 1.645
# 95% significant: >= 1.960
# 99% significant: >= 2.576
# 99.9% significant: >= 3.291

# Bin gi* scores based on whether it is past the critical cut off
pelt1$Gi_pval_TC = as.factor(pelt1$Gi_TC >= 1.96 |  pelt1$Gi_TC <= -1.96) 
pelt1$Gi_pval_UC <- as.factor(pelt1$Gi_UC >= 1.96 |  pelt1$Gi_UC <= -1.96) 

# Plot 1
#TC THg Gi* 
# Importshapefile
p1_shp = readOGR(dsn = ".", layer = "p1_Outline")
# Convert shapefile to a ggplot ready data frame
p1_shp_df <- fortify(p1_shp, region = "id")
p1_shp$id <- rownames(p1_shp@data)
p1_shp_df <- left_join(p1_shp_df, p1_shp@data,mby = "id")
pelt1 = as.data.frame(pelt1)

# Plot outline and points
p1_TC = ggplot(p1_shp_df, aes(x = p1_shp_df$long, y = p1_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$Gi_TC), size=5, shape=19) +
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, shape = factor(pelt1$Gi_pval_TC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 5),
                  guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

p1_UC = ggplot(p1_shp_df, aes(x = p1_shp_df$long, y = p1_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$Gi_UC), size=5, shape=19) +
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, shape = factor(pelt1$Gi_pval_UC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 5),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

################################################
### Pelt2 ###
#Set up spatial neighbourhood
p2_sppnt = SpatialPoints(cbind(pelt2$X_coord, pelt2$Y_coord))
p2_nb <- chooseCN(coordinates(p2_sppnt), type = 6, k = 12, plot.nb = FALSE) 
p2_lw <- nb2listw(p2_nb, style = 'W', zero.policy = TRUE)

# Gi*
pelt2$Gi_TC = as.numeric(localG(pelt2$THg_TC, listw = p2_lw))
pelt2$Gi_UC = as.numeric(localG(pelt2$THg_UC, listw = p2_lw))

# Bin gi* scores based on whether it is past the critical cut off
pelt2$Gi_pval_TC = as.factor(pelt2$Gi_TC >= 1.96 |  pelt2$Gi_TC <= -1.96) 
pelt2$Gi_pval_UC <- as.factor(pelt2$Gi_UC >= 1.96 |  pelt2$Gi_UC <= -1.96)

# Plot 2
#TC THg Gi* 
# Importshapefile
p2_shp = readOGR(dsn = ".", layer = "p2_Outline")
# Convert shapefile to a ggplot ready data frame
p2_shp_df <- fortify(p2_shp, region = "id")
p2_shp$id <- rownames(p2_shp@data)
p2_shp_df <- left_join(p2_shp_df, p2_shp@data,mby = "id")
pelt2 = as.data.frame(pelt2)

p2_TC = ggplot(p2_shp_df, aes(x = p2_shp_df$long, y = p2_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, color = pelt2$Gi_TC), size=5, shape=19) +
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, shape = factor(pelt2$Gi_pval_TC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 5),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

#UC THg Gi* 
p2_UC = ggplot(p2_shp_df, aes(x = p2_shp_df$long, y = p2_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, color = pelt2$Gi_UC), size=5, shape=19) +
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, shape = factor(pelt2$Gi_pval_UC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 5),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

############################################################
### Pelt3 ###
#Set up spatial neighbourhood
p3_sppnt = SpatialPoints(cbind(pelt3$X_coord, pelt3$Y_coord))
p3_nb <- chooseCN(coordinates(p3_sppnt), type = 6, k = 12, plot.nb = FALSE) 
p3_lw <- nb2listw(p3_nb, style = 'W', zero.policy = TRUE)

# Gi*
pelt3$Gi_TC = as.numeric(localG(pelt3$THg_TC, listw = p3_lw))
pelt3$Gi_UC = as.numeric(localG(pelt3$THg_UC, listw = p3_lw))

# Bin gi* scores based on whether it is past the critical cut off
pelt3$Gi_pval_TC = as.factor(pelt3$Gi_TC >= 1.96 |  pelt3$Gi_TC <= -1.96) 
pelt3$Gi_pval_UC <- as.factor(pelt3$Gi_UC >= 1.96 |  pelt3$Gi_UC <= -1.96) 

# Plot 3
# TC THg Gi* 
# Importshapefile
p3_shp = readOGR(dsn = ".", layer = "p3_Outline")
# Convert shapefile to a ggplot ready data frame
p3_shp_df <- fortify(p3_shp, region = "id")
p3_shp$id <- rownames(p3_shp@data)
p3_shp_df <- left_join(p3_shp_df, p3_shp@data,mby = "id")
pelt3 = as.data.frame(pelt3)

p3_TC = ggplot(p3_shp_df, aes(x = p3_shp_df$long, y = p3_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, color = pelt3$Gi_TC), size=5, shape=19) +
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, shape = factor(pelt3$Gi_pval_TC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                                          labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 6),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

#UC THg Gi* 
p3_UC = ggplot(p3_shp_df, aes(x = p3_shp_df$long, y = p3_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, color = pelt3$Gi_UC), size=5, shape=19) +
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, shape = factor(pelt3$Gi_pval_UC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                     labels = c("Not Significant", "Significant"))+
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 7),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

#########################################################
### Pelt4 ###
#Set up spatial neighbourhood
p4_sppnt = SpatialPoints(cbind(pelt4$X_coord, pelt4$Y_coord))
p4_nb <- chooseCN(coordinates(p4_sppnt), type = 6, k = 12, plot.nb = FALSE) 
p4_lw <- nb2listw(p4_nb, style = 'W', zero.policy = TRUE)

# Gi*
pelt4$Gi_TC = as.numeric(localG(pelt4$THg_TC, listw = p4_lw))
pelt4$Gi_UC = as.numeric(localG(pelt4$THg_UC, listw = p4_lw))

# Bin gi* scores based on whether it is past the critical cut off
pelt4$Gi_pval_TC = as.factor(pelt4$Gi_TC >= 1.96 |  pelt4$Gi_TC <= -1.96) 
pelt4$Gi_pval_UC <- as.factor(pelt4$Gi_UC >= 1.96 |  pelt4$Gi_UC <= -1.96) 

# Plot
#TC THg Gi* 
# Importshapefile
p4_shp = readOGR(dsn = ".", layer = "p4_Outline")
# Convert shapefile to a ggplot ready data frame
p4_shp_df <- fortify(p4_shp, region = "id")
p4_shp$id <- rownames(p4_shp@data)
p4_shp_df <- left_join(p4_shp_df, p4_shp@data,mby = "id")

pelt4 = as.data.frame(pelt4)
# Plot outline and points
p4_TC = ggplot(p4_shp_df, aes(x = p4_shp_df$long, y = p4_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, color = pelt4$Gi_TC), size=5, shape=19) +
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, shape = factor(pelt4$Gi_pval_TC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 6),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

p4_UC = ggplot(p4_shp_df, aes(x = p4_shp_df$long, y = p4_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, color = pelt4$Gi_UC), size=5, shape=19) +
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, shape = factor(pelt4$Gi_pval_UC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-5, 7),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

###################################################
# Plot all figures together
library(ggpubr)
library(sjPlot)

p1=ggarrange(p1_TC, p1_UC, p2_TC, p2_UC,p3_TC, p3_UC,p4_TC, p4_UC,
          labels = c("Pelt 1 TC", "Pelt 1 UC", "Pelt 2 TC", "Pelt 2 UC", 
                     "Pelt 3 TC", "Pelt 3 UC","Pelt 4 TC", "Pelt 4 UC"),
          vjust = 1,
          hjust = -0.75,
          ncol = 4, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")

#Plot figures with dpi=300
save_plot("individual_getis.tif", p1, width = 20, height = 30, dpi = 300,
          legend.textsize = 20, legend.titlesize = 20,
          legend.itemsize = 20)

