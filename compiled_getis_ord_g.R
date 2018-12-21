#############################################################################
# Cluster analysis for composite pelt
# Written in R Version 3.5.0
#############################################################################
# Load data
compiled = read.csv("pelts_compiled.csv")

library(ggplot2)
library(ggspatial)
library(dplyr)
library(tidyr)
library(rgeos)
library(smoothr)
library(sp)
library(rgdal)

################################
# Importshapefile
comp_shp = readOGR(dsn = ".", layer = "Comp_Outline_2")
# Convert shapefile to a ggplot ready data frame
comp_shp_df <- fortify(comp_shp, region = "id")
comp_shp$id <- rownames(comp_shp@data)
comp_shp_df <- left_join(comp_shp_df, comp_shp@data,mby = "id")

#Layered lat long
ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = compiled, aes(x = X_adj, y = Y_adj, color = factor(compiled$Pelt)), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

# Plot Normalized Hg Values for TC and UC
comp_TC = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = compiled, aes(x = X_adj, y = Y_adj, color = compiled$TC_norm), size=7, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0',
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

comp_UC = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = compiled, aes(x = X_adj, y = Y_adj, color = compiled$UC_norm), size=7, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0',
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

############################################################
#### Average Noralized values ####
adj_x_compiled = aggregate(compiled$X_adj, list(compiled$factor), mean)
adj_y_compiled = aggregate(compiled$Y_adj, list(compiled$factor), mean)
TC_compiled = aggregate(compiled$TC_norm, list(compiled$factor), mean)
UC_compiled = aggregate(compiled$UC_norm, list(compiled$factor), mean)
compiled_pelt = as.data.frame(cbind(adj_x_compiled[,2], adj_y_compiled[,2],TC_compiled[,2],UC_compiled[,2]))
colnames(compiled_pelt) <- c("x", "y", "TC", "UC")

# Plot
ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = compiled_pelt, aes(x = x, y = y, color = TC), size=4)+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()


##############################################
### Compiled Pelt Hot Spot Analysis ####
#Set up spatial neighbourhood
comp_sppnt = SpatialPoints(cbind(compiled_pelt$x, compiled_pelt$y))
comp_nb <- chooseCN(coordinates(comp_sppnt), type = 6, k = 12, plot.nb = FALSE) 
comp_lw <- nb2listw(comp_nb, style = 'W', zero.policy = TRUE)

# Gi*
compiled_pelt$Gi_TC = as.numeric(localG(compiled_pelt$TC, listw = comp_lw))
compiled_pelt$Gi_UC = as.numeric(localG(compiled_pelt$UC, listw = comp_lw))

# Gi* p-value
#Z score values for levels of statistical significance:
# 90% significant: >= 1.645
# 95% significant: >= 1.960
# 99% significant: >= 2.576
# 99.9% significant: >= 3.291

# Bin gi* scores based on whether it is past the critical cut off
compiled_pelt$Gi_pval_TC = as.factor(compiled_pelt$Gi_TC >= 1.96 |  compiled_pelt$Gi_TC <= -1.96) 
compiled_pelt$Gi_pval_UC = as.factor(compiled_pelt$Gi_UC >= 1.96 |  compiled_pelt$Gi_UC <= -1.96) 

# Plot 1
#TC THg Gi*
# Plot outline and points
hs_comp_TC = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = compiled_pelt, aes(x = x, y = y, color = compiled_pelt$Gi_TC), size=5, shape=19) +
  geom_point(data = compiled_pelt, aes(x = x, y = y, shape = factor(compiled_pelt$Gi_pval_TC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi* score significance"), 
                     labels = c("Not Significant", "Significant")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-6, 6),
                        guide = guide_colourbar(title = "Gi*score")) +
  coord_quickmap()+
  theme_void()

hs_comp_UC = ggplot(comp_shp_df, aes(x = comp_shp_df$long, y = comp_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = compiled_pelt, aes(x = x, y = y, color = compiled_pelt$Gi_UC), size=5, shape=19) +
  geom_point(data = compiled_pelt, aes(x = x, y = y, shape = factor(compiled_pelt$Gi_pval_UC)), size = 6) +
  scale_shape_manual(values = c(1,13),
                     guide = guide_legend(title = "Gi*score significance")) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue", midpoint = 0, limits = c(-6, 6),
                        guide = guide_colourbar(title = "Gi* score")) +
  coord_quickmap()+
  theme_void()

###################################################
# Plot composite figures together
library(ggpubr)
plot1=ggarrange(comp_TC, comp_UC,
          labels = c("Compiled TC", "Compiled UC"),
          vjust = 1,
          hjust = -0.5,
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = "right")

#Plot figures with dpi=300
save_plot("compsite_pelt.tif", plot1, width = 15, height = 15, dpi = 300)

#Plot composite figure hotspots together
plot2=ggarrange(hs_comp_TC, hs_comp_UC,
          labels = c("Hot Spot TC", "Hot Spot UC"),
          vjust = 1,
          hjust = -0.5,
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = "right")

#Plot figures with dpi=300
save_plot("compsite_pelt_hs.tif", plot2, width = 15, height = 15, dpi = 300)


