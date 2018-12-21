#############################################################################
# Plots of individual Hg concentrations 
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
library(dplyr)
library(ggplot2)
library(ggpubr)
library(sjPlot)

#############################################################################
# Plot 1
# Outline
# Importshapefile
p1_shp = readOGR(dsn = ".", layer = "p1_Outline")
# Convert shapefile to a ggplot ready data frame
p1_shp_df <- fortify(p1_shp, region = "id")
p1_shp$id <- rownames(p1_shp@data)
p1_shp_df <- left_join(p1_shp_df, p1_shp@data,mby = "id")

# Plot
p1_TC = ggplot(p1_shp_df, aes(x = p1_shp_df$long, y = p1_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$THg_TC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0',
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  theme(legend.key.size=10)+
  coord_quickmap()+
  theme_void()

p1_UC = ggplot(p1_shp_df, aes(x = p1_shp_df$long, y = p1_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$THg_UC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0',
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

#############################################################################
# Plot 2
# Outline
# Importshapefile
p2_shp = readOGR(dsn = ".", layer = "p2_Outline")
# Convert shapefile to a ggplot ready data frame
p2_shp_df <- fortify(p2_shp, region = "id")
p2_shp$id <- rownames(p2_shp@data)
p2_shp_df <- left_join(p2_shp_df, p2_shp@data,mby = "id")

# Plot
p2_TC = ggplot(p2_shp_df, aes(x = p2_shp_df$long, y = p2_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, color = pelt2$THg_TC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0', 
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

p2_UC = ggplot(p2_shp_df, aes(x = p2_shp_df$long, y = p2_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, color = pelt2$THg_UC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0', 
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

#############################################################################
# Plot 3
# Outline
# Importshapefile
p3_shp = readOGR(dsn = ".", layer = "p3_Outline")
# Convert shapefile to a ggplot ready data frame
p3_shp_df <- fortify(p3_shp, region = "id")
p3_shp$id <- rownames(p3_shp@data)
p3_shp_df <- left_join(p3_shp_df, p3_shp@data,mby = "id")

# Plot
p3_TC = ggplot(p3_shp_df, aes(x = p3_shp_df$long, y = p3_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, color = pelt3$THg_TC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0', 
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

p3_UC = ggplot(p3_shp_df, aes(x = p3_shp_df$long, y = p3_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, color = pelt3$THg_UC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0', 
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

#############################################################################
# Plot 4
# Outline
# Importshapefile
p4_shp = readOGR(dsn = ".", layer = "p4_Outline")
# Convert shapefile to a ggplot ready data frame
p4_shp_df <- fortify(p4_shp, region = "id")
p4_shp$id <- rownames(p4_shp@data)
p4_shp_df <- left_join(p4_shp_df, p4_shp@data,mby = "id")

# Plot
p4_TC = ggplot(p4_shp_df, aes(x = p4_shp_df$long, y = p4_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, color = pelt4$THg_TC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0',
                       guide = guide_colourbar(title = "Total Hg (ug/g)"))+
  coord_quickmap()+
  theme_void()

p4_UC = ggplot(p4_shp_df, aes(x = p4_shp_df$long, y = p4_shp_df$lat)) +
  geom_polygon(colour='black', fill='white') +
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, color = pelt4$THg_UC), size=4, shape=19) +
  scale_color_gradient(low='#ece7f2', high='#0570b0',
                       guide = guide_colourbar(title = "Total Hg (ug/g)")) +
  coord_quickmap()+
  theme_void()

#############################################################################
# Plot all figures together
plot1=ggarrange(p1_TC, p1_UC, p2_TC, p2_UC,p3_TC, p3_UC,p4_TC, p4_UC,
          labels = c("Pelt 1 TC", "Pelt 1 UC", "Pelt 2 TC", "Pelt 2 UC", "Pelt 3 TC", "Pelt 3 UC","Pelt 4 TC", "Pelt 4 UC"),
          vjust = 1,
          hjust = -0.5,
          ncol = 4, nrow = 2,
          common.legend = FALSE,
          legend = "right")

#Plot figures with dpi=300
save_plot("individual_hg.tif", plot1, width = 20, height = 20, dpi = 300)
