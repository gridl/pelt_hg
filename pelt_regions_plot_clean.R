#############################################################################
# Figures illustrating different sample regions 
# Written in R Version 3.5.0
#############################################################################
# Load data
pelt1 = read.csv("PELT1.csv")
pelt2 = read.csv("PELT2.csv")
pelt3 = read.csv("PELT3.csv")
pelt4 = read.csv("PELT4.csv")

#Load Libraries
library(ggplot2)
library(ggspatial)
library(dplyr)
library(tidyr)
library(tidyverse)
library(rgeos)
library(smoothr)
#############################################################################
##### Pelt 1 Outline ####
# Import shapefile created in "pelt_outlines_clean.R"
p1_shp = readOGR(dsn = ".", layer = "p1_Outline")

# Convert shapefile to a ggplot ready data frame
p1_shp_df <- fortify(p1_shp, region = "id")
p1_shp$id <- rownames(p1_shp@data)
p1_shp_df <- left_join(p1_shp_df, p1_shp@data,mby = "id")

# Plot outline and points
p1_fur = ggplot(p1_shp_df, aes(x = p1_shp_df$long, y = p1_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$Fur_region), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

#############################################################################
##### Pelt 2 Outline ####
# Importshapefile
p2_shp = readOGR(dsn = ".", layer = "p2_Outline")

# Convert shapefile to a ggplot ready data frame
p2_shp_df <- fortify(p2_shp, region = "id")
p2_shp$id <- rownames(p2_shp@data)
p2_shp_df <- left_join(p2_shp_df, p2_shp@data,mby = "id")

# Plot outline and points
p2_fur = ggplot(p2_shp_df, aes(x = p2_shp_df$long, y = p2_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, color = pelt2$Fur_region), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

#############################################################################
##### Pelt 3 Outline ####
# Importshapefile
p3_shp = readOGR(dsn = ".", layer = "p3_Outline")

# Convert shapefile to a ggplot ready data frame
p3_shp_df <- fortify(p3_shp, region = "id")
p3_shp$id <- rownames(p3_shp@data)
p3_shp_df <- left_join(p3_shp_df, p3_shp@data,mby = "id")

# Plot outline and points
p3_fur = ggplot(p3_shp_df, aes(x = p3_shp_df$long, y = p3_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, color = pelt3$Fur_region), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

#############################################################################
##### Pelt 4 Outline ####
# Importshapefile
p4_shp = readOGR(dsn = ".", layer = "p4_Outline")

# Convert shapefile to a ggplot ready data frame
p4_shp_df <- fortify(p4_shp, region = "id")
p4_shp$id <- rownames(p4_shp@data)
p4_shp_df <- left_join(p4_shp_df, p4_shp@data,mby = "id")

# Plot outline and points
p4_fur = ggplot(p4_shp_df, aes(x = p4_shp_df$long, y = p4_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, color = pelt4$Fur_region), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

#############################################################################
#### Plot all fur region figures together ####
library(ggpubr)
library(sjPlot)

plot1=ggarrange(p1_fur, p2_fur, p3_fur, p4_fur, 
          labels = c("Pelt 1", "Pelt 2", "Pelt 3", "Pelt 4"),
          ncol = 4, nrow = 1,
          common.legend = TRUE,
          legend = "bottom")

#Plot figures with dpi=300
save_plot("fur_region.tif", plot1, width = 30, height = 20, dpi = 300,
          legend.textsize = 20, legend.titlesize = 20,
          legend.itemsize = 20)

############################################################################
##### Anatomical Regions ####
p1_ana = ggplot(p1_shp_df, aes(x = p1_shp_df$long, y = p1_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$Anatomical_region), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

p2_ana = ggplot(p2_shp_df, aes(x = p2_shp_df$long, y = p2_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, color = pelt2$Anatomical), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

p3_ana = ggplot(p3_shp_df, aes(x = p3_shp_df$long, y = p3_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, color = pelt3$Anatomical), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

p4_ana = ggplot(p4_shp_df, aes(x = p4_shp_df$long, y = p4_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, color = pelt4$Anatomical_region), size=4, shape=19)+
  scale_colour_manual(values = c("#225ea8", "#FDD835", "#41b6c4", "#a1dab4"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

#############################################################################
#### Plot all figures together ####
plot2=ggarrange(p1_ana, p2_ana, p3_ana, p4_ana, 
                labels = c("Pelt 1", "Pelt 2", "Pelt 3", "Pelt 4"),
                ncol = 4, nrow = 1,
                common.legend = TRUE,
                legend = "bottom")

#Plot figures with dpi=300
save_plot("anatomical_region.tif", plot2, width = 30, height = 20, dpi = 300)
