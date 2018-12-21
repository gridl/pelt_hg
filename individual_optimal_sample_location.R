#############################################################################
# Otter fur THg to organ THg otter sampling spot optimization
# Written in R Version 3.5.0
#############################################################################
# Load data
pelt1 = read.csv("PELT1.csv")
pelt2 = read.csv("PELT2.csv")
pelt3 = read.csv("PELT3.csv")
pelt4 = read.csv("PELT4.csv")


# Library
library(ggplot2)
library(Metrics)
library(rgdal)
library(dplyr)

#ADD COLUMNS PREDICTING [ORGAN THg] from [Fur THg]####
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

# Pelt1
p1_shp = readOGR(dsn = ".", layer = "p1_Outline")
# Convert shapefile to a ggplot ready data frame
p1_shp_df <- fortify(p1_shp, region = "id")
p1_shp$id <- rownames(p1_shp@data)
p1_shp_df <- left_join(p1_shp_df, p1_shp@data,mby = "id")

# Pelt2
p2_shp = readOGR(dsn = ".", layer = "p2_Outline")
# Convert shapefile to a ggplot ready data frame
p2_shp_df <- fortify(p2_shp, region = "id")
p2_shp$id <- rownames(p2_shp@data)
p2_shp_df <- left_join(p2_shp_df, p2_shp@data,mby = "id")

# Pelt3
p3_shp = readOGR(dsn = ".", layer = "p3_Outline")
# Convert shapefile to a ggplot ready data frame
p3_shp_df <- fortify(p3_shp, region = "id")
p3_shp$id <- rownames(p3_shp@data)
p3_shp_df <- left_join(p3_shp_df, p3_shp@data,mby = "id")

# Pelt4
p4_shp = readOGR(dsn = ".", layer = "p4_Outline")
# Convert shapefile to a ggplot ready data frame
p4_shp_df <- fortify(p4_shp, region = "id")
p4_shp$id <- rownames(p4_shp@data)
p4_shp_df <- left_join(p4_shp_df, p4_shp@data,mby = "id")

#######################################
# pelt 1
# Brain
pelt1$est_brain= 0.15*pelt1$THg_TC
pelt1$resid_brain=pelt1$est_brain-pelt1$Brain
pelt1$perc_brain= abs(((pelt1$est_brain-pelt1$Brain)/(pelt1$Brain))*100)

# Plot outline and points
p1_brain = ggplot(p1_shp_df, aes(x = p1_shp_df$long, y = p1_shp_df$lat))+
  geom_polygon(colour='black', fill='white')+
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$perc_brain), size=4, shape=19)+
  #scale_colour_fill(values = c("red", "white"))+
  labs(colour=" ")+
  theme(legend.title = element_blank())+
  coord_quickmap()+
  theme_void()

# Liver
pelt1$est_liver= 0.70*pelt1$THg_TC
pelt1$resid_liver=pelt1$est_liver-pelt1$Liver
pelt1$perc_liver= abs(((pelt1$est_liver-pelt1$Liver)/(pelt1$Liver))*100)
ggplot() + geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, colour =perc_liver), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Kidney
pelt1$est_kidney= 0.46*pelt1$THg_TC
pelt1$resid_kidney=pelt1$est_kidney-pelt1$Kidney
pelt1$perc_kidney= abs(((pelt1$est_kidney-pelt1$Kidney)/(pelt1$Kidney))*100)
ggplot() + geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, colour =perc_kidney), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Muscle
pelt1$est_muscle = 0.62*pelt1$THg_TC
pelt1$resid_muscle = pelt1$est_muscle-pelt1$Muscle
pelt1$perc_muscle = abs(((pelt1$est_muscle-pelt1$Muscle)/(pelt1$Muscle))*100)
ggplot() + geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, colour =perc_muscle), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Average error
pelt1$p1_average_error = (pelt1$perc_brain + pelt1$perc_liver + pelt1$perc_kidney + pelt1$perc_muscle)/4

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
  geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, color = pelt1$p1_average_error), size=5, shape=19) +
  scale_color_gradient2(high = "red", low = "yellow", limits = c(0, 100),
                        guide = guide_colourbar(title = "Error")) +
  coord_quickmap()+
  theme_void()

################################################
# Pelt2
# Brain
pelt2$est_brain= 0.15*pelt2$THg_TC
pelt2$resid_brain=pelt2$est_brain-pelt2$Brain
pelt2$perc_brain= abs(((pelt2$est_brain-pelt2$Brain)/(pelt2$Brain))*100)
ggplot() + geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, colour =perc_brain), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Liver
pelt2$est_liver= 0.70*pelt2$THg_TC
pelt2$resid_liver=pelt2$est_liver-pelt2$Liver
pelt2$perc_liver= abs(((pelt2$est_liver-pelt2$Liver)/(pelt2$Liver))*100)
ggplot() + geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, colour =perc_liver), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Kidney
pelt2$est_kidney= 0.46*pelt2$THg_TC
pelt2$resid_kidney=pelt2$est_kidney-pelt2$Kidney
pelt2$perc_kidney= abs(((pelt2$est_kidney-pelt2$Kidney)/(pelt2$Kidney))*100)
ggplot() + geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, colour =perc_kidney), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Muscle
pelt2$est_muscle = 0.62*pelt2$THg_TC
pelt2$resid_muscle = pelt2$est_muscle-pelt2$Muscle
pelt2$perc_muscle = abs(((pelt2$est_muscle-pelt2$Muscle)/(pelt2$Muscle))*100)
ggplot() + geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, colour =perc_muscle), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Average error
pelt2$p2_average_error = (pelt2$perc_brain + pelt2$perc_liver + pelt2$perc_kidney + pelt2$perc_muscle)/4
ggplot() + geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, colour =p2_average_error), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

###################################################################################
# Pelt 3
# Brain
pelt3$est_brain= 0.15*pelt3$THg_TC
pelt3$resid_brain=pelt3$est_brain-pelt3$Brain
pelt3$perc_brain= abs(((pelt3$est_brain-pelt3$Brain)/(pelt3$Brain))*100)

# Liver
pelt3$est_liver= 0.70*pelt3$THg_TC
pelt3$resid_liver=pelt3$est_liver-pelt3$Liver
pelt3$perc_liver= abs(((pelt3$est_liver-pelt3$Liver)/(pelt3$Liver))*100)
ggplot() + geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, colour =perc_liver), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Kidney
pelt3$est_kidney= 0.46*pelt3$THg_TC
pelt3$resid_kidney=pelt3$est_kidney-pelt3$Kidney
pelt3$perc_kidney= abs(((pelt3$est_kidney-pelt3$Kidney)/(pelt3$Kidney))*100)
ggplot() + geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, colour =perc_kidney), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Muscle
pelt3$est_muscle = 0.62*pelt3$THg_TC
pelt3$resid_muscle = pelt3$est_muscle-pelt3$Muscle
pelt3$perc_muscle = abs(((pelt3$est_muscle-pelt3$Muscle)/(pelt3$Muscle))*100)
ggplot() + geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, colour =perc_muscle), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Average error
pelt3$p3_average_error = (pelt3$perc_brain + pelt3$perc_liver + pelt3$perc_kidney + pelt3$perc_muscle)/4
ggplot() + geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, colour =p3_average_error), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

############################################################
# Pelt 4
# Brain
pelt4$est_brain= 0.15*pelt4$THg_TC
pelt4$resid_brain=pelt4$est_brain-pelt4$Brain
pelt4$perc_brain= abs(((pelt4$est_brain-pelt4$Brain)/(pelt4$Brain))*100)
ggplot() + geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, colour =perc_brain), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Liver
pelt4$est_liver= 0.70*pelt4$THg_TC
pelt4$resid_liver=pelt4$est_liver-pelt4$Liver
pelt4$perc_liver= abs(((pelt4$est_liver-pelt4$Liver)/(pelt4$Liver))*100)
ggplot() + geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, colour =perc_liver), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Kidney
pelt4$est_kidney= 0.46*pelt4$THg_TC
pelt4$resid_kidney=pelt4$est_kidney-pelt4$Kidney
pelt4$perc_kidney= abs(((pelt4$est_kidney-pelt4$Kidney)/(pelt4$Kidney))*100)
ggplot() + geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, colour =perc_kidney), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Muscle
pelt4$est_muscle = 0.62*pelt4$THg_TC
pelt4$resid_muscle = pelt4$est_muscle-pelt4$Muscle
pelt4$perc_muscle = abs(((pelt4$est_muscle-pelt4$Muscle)/(pelt4$Muscle))*100)
ggplot() + geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, colour =perc_muscle), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Average error
pelt4$p4_average_error = (pelt4$perc_brain + pelt4$perc_liver + pelt4$perc_kidney + pelt4$perc_muscle)/4
ggplot() + geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, colour =p4_average_error), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

#########################################################################
# Compiled
# Brain
compiled$est_brain= 0.15*compiled$THg_TC
compiled$resid_brain=compiled$est_brain-compiled$Brain
compiled$perc_brain= abs(((compiled$est_brain-compiled$Brain)/(compiled$Brain))*100)
ggplot() + geom_point(data = compiled, aes(x = X_coord, y = Y_coord, colour =perc_brain), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Liver
compiled$est_liver= 0.70*compiled$THg_TC
compiled$resid_liver=compiled$est_liver-compiled$Liver
compiled$perc_liver= abs(((compiled$est_liver-compiled$Liver)/(compiled$Liver))*100)
ggplot() + geom_point(data = compiled, aes(x = X_coord, y = Y_coord, colour =perc_liver), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Kidney
compiled$est_kidney= 0.46*compiled$THg_TC
compiled$resid_kidney=compiled$est_kidney-compiled$Kidney
compiled$perc_kidney= abs(((compiled$est_kidney-compiled$Kidney)/(compiled$Kidney))*100)
ggplot() + geom_point(data = compiled, aes(x = X_coord, y = Y_coord, colour =perc_kidney), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Muscle
compiled$est_muscle = 0.62*compiled$THg_TC
compiled$resid_muscle = compiled$est_muscle-compiled$Muscle
compiled$perc_muscle = abs(((compiled$est_muscle-compiled$Muscle)/(compiled$Muscle))*100)
ggplot() + geom_point(data = compiled, aes(x = X_coord, y = Y_coord, colour =perc_muscle), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

# Average error
compiled$comp_average_error = (compiled$perc_brain + compiled$perc_liver + compiled$perc_kidney + compiled$perc_muscle)/4
ggplot() + geom_point(data = compiled, aes(x = X_coord, y = Y_coord, colour =comp_average_error), size=5)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

########################################################################
Plots

p1 = ggplot() + geom_point(data = pelt1, aes(x = X_coord, y = Y_coord, colour =p1_average_error), size=4)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

p2=  ggplot() + geom_point(data = pelt2, aes(x = X_coord, y = Y_coord, colour =p2_average_error), size=4)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

p3 = ggplot() + geom_point(data = pelt3, aes(x = X_coord, y = Y_coord, colour =p3_average_error), size=4)+
  scale_colour_gradient(low = "white", high = "red")+
  coord_quickmap()+
  theme_void()

p4 = ggplot() + geom_point(data = pelt4, aes(x = X_coord, y = Y_coord, colour =p4_average_error), size=4)+
  scale_colour_gradient(low = "white", high = "red")+
  guides(fill=guide_legend(title="Average Error"))+
  coord_quickmap()+
  theme_void()

library(gridExtra)
grid.arrange(grobs = list(p1, p2, p3, p4), labels = c("A", "B", "C"), ncol = 2)

# Composite error image
#RMSE
rmse_brain= (rmse(pelt1$Brain,pelt1$est_brain)) 
rmse_liver= rmse(pelt1$Liver,pelt1$est_liver) 
rmse_kidney= rmse(pelt1$Kidney,pelt1$est_kidney)
rmse_muscle= rmse(pelt1$Muscle,pelt1$est_muscle)


