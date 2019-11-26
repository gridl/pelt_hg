#############################################################################
# Create pelt outlines
# Written in R Version 3.5.0

# This script must be run first
#############################################################################
# Load Libraries
pelt1 = read.csv("PELT1.csv")
pelt2 = read.csv("PELT2.csv")
pelt3 = read.csv("PELT3.csv")
pelt4 = read.csv("PELT4.csv")
compiled = read.csv("pelts_compiled.csv")

# Load libraries
library(smoothr)
library(rgeos)
#############################################################################
#### Pelt 1 Outline ####
# set up coordinates
coords = coordinates(cbind(pelt1$X_coord, pelt1$Y_coord))
p1 = SpatialPoints(coords)
plot(p1, axes = TRUE)

# Buffer
buf1 = gBuffer(p1, byid=FALSE, width=5, capStyle="SQUARE")
plot(buf1)

#Smooth edges
p1_smooth = smooth(buf1, method = "chaikin")
plot(p1_smooth)

#Plot final points and polygon
plot(p1_smooth)
points(coords)

#Write spatialpolygon to dataframe and then to shapefile
df = data.frame(id = getSpPPolygonsIDSlots(p1_smooth))
row.names(df) = getSpPPolygonsIDSlots(p1_smooth)
p1_spdf <- SpatialPolygonsDataFrame(p1_smooth, data =df)
writeOGR(p1_spdf, dsn = "C:/Users/keccl081/Dropbox/pelt", driver="ESRI Shapefile", layer = "p1_Outline")

#############################################################################
##### Pelt 2 Outline ####
# set up coordinates
coords = coordinates(cbind(pelt2$X_coord, pelt2$Y_coord))
p2 = SpatialPoints(coords)
plot(p2, axes = TRUE)

# Buffer
buf1 = gBuffer(p2, byid=FALSE, width=5, capStyle="SQUARE")
plot(buf1)

#Smooth edges
p2_outline = smooth(buf1, method = "chaikin")
plot(p2_outline)
#Plot final points and polygon
p2 = plot(p2_outline)
points(coords, col = pelt2$Anatomical_region, pch=16)

#Write spatialpolygon to dataframe and then to shapefile
df = data.frame(id = getSpPPolygonsIDSlots(p2_outline))
row.names(df) = getSpPPolygonsIDSlots(p2_outline)
p2_spdf = SpatialPolygonsDataFrame(p2_outline, data =df)
writeOGR(p2_spdf, dsn = "C:/Users/keccl081/Dropbox/pelt", driver="ESRI Shapefile", layer = "p2_Outline")

#############################################################################
##### Pelt 3 Outline ####
# set up coordinates
coords = coordinates(cbind(pelt3$X_coord, pelt3$Y_coord))
p3 = SpatialPoints(coords)
plot(p3, axes = TRUE)

# Buffer
buf1 = gBuffer(p3, byid=FALSE, width=5, capStyle="SQUARE")
plot(buf1)

#Smooth edges
p3_outline = smooth(buf1, method = "chaikin")
plot(p3_outline)
#Plot final points and polygon
p3 = plot(p3_outline)
points(coords, col = pelt3$Anatomical_region, pch=16)

#Write spatialpolygon to dataframe and then to shapefile
df = data.frame(id = getSpPPolygonsIDSlots(p3_outline))
row.names(df) = getSpPPolygonsIDSlots(p3_outline)
p3_spdf = SpatialPolygonsDataFrame(p3_outline, data =df)
writeOGR(p3_spdf, dsn = "C:/Users/keccl081/Dropbox/pelt", driver="ESRI Shapefile", layer = "p3_Outline")

#############################################################################
##### Pelt 4 Outline ####
# set up coordinates
coords = coordinates(cbind(pelt4$X_coord, pelt4$Y_coord))
p4 = SpatialPoints(coords)
plot(p4, axes = TRUE)

# Buffer
buf1 = gBuffer(p4, byid=FALSE, width=5, capStyle="SQUARE")
plot(buf1)

#Smooth edges
p4_outline = smooth(buf1, method = "chaikin")
plot(p4_outline)
#Plot final points and polygon
p4 = plot(p4_outline)
points(coords, col = pelt4$Anatomical_region, pch=16)

#Write spatialpolygon to dataframe and then to shapefile
df = data.frame(id = getSpPPolygonsIDSlots(p4_outline))
row.names(df) = getSpPPolygonsIDSlots(p4_outline)
p4_spdf = SpatialPolygonsDataFrame(p4_outline, data =df)
writeOGR(p4_spdf, dsn = "C:/Users/keccl081/Dropbox/pelt", driver="ESRI Shapefile", layer = "p4_Outline")

#############################################################################
#### Pelt Compiled Outline ####
# set up coordinates
coords = coordinates(cbind(compiled$X_adj, compiled$Y_adj))
comp_p = SpatialPoints(coords)
plot(comp_p, axes = TRUE)

# Buffer
buf_comp = gBuffer(comp_p, byid=FALSE, width=.75, capStyle="SQUARE")
plot(buf_comp)

#Smooth edges
comp_smooth = smooth(buf_comp, method = "chaikin")
plot(comp_smooth)

#Plot final points and polygon
plot(comp_smooth)
points(coords)

#Write spatialpolygon to dataframe and then to shapefile
df = data.frame(id = getSpPPolygonsIDSlots(comp_smooth))
row.names(df) = getSpPPolygonsIDSlots(comp_smooth)
comp_spdf = SpatialPolygonsDataFrame(comp_smooth, data =df)
writeOGR(comp_spdf, dsn = "C:/Users/keccl081/Dropbox/pelt", driver="ESRI Shapefile", layer = "Comp_Outline_2")
