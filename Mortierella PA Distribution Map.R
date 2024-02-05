#This makes pies of species at each location

#######################################################################
#Bring in required items (These are all mapping libraries, not all used)
#######################################################################

#Load installed packages
library(tidyverse)
library(readr)
library(maps)
library(ggmap)
library(raster)
library(dismo)
library(geodata)
library(osmdata)
library(letsR)
library(ggplot2)
library(scatterpie)
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)
library(scales)


#Read in presence data, then verify everything looks okay
obs_data <- read.csv("inputfiles/Mortotpies.csv")
summary(obs_data)

#######################################################################
#Basic Mapping of Species in PA
#######################################################################

#Create PA Map with Terrain, then show in plots
pa_map_terrain <- get_map(getbb('Pennsylvania'),  maptype='terrain-background', source="stamen")
pa_map_white <- get_map(getbb('Pennsylvania'),  maptype='toner', source="stamen")

ggmap(pa_map_terrain)
ggmap(pa_map_white)

#Obtain PA County Map
counties <- map_data("county")
pa_county <- subset(counties, region == "pennsylvania")

#Map county lines on PA terrain [ggmap() + geom_polygon() fxns)], then view map
#Note: 'group' in polygon function prevents odd county lines, although doesn't fix Erie area line

pa_basemap <- ggmap(pa_map_terrain) + 
  geom_polygon(data = pa_county, aes(x=long, y=lat, group = group), fill = NA, color = "black", size=.25) 

pa_basemap

# Plot pies!
p <- pa_basemap + geom_scatterpie(aes(x=longitude, y=latitude, group = region), 
                             data = obs_data, cols = colnames(obs_data[,c(5:13)])) + 
                             labs(x= "Longitude", y= "Latitude", fill = "Species", title = ~italic("Mortierella")~ "in Pennsylvania")

p <- pa_basemap +
  geom_scatterpie(aes(x=longitude, y=latitude, group = region, r=0.005*total), data = obs_data, cols = colnames(obs_data[,c(5:13)])) +
  geom_scatterpie_legend(radius=obs_data$total*.005,x=-75.5, y=42.25, labeller=function(x) x*200) +
  labs(x= "Longitude", y= "Latitude", fill = "Species", title = ~italic("Mortierella")~ "in Pennsylvania")

p
p
  