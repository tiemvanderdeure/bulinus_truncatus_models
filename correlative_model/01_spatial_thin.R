##############################
###### Spatial thinning ######
##############################


# Load libraries
library(spThin)
library(ggplot2)
library(rnaturalearth)

# Set your paths
source("paths.R")

#read in the data
data <- read.csv(path2allpresences, header=T)

#plot the species presence points
plot(data$x, data$y)
data$Species<-"target_group"
#thin the data
thinned_data_full <- thin(loc.data=data, long.col = "x", lat.col="y", spec.col="Species", 
                          thin.par=5, reps=100, # thin.par specifies the minimum distance between points in kilometers
                          locs.thinned.list.return = TRUE,
                          write.files=FALSE,write.log.file=FALSE)

#All datasets seem to retain the same number of occurences. Plot one of them above the previous plot to show which points are retained.
thinned<-thinned_data_full[[1]]
points(thinned$Longitude, thinned$Latitude, col="red")# plots the retained points over the plot previously made 

names(thinned) <- c("x", "y")

# now retrieve all the columns from the data dataframe and merge them with the columns from the thinned dataframe
spt<-merge(thinned,data,by= c("row.names", "x", "y"))

#Plot the datapoints on a map:
sf_world <- ne_countries(returnclass = "sf")
tibble::glimpse(sf_world)
ggplot(sf_world) +
  geom_sf()+
  geom_point(data=thinned, aes(x=x, y=y))

#export the dataframe with the selected points
write.csv(thinned,file=path2presence,row.names=FALSE)

#Transform this data into a shapefile for easy visualisation in gis software
# WGScoor <- spt
# coordinates(WGScoor)=~x+y
# proj4string(WGScoor)<-CRS("+proj=longlat +datum=WGS84")
# raster::shapefile(WGScoor,"thinned_data.shp",overwrite=TRUE)