################################################
##### Prepare target group background data #####
################################################

# load the libraries
library(rgbif)
library(geodata)

library(Taxonstand)
library(CoordinateCleaner)
library(maps)

library(spThin)
library(dplyr)
library(readr)

# Set your paths
source("paths.R")

# get the data
species<- "Bulinus tropicus"
occs<-sp_occurrence(genus="Bulinus",species="tropicus",removeZeros=TRUE,download=TRUE,geo=T)
nrow(occs)
colnames(occs)
head(occs)

#write the raw data to a csv file
write.csv(occs,"./raw_data/raw_data_Bt.csv", row.names=FALSE)

# Check the species' taxonomy
sort(unique(occs$scientificName))

#show the currently accepted taxonomy
table(occs$taxonomicStatus)

#Check the species' coordinates
plot(lat~lon,data=occs)
map(,,,add=TRUE)

#clean the coordinates:
occs.coord <- occs[!is.na(occs$lat)& !is.na(occs$lon),] # remove observations with NA values 
geo.clean <- clean_coordinates(x=occs.coord,lon="lon",lat="lat", species="species",value="clean")

# plot the cleaned data:
par(mfrow=c(1,2))
plot(lat~lon, data=occs)
map(,,,add=TRUE)
plot(lat~lon,data=geo.clean)
map(,,,add=TRUE)
par(mfrow=c(1,1))

head(geo.clean)

# now thin the data in the same way as we did with the presence data points
thinned_data <- thin(loc.data=geo.clean, long.col = "lon", lat.col="lat", spec.col="species", 
                     thin.par=5, reps=10, # thin.par specifies the minimum distance between points in kilometers
                     locs.thinned.list.return = TRUE,
                     write.files=FALSE,write.log.file=FALSE)

#some datasets retain more occurences than others. find one with nrow that has the most points retained and plot this layer.
nrow(thinned_data[[1]])
plot(lat~lon,data=geo.clean)
map(,,,add=TRUE)
points(thinned_data[[1]]$Longitude, thinned_data[[1]]$Latitude, col="red")
thinned<-thinned_data[[1]]

names(thinned)[names(thinned)=="Longitude"]<-"x" 
names(thinned)[names(thinned)=="Latitude"]<-"y"
head(thinned)

#clip the thinned layer to the area of interest if necessary:
shp<-readOGR(path2clipshp)

WGScoor<-thinned
coordinates(WGScoor)<-~Longitude+Latitude
proj4string(WGScoor)<-CRS("+proj=longlat +datum=WGS84")
plot(shp)
points(WGScoor)

clip.r<-crop(WGScoor,shp,snap="in")
plot(shp)
points(clip.r)
clip.r
thinned<-as.data.frame(clip.r)

#write this thinned data to a csv file
write.csv(thinned,"./cleaned_data_target_group.csv",row.names=FALSE)

##############################################
###### Creating biased background layer ######
##############################################


#load libraries
library(raster)
library(dplyr)
library(ks)
library(spatialEco)
library(terra)
library(dismo)

#set up paths: 
setwd("/path_2_working_directory")
#path to environmental variables
path2layers<-"./path_2_environmental_training_layers"
#path to occurence data target group species
path2sp<-"./target_group_presences.shp"
#path to output bias file:
path2output<-"./TargetGroup_biasfile.tif"
#path to presence points
path2presence<-"./thinned_data.csv"
#path to output file for background points
path2bckgroutput<-"./biased_background_points.csv"


### read in data ###


#read in environmental rasters
layers<-stack(list.files(path=path2layers,pattern='asc',full.names=TRUE ))
bio_stack<-layers[[c("bio4","max_NDVI","bio17","bio2","bio5","bio6","bio1","bio12","river_slope")]]

# read in occurence data:
all_snails<-shapefile(path2sp)
all_presences <- as.data.frame(all_snails)
head(all_presences)
all_presences<-all_presences[,c(14,15)]
names(all_presences)<-c("x","y")

#set seed.
set.seed(12345)


### Create target group bias file ###

# aggregate target group presences by coordinates
sum_records<-as.data.frame(dplyr::count(all_presences,x,y))

#extract coordinates
coords<-cbind(sum_records[,1],sum_records[,2])

#create scale
scale<-length(sum_records[,3])/sum(sum_records[,3])
sum_records[,4]<-sum_records[,3]*scale

#do a 2d kernel density estimation.
target_density<-kde(coords,w=sum_records[,4])

#create raster
target_raster<-raster(target_density)
plot(target_raster)

# assign a coordinate reference system
crs(target_raster)<-'+init=EPSG:4326'

# Clip data to the same resolution and extent as the environmental layers:
target_raster<-resample(target_raster,bio_stack,method="bilinear")

#mask bias file
for (x in 1:nlayers(bio_stack)){
  target_raster<-mask(target_raster,bio_stack[[x]])
}

#normalize bias file between 0 and 1
target_raster<-target_raster-minValue(target_raster)
target_raster<-rast(target_raster) # convert to a spatial raster
target_raster<-raster.transformation(target_raster,trans="norm")
plot(target_raster)

#export the raster
writeRaster(target_raster,path2output,overwrite=T)

target_biasfile<-raster(path2output)

#select background points:
presences<-read.csv(path2presence)
presences<-presences[,c(1:2)]
head(presences)
background_points<-randomPoints(mask=target_biasfile,n=8000,p=presences,,excludep=TRUE,prob=TRUE)
plot(target_raster)
points(background_points)
write.csv(background_points,file=path2bckgroutput)

#select additional background points from the sahara and the middle east: 
path2clipshp<-"./clipping_area_only_Sahara.shp"
shp<-shapefile(path2clipshp)
plot(shp)
#now clip the raster layer to the clipping area extent
clip.r<-crop(bio_stack$bio5,shp,snap="in")
# mask every cel outside the clipping area
mask.r<-mask(clip.r,shp)
plot(mask.r)
bckgr<-as.data.frame(background_points)
exclude_points<-rbind(bckgr,presences)
background_points_SAH<-randomPoints(mask=mask.r,n=2000,p=exclude_points,,excludep=TRUE,prob=FALSE)
points(background_points_SAH)
all_background<-rbind(background_points,background_points_SAH)
write.csv(all_background,file=path2bckgroutput)