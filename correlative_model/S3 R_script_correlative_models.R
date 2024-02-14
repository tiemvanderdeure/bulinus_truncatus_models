#########################################
############### BIOMOD2 #################
#########################################


##############################
###### Spatial thinning ######
##############################


# Load libraries
library(spThin)
library(raster)
library(maptools)
library(rgdal)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)

# Set your paths

#set working directory
setwd("/path_2_species_presence_datafile")
path2data<-"/path_2_species_presence_datafile/species_presence_datafile.csv"
#read in the data
data=read.csv(path2data, header=T, sep=",",dec=".")
head(data)

#plot the species presence points
plot(data$x.1, data$y.1)
data$Species<-"target_group"
#thin the data
thinned_data_full <- thin(loc.data=data, long.col = "x.1", lat.col="y.1", spec.col="Species", 
                          thin.par=5, reps=100, # thin.par specifies the minimum distance between points in kilometers
                          locs.thinned.list.return = TRUE,
                          write.files=FALSE,write.log.file=FALSE)

#All datasets seem to retain the same number of occurences. Plot one of them above the previous plot to show which points are retained.
points(thinned_data_full[[1]]$x.1, thinned_data_full[[1]]$y.1, col="red")# plots the retained points over the plot previously made 
nrow(thinned_data_full[[1]])
View(thinned_data_full[[1]])
thinned<-thinned_data_full[[1]]

names(thinned)[names(thinned)=="Longitude"]<-"x" 
names(thinned)[names(thinned)=="Latitude"]<-"y"
head(thinned)
nrow(thinned)

# now retrieve all the columns from the data dataframe and merge them with the columns from the thinned dataframe
spt<-merge(thinned,data,by="row.names")

#clean up this new dataframe:
head(spt)
spt<-spt[,!names(spt) %in% c("Row.names","y.y","x.y")]
spt<-rename(spt,x=x.x,y=y.x)
head(spt)

#Plot the datapoints on a map:

sf_world <- ne_countries(returnclass = "sf")
tibble::glimpse(sf_world)
ggplot(sf_world) +
  geom_sf()+
  geom_point(data=thinned, aes(x=x, y=y))

#export the dataframe with the selected points
write.csv(thinned,file=paste("thinned_data.csv",sep=''),row.names=FALSE)

#Transform this data into a shapefile for easy visualisation in gis software
WGScoor <- spt
coordinates(WGScoor)=~x+y
proj4string(WGScoor)<-CRS("+proj=longlat +datum=WGS84")
raster::shapefile(WGScoor,"thinned_data.shp",overwrite=TRUE)



################################################
##### Prepare target group background data #####
################################################

# load the libraries
library(rgbif)
library(Taxonstand)
library(CoordinateCleaner)
library(maps)
library(geodata)
library(spThin)
library(dplyr)
library(readr)

# Set your paths
setwd("path_to_Working_directory")

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
path2clipshp<-"/clipping_area_training_layers.shp"
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

###################################
# Prepare environmental variables #
###################################



########################## CLIPPING RASTERS #############################

setwd("/path_2_environmental layers")

library(maps)
library(raster)
library(maptools)
library(rgdal)
library(terra)

path2layers <- "/path_2_environmental_layers" # path to folder with all environmental variables of interest
path2clipshp<-"/path_2_clip_shapefile/clipping_area.shp"
path2cliplayers <- "/path_2_clipped_layers"
path2resamplayers<-"/path_2_resample_template/env_layer.asc"

# Load environmental variables
layers<-list.files(path=path2layers,pattern='asc',full.names=TRUE)
# Join all rasters in a single object
egv<-stack(layers)
# check the structure of the rasters
print(egv)


# Load clipping area shapefile
shp<-readOGR(path2clipshp)
# Plot the shapefile
plot(shp)
map(,,,add=T)

# Loop to clip each raster sequentially
for (i in layers) { # look for each ascii file in layers
  #load raster
  egv.i<-raster(i)
  # clip raster
  # the raster to be clipped; the clipping feature; snap determines in which direction the extent should be aligned
  # to the nearest border, inwards or outwards
  clip.r<-crop(egv.i, shp, snap="in")
  # plot the clipped raster
  # X11()
  #  plot(clip.r,main=clip.r@data@names) #main= title
  #mask raster with the shapefile to create NoData 
  mask.r<-mask(clip.r, shp)
  # plot the masked raster
  # X11()
  #  plot(mask.r,main=mask.r@data@names)
  # Get name of the band
  name<-mask.r@data@names
  # resample the layer if necessary:
  template=raster(path2resamplayers)
  resampled.mask<-resample(mask.r,template,method="bilinear")
  # export the raster to a folder using paste
  writeRaster(resampled.mask,paste(path2cliplayers,"clipped_",name,sep=''),format="ascii", overwrite=TRUE)
}

##################### ASSESSING CORRELATION AMONG VARIABLES ########################

# Load libraries
library(raster)
library(maptools)
library(rgdal)

# set paths
path2cliplayers <- "/path_2_clipped_layers"
# path to correlation output folder
path2corr<-"./04_correlation_results/"


# Load environmental variables
layers<-list.files(path=path2cliplayers,pattern='asc',full.names=TRUE)
layers
# Join all rasters in a single object
egv<-stack(layers)
# check the structure of the rasters
print(egv)
names(egv)
# check descriptive statistics of the rasters
summary(egv)
# Plot the stack
plot(egv)

#IF APPLICABLE: make a subset of the raster stack with the variables that you want to carry out a correlation analysis on
egv.subset<-egv[[c("max_NDVI","bio5","bio6","bio4","bio12","soil_water_infiltration_flux","sm_surface","river_slope","elevation","land_use_2016","land_fraction_saturated","vegetation_greenness_fraction","bio17","bio2","slope","soil_calcium","bio1")]]
names(egv.subset)

# Calculate Pearson's correlation among variables
# Pearson's correlation is parametric
pearson<-layerStats(egv.subset,'pearson',na.rm=T)
# Check correlation results
print(pearson)

# Transform correlation results in a dataframe
cor.df<-as.data.frame(pearson)
# Check dataframe structure
head(cor.df)
cor.df<-cor.df[ , -which(names(cor.df) %in% c("mean"))]
#cor.df<-sub("*.clip_","",cor.df[,1])
# Export dataframe
setwd(path2corr)
write.table(cor.df,paste('correlations.csv',sep=''),sep=",") 
# Transform correlation matrix to distances
var.dist <- abs(as.dist(cor.df))
# Calculate dendrogram based on distance (less distance = more correlation)
var.cluster <- hclust(1-var.dist)
# Plot dendrogram
plot(var.cluster)
abline(h=0.25, lty=2, lwd=2, col="red") # variables that have a correlation < 0.75
# Export correlation graph
jpeg(paste('correlation_graph.jpg',sep='')) #better to use export function in R studio
plot(var.cluster)
abline(h=0.25, lty=2, lwd=2) # variables that have a correlation < 0.75
dev.off()

########################## Plot the pca #############################

# Load libraries
library(ade4)
library(raster)
library(gridExtra)

# path to clipped layers
path2cliplayers <- "/path_2_clipped_layers"

#path to presence points
path2presence<-"./thinned_data.csv"
#path background points
path2bckgroutput<-"./biased_background_points.csv"

#stack the environmental layers
layers<-stack(list.files(path=path2cliplayers,pattern='asc',full.names=TRUE ))
names(layers)
myExpl<-layers[[c("bio4","bio12","max_NDVI","bio17","bio2","bio5","sm_surface","soil_water_infiltration_flux","bio6","bio1","elevation","river_slope")]]

# load our species data
FullDataSpecies <- read.csv(path2sp,sep=',')
# Check species data
head(FullDataSpecies)
str(FullDataSpecies)
nrow(FullDataSpecies)
#keep only the columns with the species name, x and y coordinates and place the column with the presence or absence of the species first
DataSpecies <- FullDataSpecies[,c(2:3)]
names(DataSpecies)<-c("decimalLongitude","decimalLatitude")
head(DataSpecies)

#obtain cell identifiers of the cells where B. truncatus occurs
cell_id<-cellFromXY(subset(myExpl,1),DataSpecies)

#remove non-defined area from the environmental layers:
layers_df<-na.omit(as.data.frame(myExpl))
head(layers_df)

#calculate a pca from the environmental layers:
pca.env<-dudi.pca(layers_df,scannf=F,nf=2)

#investigate the distribution of the target species in environmental space
par(mfrow=c(1,2))
s.class(pca.env$li[,1:2],fac=factor(rownames(layers_df) %in% cell_id, levels=c("FALSE","TRUE"),
                                    labels=c("background","B.truncatus")), col=c("red","blue"),
        csta=0, cellipse=2,cpoint=.3, pch=16)
mtext("(a)",side=3,line=3,adj=0)
s.corcircle(pca.env$co,clabel=.5)
mtext("(b)",side=3,line=3,adj=0)

########### DENSITY PLOTS BACKGROUND/PRESENCE DATA + JOIN PRESENCE AND BACKGROUND DATAFRAMES ##########

# Load libraries
library(raster)
library(maptools)
library(rgdal)
library(maps)
library(dismo)


# path to clipped layers
path2cliplayers <- "/path_2_clipped_layers"
#path to presence points
path2presence<-"./thinned_data.csv"
#path background points
path2bckgroutput<-"./biased_background_points.csv"

# path to output dataframe including both presence and background data
path2spdb<-"/path_2_species_background_file/my_presence_absence_background_file.csv"

# Load environmental variables
layers<-list.files(path=path2cliplayers,pattern='asc',full.names=TRUE)
# make the stack
egv<-stack(layers)
egv
#make a subset of the variables that you want to use
egv<-egv[[c("bio4","bio12","max_NDVI","bio17","bio2","bio5","sm_surface","soil_water_infiltration_flux","bio6","bio1","elevation","river_slope")]]
# check rasters
print(egv)
# Descriptive statistics for environmental layers
summary(egv)
# Plot environmental variables on geographical space
plot(egv)

# Load species records
data.sp<-read.csv(path2spcsv,sep=",",dec=".",header=T)
#data.sp<-shapefile(path2spshp)
data.sp<-as.data.frame(data.sp)
# Check the structure of the species dataframe
head(data.sp)
# Select coordinate fields on species records
data.sp<-data.sp[c(1,2)]
# Check the structure of the species dataframe
head(data.sp)
# Plot species records on geographical space
plot(data.sp)
map(,,,add=TRUE)

# Extract background data
# Create a cloud of random points
#backgr<-as.data.frame(randomPoints(egv,500))
# or load the background dataframe created earlier:
backgr<-read.csv(path2back,header=TRUE)
#backgr<-as.data.frame(shapefile(path2back)) #when the file is a shapefile
# Check background dataframe
head(backgr)
backgr<-backgr[,c(2,3)]
head(backgr)
names(backgr)<-c("x","y")
str(backgr)
nrow(backgr)
plot(backgr)
map(,,,add=TRUE)

# Extract environmental data on species records
data.pres<-extract(egv,data.sp)
# Check presence dataframe
str(data.pres)
# Join coordinates to presence records
data.pres<-cbind(data.pres,data.sp)
# Delete NA values from data frame
data.pres<-na.omit(data.pres)
# Check presence dataframe
str(data.pres)

# Extract environmental data on background records
data.backgr<-extract(egv,backgr)
# Delete NA values from data frame
data.backgr_noNA<-na.omit(data.backgr)
# Check background dataframe
str(data.backgr_noNA)
data.backgr_noNA
#delete the omitted NA background points also from the backgr data
backgr<-backgr[c(-8219,-9090,-9986),]
# Join coordinates to background records
data.backgr<-cbind(data.backgr_noNA,backgr)
# Check background dataframe
str(data.backgr)

# It is necessary to change the names of fields x/y to long/lat
# If not, colum names will not match
# Change number of colums depeding on the number of egv
colnames(data.backgr)[13]<-"long"
colnames(data.backgr)[14]<-"lat"
data.backgr$present<-"0"
# Check background dataframe
str(data.backgr) #kan ook met view(data.backgr) controleer of de x/y kolomnamen vervangen zijn door long/lat
# It is necessary to change the names of fields x/y to long/lat
# If not, column names will not match
# Change number of columns depending on the number of columns in egv
colnames(data.pres)[13]<-"long"
colnames(data.pres)[14]<-"lat"
data.pres$present<-"1"
# Check background dataframe
str(data.pres)#kan ook met view(data.pres)


# Join presence and background dataframes
data.pb<-data.frame(rbind(data.pres,data.backgr))
# Check presence+background dataframe
str(data.pb)
# Export presence and background dataframe
write.table(data.pb,file=path2spdb,sep=',')

# Represent species records against background
# Presences are red
par(mfrow=c(2,2))
for (j in names(egv)) {
  hist(data.backgr[,j],freq=T,main=j) # returns the density data
  hist(data.pres[,j],freq=T,col='red',add=T) # returns the density data
}
par(mfrow=c(1,1))




##############################
# Correlative model building #
##############################

# load the necessary libraries

library(biomod2)
library(raster)
library(dismo)
library(tidyterra)
library(mgcv)
library(unixtools)

# Set your paths
setwd("/my_working_directory")
wd<-"/my_working_directory"
# path to species presence/background files
path2sp<-"/path_2_species_background_file/my_presence_absence_background_file.csv"
# path to training layers
path2layers<-"/path_2_environmental_training_layers/"
# path to maxent jar file
path2maxent<-'/path_2_maxent.jar_file'
# path to output directory
path2work<-'/path_2_output_directory'
# path to future environmental layers
path2future90ssp370<-"/path_2_future_2090-ssp370"
path2future90ssp126<-"/path_2_future_2090-ssp126"
path2future50ssp370<-"/path_2_future_2050-ssp370"
path2future50ssp126<-"/path_2_future_2050-ssp126"

# path to current environmental layers covering the full map:
path2current<-"/path_2_environmental_projection_layers"

modelID<-"your_model_output_ID" # name of the individual model projections
ensemblemodelID<-"your_ensemble_model_output_ID"# name of the ensemble model projections

#change the temporary directory to a place with sufficient memory to store all the temporary files:
set.tempdir("path_2_tempdir")


########################
# Prepare all the data #
########################


setwd(wd)
# load the species data
FullDataSpecies <- read.csv(path2sp,sep=',')

# Check species data
head(FullDataSpecies)
str(FullDataSpecies)
nrow(FullDataSpecies)

#keep only the columns with the species name, x and y coordinates and place the column with the presence or absence of the species first
DataSpecies <- FullDataSpecies[,c(4:6)]
names(DataSpecies)<-c("long","lat","B.truncatus")
head(DataSpecies)

#order the data so that presences are ranked on the upper half (needed when setting the BIOMOD_formatting as only the first data in each cell will be kept)
DataSpecies$B.truncatus<-as.numeric(DataSpecies$B.truncatus)
DataSpecies<-DataSpecies[order(DataSpecies$B.truncatus,decreasing=TRUE),]
head(DataSpecies)
nrow(DataSpecies)

# the name of studied species
myRespName <- 'B.truncatus'

# load the environmental raster layers
layers<-stack(list.files(path=path2layers,pattern='asc',full.names=TRUE ))
# Select variables to include in the model
names(layers)
myExpl<-layers[[c("bio4","max_NDVI","bio17","bio2","bio5","bio6","bio1","bio12","river_slope")]]

# check rasters
myExpl

# Descriptive statistics for environmental layers
summary(myExpl$max_NDVI)
# Plot environmental variables on geographical space
plot(myExpl$max_NDVI)

# Define the presence/absence data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("long","lat")]

# Format the data so it can be used by BIOMOD2
myBiomodData <- BIOMOD_FormatingData(
  resp.var = myResp, 
  expl.var = myExpl, 
  resp.xy = myRespXY, 
  resp.name = myRespName,
  filter.raster=F
  )

# print_formating_data
myBiomodData

#: plot_formating_data
plot(myBiomodData)

#####################################
# Defining MAXENT Modelling options #
#####################################

# Consult help
?BIOMOD_ModelingOptions

myBiomodOption <- BIOMOD_ModelingOptions(MAXENT= list(path_to_maxent.jar=path2maxent, 
                         memory_allocated=1024,
                         maximumiterations = 500, 
                         visible = TRUE, 
                         linear = TRUE,
                         quadratic = TRUE,
                         product = TRUE, 
                         threshold = TRUE, 
                         hinge = FALSE, 
                         lq2lqptthreshold = 80, 
                         l2lqthreshold = 10, 
                         hingethreshold = 15, 
                         beta_threshold = 1, 
                         beta_categorical = 1, 
                         beta_lqp = 1, 
                         beta_hinge = 1, 
                         betamultiplier = 1,
                         defaultprevalence = 0.5),
                         #GLM settings
                         GLM=list(type="polynomial",
                                  interaction.level=1,
                                  test="AIC",
                                  control=glm.control(epsilon=1e-10, 
                                                      maxit=1000, 
                                                      trace=FALSE)), 
                         #GBM settings
                         GBM=list(n.trees=100,
                         interaction.depth=2,
                         n.minobsinnode=5,
                         shrinkage=0.1,
                         bag.fraction=0.5,
                         train.fraction=1,
                         cv.folds=3,
                         keep.data=F,
                         verbose=F,
                         n.cores=3),
                         #GAM settings
                         GAM=list(algo="GAM_mgcv",
                                  control=gam.control(efs.tol="REML"),
                                  k=4
                                  )
                         )


########################
# Computing the models #
########################

# Setting working directory
setwd(path2work)

# See help
?BIOMOD_Modeling
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                     models = c('GAM','MAXENT',"GBM","GLM"), 
                                     bm.options = myBiomodOption, 
                                     CV.strategy="random",
                                     CV.nb.rep =100, # number of evaluation runs
                                     CV.perc=0.7, # % of data used as calibration
                                     var.import=4, # number of permutations to estimate variable importance
                                     metric.eval = c("KAPPA","TSS","ROC"), # evaluation metrics
                                     save.output = TRUE,# keep all results and outputs on hard drive
                                     scale.models = FALSE,
                                     modeling.id = modelID,
                                     nb.cpu=3, # the number of cores available on your machine -1
                                     do.progress=TRUE)

save(myBiomodModelOut,file="myBiomodModelOut.RData")
#load the saved results from the file:
load("myBiomodModelOut.RData")

# modeling_summary
myBiomodModelOut

# get all models evaluation 
myBiomodModelEval<-get_evaluations(myBiomodModelOut) 
myBiomodModelEval
write.table(myBiomodModelEval, file="individual_model_evaluation.txt")

model_eval_metric<-myBiomodModelEval[myBiomodModelEval$full.name %in% c("B.truncatus_allData_allRun_GLM","B.truncatus_allData_allRun_MAXENT","B.truncatus_allData_allRun_GBM","B.truncatus_allData_allRun_RF","B.truncatus_allData_allRun_GAM"),c(1,5,7:9)]
model_eval_metric
write.table(model_eval_metric,file="average_individual_model_evaluation.txt")
model_evalsens<-aggregate(myBiomodModelEval$sensitivity,list(myBiomodModelEval$algo),FUN=mean)
names(model_evalsens)<-c("algorithm","mean_sensitivity")
model_evalspec<-aggregate(myBiomodModelEval$specificity,list(myBiomodModelEval$algo),FUN=mean)
names(model_evalspec)<-c("algorithm","mean_specificity")
model_evalcali<-aggregate(myBiomodModelEval$calibration,list(myBiomodModelEval$algo),FUN=mean)
names(model_evalspec)<-c("algorithm","mean_calibration")
model_eval<-cbind(model_evalsens[,c(1,2)],model_evalspec[,2],model_evalcali[,2])
names(model_eval)<-c("algorithm","mean_specificity","mean_sensitivity","mean_calibration")
model_eval
bm_PlotEvalBoxplot(myBiomodModelOut,dataset="calibration",group.by=c("algo","run"),do.plot=TRUE)


#compare the different algorithms:
bm_PlotEvalMean(bm.out=myBiomodModelOut, dataset="calibration")
bm_PlotEvalMean(bm.out=myBiomodModelOut, dataset="validation")
bm_PlotEvalMean(bm.out=myBiomodModelOut,group.by="algo")

# print the ROC scores of all selected models
myBiomodModelEval["sensitivity"]

# print variable importances
myBiomodModel_var_import<-get_variables_importance(myBiomodModelOut)
write.table(myBiomodModel_var_import,file="variable_importance_all_models.txt")
myBiomodModel_var_import[,c(4,5,7)]

var_import<-aggregate(myBiomodModel_var_import$var.imp,list(myBiomodModel_var_import$expl.var),FUN=mean)
names(var_import)<-c("variable","mean_importance")
var_import
write.table(myBiomodModel_var_import,file="average_variable_importance_all_models.txt")

# print response curves for each variable 
eval_strip<-bm_PlotResponseCurves(
  myBiomodModelOut,
  models.chosen=c("B.truncatus_allData_allRun_GAM","B.truncatus_allData_allRun_GLM","B.truncatus_allData_allRun_MAXENT","B.truncatus_allData_allRun_GBM"),
  new.env=get_formal_data(myBiomodModelOut,"expl.var"),
  show.variables=get_formal_data(myBiomodModelOut,"expl.var.names"),
  fixed.var="mean",
  do.plot=T
)


# save models evaluation scores and variables importance on hard drive
capture.output(get_evaluations(myBiomodModelOut),
               file=file.path(myRespName, 
                              paste(myRespName,"_formal_models_evaluation.txt", sep="")))

capture.output(get_variables_importance(myBiomodModelOut),
               file=file.path(myRespName, 
                              paste(myRespName,"_formal_models_variables_importance.txt", sep="")))  


#############################
# Construct ensemble models #
#############################

B.t_ensemble_models<-BIOMOD_EnsembleModeling(bm.mod=myBiomodModelOut,
                                             models.chosen="all",
                                             metric.select ="ROC",
                                             metric.select.thresh = c(0.7),
                                             metric.eval=c("TSS","ROC","KAPPA"),
                                             em.by="all",
                                             em.algo=c("EMmean","EMca","EMcv"),
                                             var.import=4, 
                                             nb.cpu=3, 
                                             do.progress=TRUE)

B.t_ensemble_models
save(B.t_ensemble_models,file="B.t_ensemble_models.RData")
#load the saved results from the file:
load("B.t_ensemble_models.RData")

# Evaluate the models

myEMBiomodModelEval<-get_evaluations(B.t_ensemble_models) 
myEMBiomodModelEval
write.table(myEMBiomodModelEval,file="ensemble_model_evaluation")
myEMBiomodModelEval[myEMBiomodModelEval$metric.eval %in% c("ROC"),c(1,8:11)]
myEMBiomodModel_var_import<-get_variables_importance(B.t_ensemble_models)
EMvar_import<-aggregate(myEMBiomodModel_var_import$var.imp,list(myEMBiomodModel_var_import$expl.var),FUN=mean)
names(EMvar_import)<-c("variable","mean_importance")
EMvar_import
write.table(EMvar_import,file="ensemble_model_variable_importance")

# print response curves for each variable 
EMeval_strip<-bm_PlotResponseCurves(
  B.t_ensemble_models,
  models.chosen="all",
  new.env=get_formal_data(B.t_ensemble_models,"expl.var"),
  show.variables=get_formal_data(B.t_ensemble_models,"expl.var.names"),
  fixed.var="mean",
  do.plot=T
)


#########################################
# Project the model over the study area #
#########################################

# environmental variables covering the same extent as the future projections:
current<-stack(list.files(path=path2current,pattern="asc",full.names=TRUE))
names(current)

myExplcurrent<-current[[c("bio4","max_NDVI","bio17","bio2","bio5","bio6","bio1","bio12","river_slope")]]
myExplcurrent

# project the ensemble model over the study area:
myEMBiomodProj <- BIOMOD_EnsembleForecasting(
  bm.em = B.t_ensemble_models,# modelling results
  new.env = myExplcurrent, # environmental variables 
  proj.name = ensemblemodelID, # name of the projections
  models.chosen="all",# the models to project
  metric.binary = 'all', # a vector of a subset of models evaluation method computed before
  output.format = '.img') # the format of the GIS files: .RData, .grd or .img


plot(myEMBiomodProj)
output<-raster("proj_ensemble_B.truncatus_ensemble.img")
plot(output)
writeRaster(output,file="B.truncatus_ensemble",format="GTiff") # save the projection as a geotiff file

save(myEMBiomodProj,file="myEMBiomodProj.RData")
#load the saved results from the file:
load("myEMBiomodProj.RData")

##################################
# Project to model to the future #
##################################

# load the future environmental data layers:

future50370<-stack(list.files(path=path2future50ssp370,pattern='asc',full.names=TRUE ))
future50126<-stack(list.files(path=path2future50ssp126,pattern='asc',full.names=TRUE ))
future90370<-stack(list.files(path=path2future90ssp370,pattern='asc',full.names=TRUE ))
future90126<-stack(list.files(path=path2future90ssp126,pattern='asc',full.names=TRUE ))

# select the different climate models to use:

models<-c("BCC.CSM2.MR","IPSL.CM6A.LR","MIROC6","MRI.ESM2.0") 


# the variables to include in the model:
variables<-c("bio4","max_NDVI","bio17","bio2","bio5","bio6","bio1","bio12","river_slope") 
#the climate change scenario to run the model for
scenario50370<-"future50370"
scenario50126<-"future50126"
scenario90370<-"future90370"
scenario90126<-"future90126"


#50370
  for (model in models){ #loop over all models
    modvarlist<-c() # make an empty list
    for (variable in variables) { #loop over all variables of interest
      modvar<-paste(variable,model, sep="_") #construct new names with the variable and the model name
      modvarlist<-c(modvarlist,modvar) #add each new name to the list
      modvarlist<-gsub("max_NDVI.*","max_NDVI",modvarlist) #look for the pattern "max_NDVI" in the list and only keep the max_NDVI part 
      modvarlist<-gsub("river_slope.*","river_slope",modvarlist) #look for the pattern "river_slope" in the list and only keep the river_slope part 
       }
    print(modvarlist[1:9])
    modelID<-paste(scenario50370,model,sep="_")
    myExpl<-future50370[[modvarlist]]
    names(myExpl)<-variables
    futureID<-paste("future50370",model,sep="_")
    print(futureID)
    outputfilename<-paste(futureID,"RData",sep=".")
    myEMBiomodProj <- BIOMOD_EnsembleForecasting(
      bm.em = B.t_ensemble_models,# modelling results
      new.env = myExpl, # environmental variables 
      proj.name = futureID, # name of the projections
      models.chosen="all",
      metric.binary = 'all',
      output.format = '.img') 
    save(myEMBiomodProj,file=outputfilename)
}



#50126
for (model in models){
  modvarlist<-c()
  for (variable in variables) {
    modvar<-paste(variable,model, sep="_") 
    modvarlist<-c(modvarlist,modvar)
    modvarlist<-gsub("max_NDVI.*","max_NDVI",modvarlist) 
    modvarlist<-gsub("river_slope.*","river_slope",modvarlist) 
    }
  modelID<-paste(scenario50126,model,sep="_")
  myExpl<-future50126[[modvarlist]]
  names(myExpl)<-variables
  futureID<-paste("future50126",model,sep="_")
  print(futureID)
  outputfilename<-paste(futureID,"RData",sep=".")
  myEMBiomodProj <- BIOMOD_EnsembleForecasting(
    bm.em = B.t_ensemble_models,
    new.env = myExpl, 
    proj.name = futureID, 
    models.chosen="all",
    metric.binary = 'all', 
    output.format = '.img') 
  save(myEMBiomodProj,file=outputfilename)
}

#90126
for (model in models){
  modvarlist<-c()
  for (variable in variables) {
    modvar<-paste(variable,model, sep="_") 
    modvarlist<-c(modvarlist,modvar)
    modvarlist<-gsub("max_NDVI.*","max_NDVI",modvarlist) 
    modvarlist<-gsub("river_slope.*","river_slope",modvarlist) 
  }
  modelID<-paste(scenario90126,model,sep="_")
  myExpl<-future90126[[modvarlist]]
  names(myExpl)<-variables
  futureID<-paste("future90126",model,sep="_")
  print(futureID)
  outputfilename<-paste(futureID,"RData",sep=".")
  myEMBiomodProj <- BIOMOD_EnsembleForecasting(
    bm.em = B.t_ensemble_models,
    new.env = myExpl, 
    proj.name = futureID, 
    models.chosen="all",
    metric.binary = 'all', 
    output.format = '.img') 
  save(myEMBiomodProj,file=outputfilename)
}

#90370
for (model in models){
  modvarlist<-c()
  for (variable in variables) {
    modvar<-paste(variable,model, sep="_") 
    modvarlist<-c(modvarlist,modvar)
    modvarlist<-gsub("max_NDVI.*","max_NDVI",modvarlist) 
    modvarlist<-gsub("river_slope.*","river_slope",modvarlist) 
  }
  modelID<-paste(scenario90370,model,sep="_")
  myExpl<-future90370[[modvarlist]]
  names(myExpl)<-variables
  futureID<-paste("future90370",model,sep="_")
  print(futureID)
  outputfilename<-paste(futureID,"RData",sep=".")
  myEMBiomodProj <- BIOMOD_EnsembleForecasting(
    bm.em = B.t_ensemble_models,
    new.env = myExpl, 
    proj.name = futureID, 
    models.chosen="all",
    metric.binary = 'all', 
    output.format = '.img') 
  save(myEMBiomodProj,file=outputfilename)
}

