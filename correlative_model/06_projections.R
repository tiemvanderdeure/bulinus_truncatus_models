#########################################
# Project the model over the study area #
#########################################
source("paths.R")

# load the necessary libraries
library(biomod2)
library(dplyr)
library(mgcv)
library(terra)
library(doParallel)

# set number of cores for parallel computing
ncores = detectCores() - 1

vars <- paste0("bio_", c(1,2,4,5,6,9,12,17,18,19))
set <- "set2"

# environmental variables covering the same extent as the future projections:
current<-rast(path2current)

ensemble <- load(paste0("B.t_ensemble_models_", set, ".RData"))

# project the ensemble model over the study area:
myEMBiomodProj <- BIOMOD_EnsembleForecasting(
  bm.em = ensemble,# modelling results
  new.env = current[[vars]], # environmental variables 
  proj.name = ensemblemodelID, # name of the projections
  models.chosen="all",# the models to project
  metric.binary = 'all', # a vector of a subset of models evaluation method computed before
  output.format = '.img',
  nb.cpu = ncores
) # the format of the GIS files: .RData, .grd or .img

plot(myEMBiomodProj)
writeRaster(myEMBiomodProj,file="B.truncatus_ensemble.tif") # save the projection as a geotiff file




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

#the climate change scenario to run the model for
scenario50370<-"future50370"
scenario50126<-"future50126"
scenario90370<-"future90370"
scenario90126<-"future90126"

#50370
for (model in models){ #loop over all models
  modelID<-paste(scenario50370,model,sep="_")
  myExpl<-future50370[[vars]]
  futureID<-paste("future50370",model,sep="_")
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
  modelID<-paste(scenario50126,model,sep="_")
  myExpl<-future50126[[vars]]
  futureID<-paste("future50126",model,sep="_")
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
  modelID<-paste(scenario90126,model,sep="_")
  myExpl<-future90126[[vars]]
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
  modelID<-paste(scenario90370,model,sep="_")
  myExpl<-future90370[[vars]]
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
