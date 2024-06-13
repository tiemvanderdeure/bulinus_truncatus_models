# If needed, change the working directory here. Scripts assume the working directory is the 'correlative_model' folder that these scripts are in.
# setwd("correlative_model")

#### Paths to define manually!
# path to maxent jar file
path2maxent<-"~/Desktop/maxent/maxent.jar"


#### Other paths
path2layers <- "env_layers" # path to folder with all environmental variables of interest
path2clipshp<-"roi/clipping_area.shp"
path2clippedlayers <- paste0(path2layers, "/clipped.tif")
# path to current environmental layers covering the full map:
path2current<-"env_layers/cropped.tif"

#path to presence points
path2allpresences<-"../occurrence_data/all_presences.csv"
#path to presence points
path2presence<-"../occurrence_data/thinned_presences.csv"
#path background points
path2background<-"../occurrence_data/background_points.csv"
# path to output dataframe including both presence and background data
# path2spdb<-"/path_2_species_background_file/my_presence_absence_background_file.csv"

# path to output directory
path2work<-'outputs/'
# path to future environmental layers
path2future90ssp370<-"/path_2_future_2090-ssp370"
path2future90ssp126<-"/path_2_future_2090-ssp126"
path2future50ssp370<-"/path_2_future_2050-ssp370"
path2future50ssp126<-"/path_2_future_2050-ssp126"
