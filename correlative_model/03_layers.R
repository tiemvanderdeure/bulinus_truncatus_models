###################################
# Prepare environmental variables #
###################################
source("paths.R")

########################## CLIPPING RASTERS #############################
library(geodata)
library(terra)

# Load environmental variables
bioclim <- geodata::worldclim_global(var = "bio", res = "2.5", path = path2layers)
bioclim <- crop(bioclim, ext(-20,60,-40,60))

#### To get all layers, also run both .ipynb notebooks and download elevation from worldclim
layers <- list.files(path=path2layers,pattern=c(".tif$|.asc$"),full.names=TRUE)
layers <- layers[!grepl("clipped|cropped", layers)]
# Join all rasters in a single object
egv <- bioclim
for (layer in layers){
  print(layer)
    # resample to bioclim to make sure all rasters align
  new <- resample(rast(layer),bioclim,method="bilinear")
  egv <- c(egv, new)
}

names(egv) <- gsub("wc2.1_2.5m_", "", names(egv))

# Load clipping area shapefile
shp<-vect(path2clipshp)

cropped <- crop(egv, shp)
masked <- mask(egv, shp)

# write the rasters
writeRaster(cropped, path2current, overwrite = TRUE)
writeRaster(masked, path2clippedlayers, overwrite = TRUE)
