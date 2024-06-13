##################### ASSESSING CORRELATION AMONG VARIABLES ########################
source("paths.R")

# Load libraries
library(terra)
library(maptools)
library(rgdal)

# Load environmental variables
egv <- rast(path2clippedlayers)

# Calculate Pearson's correlation among variables
# Pearson's correlation is parametric
pearson<-layerCor(egv[[c("bio_1", "bio_2")]],'pearson')
# Check correlation results
print(pearson)

# Transform correlation results in a dataframe
cor.df<-as.data.frame(pearson$correlation)
# Check dataframe structure
head(cor.df)
# Export dataframe
write.csv(cor.df,'outputs/correlations.csv') 
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

#stack the environmental layers
myExpl <- rast(path2clippedlayers)

# load our species data
FullDataSpecies <- read.csv(path2presence,sep=",",dec=".",header=T)

#keep only the columns with the species name, x and y coordinates and place the column with the presence or absence of the species first
DataSpecies <- FullDataSpecies[,c(2:3)]
names(DataSpecies)<-c("decimalLongitude","decimalLatitude")

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

