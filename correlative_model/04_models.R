##############################
# Correlative model building #
##############################
source("paths.R")

# load the necessary libraries
library(biomod2)
library(dplyr)
library(mgcv)
library(terra)
library(doParallel)

# set number of cores for parallel computing
ncores = detectCores() - 1


sets <- list(
  c("max_NDVI","river_slope", paste0("bio_", c(1,2,4,5,6,9,12,17))),
  paste0("bio_", c(1,2,4,5,6,9,12,17,18,19))
)
modelIDs <- c("set1", "set2") # name of the individual model projections
ensemblemodelIDs<- paste0(modelIDs, "_ensemble")# name of the ensemble model projections

#change the temporary directory to a place with sufficient memory to store all the temporary files:
#set.tempdir("path_2_tempdir")

########################
# Prepare all the data #
########################
# Load environmental variables
egv <- rast(path2clippedlayers)

# Load species records
data.sp<-read.csv(path2presence,sep=",",dec=".",header=T)
data.sp[,1] <- NULL
data.sp$present <- 1
# Extract background data
backgr<-read.csv(path2background,header=TRUE)
backgr[,1] <- NULL
backgr$present <- 0
# bind records and background
pr_bg_points <- vect(rbind(data.sp, backgr), geom = c("x", "y"))
truncatus <- extract(egv, pr_bg_points, bind = TRUE)
myResp <- truncatus$present # for some reasons this needs to be defined for BIOMOD_PresenceOnly to work here

#####################################
# Defining MAXENT Modelling options #
#####################################

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
                                                  n.cores=1),
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

for(i in c(1,2)){
  
  set <- sets[[i]]
  modelID <- modelIDs[i]
  ensemblemodelID <- ensemblemodelIDs[i]
  
  myBiomodData <- BIOMOD_FormatingData(
    resp.var = truncatus$present, 
    expl.var = truncatus[, set], 
    resp.name = "B.truncatus",
    resp.xy = geom(truncatus)[, c("x", "y")],
  )

  # Run the models
  myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
                                       models = c("MAXENT", "GLM", "GAM", "GBM"), 
                                       bm.options = myBiomodOption, 
                                       CV.strategy="random",
                                       CV.nb.rep =100, # number of evaluation runs
                                       CV.perc=0.7, # % of data used as calibration
                                       var.import=4, # number of permutations to estimate variable importance
                                       metric.eval = c("KAPPA","TSS","ROC"), # evaluation metrics
                                       #save.output = TRUE,# keep all results and outputs on hard drive
                                       scale.models = FALSE,
                                       modeling.id = modelID,
                                       nb.cpu=ncores, # the number of cores available on your machine -1
                                       do.progress=TRUE)

  # save the results
  save(myBiomodModelOut,file=paste0("myBiomodModelOut_set", i, ".RData"))
  # get all models evaluation 
  myBiomodModelEval<-biomod2::get_evaluations(myBiomodModelOut)

  write.csv(myBiomodModelEval, paste0("ensemble_evaluation_set", i, ".csv"))

  model_eval_metric <- myBiomodModelEval %>% group_by(algo, metric.eval) %>% summarise_at(vars(cutoff, calibration, validation), list( mean=  mean))
  write.csv(model_eval_metric, paste0("average_individual_model_evaluation_set", i, ".csv"))

  # print variable importances
  myBiomodModel_var_import<-get_variables_importance(myBiomodModelOut)
  write.csv(myBiomodModel_var_import,file=paste0("variable_importance_all_models_set", i, ".csv"))

  var_import<-aggregate(myBiomodModel_var_import$var.imp,list(myBiomodModel_var_import$expl.var),FUN=mean)
  names(var_import)<-c("variable","mean_importance")
  write.csv(var_import,file=paste0("average_variable_importance_all_models_set", i, ".csv"))

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
                                               nb.cpu=ncores, 
                                               do.progress=TRUE)

  save(B.t_ensemble_models,file= paste0("B.t_ensemble_models_set", i, ".RData"))

  # Evaluate the models
  myEMBiomodModelEval<-biomod2::get_evaluations(B.t_ensemble_models)

  ### Additionally calculate the continuous boyce index
  pred <- get_predictions(B.t_ensemble_models)
  predmean <- filter(pred, algo == "EMmean")
  predca <- filter(pred, algo == "EMca")
  boyce_mean <- ecospat::ecospat.boyce(predmean$pred, predmean$pred[(myResp == 1)], method = "pearson")$cor
  boyce_ca <- ecospat::ecospat.boyce(predca$pred, predca$pred[(myResp == 1)], method = "pearson")$cor
  boyce_df <- cbind(
    rbind(predmean[1,1:6], predca[1,1:6]),
    data.frame(metric.eval = "BOYCE", cutoff = NA, sensitivity = NA, specificity =NA, 
               calibration = c(boyce_mean, boyce_ca), validation = NA, evaluation = NA)
  )
  eval_with_boyce <- rbind(myEMBiomodModelEval, boyce_df)
  
  write.csv(eval_with_boyce,file=paste0("ensemble_model_evaluation_set", i, ".csv"), row.names = FALSE)

  # variable importances
  myEMBiomodModel_var_import<-get_variables_importance(B.t_ensemble_models)
  EMvar_import<-aggregate(myEMBiomodModel_var_import$var.imp,list(myEMBiomodModel_var_import$expl.var),FUN=mean)
  
  names(EMvar_import)<-c("variable","mean_importance")
  write.csv(EMvar_import,file=paste0("ensemble_model_variable_importance_set", i, ".csv"))
}
  
# print response curves for each variable 
# EMeval_strip<-bm_PlotResponseCurves(
#   B.t_ensemble_models,
#   models.chosen="all",
#   new.env=get_formal_data(B.t_ensemble_models,"expl.var"),
#   show.variables=get_formal_data(B.t_ensemble_models,"expl.var.names"),
#   fixed.var="mean",
#   do.plot=T
# )
