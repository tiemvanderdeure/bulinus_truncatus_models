#############################
## Life history trait fits ##
#############################

## Code by: Tiem van der Deure, tvd@sund.ku.dk

# Load libraries
library(R2jags)
library(data.table)

# Set seed for reproducibility
set.seed(1234)

#####################
#### Read data ######
#####################

# Read in and clean data collected for hatching and maturation
traits_data <- fread(
    "trait_fits/data/traits_data.csv", 
    sep = ";", dec = ","
)

## Read in cleaned survival data
survival_data <- fread("trait_fits/data/survival_data.csv")
# Read in cleaned egg laying data
clutch_size <- fread(file = "trait_fits/data/clutch_size.csv")
n_clutches <- fread(file = "trait_fits/data/number_of_clutches.csv")

#####################
### JAGS settings ###
#####################
# These are the same for all model runs
nthin <- 8
nchains <- 5
nburnin <- 5000
niter <- 25000

############################################
# Set up and run JAGS model for each trait #
############################################

##########################
#### Egg-laying rate #####
##########################

### Little bit of data cleaning
populations <- unique(n_clutches[, .(country, population_id)])
experiments <- unique(n_clutches[, .(temperature, country, experiment_id, population_id)])
snails <- unique(n_clutches[, .(individual, experiment_id, snail_id)])


## Choose inital values
egg_inits <- function(){
  list(
    q_mu = 0.3,
    Tm_mu = 34,
    T0_mu = 12,
    q_sigma = 0.05,
    T0_sigma = 1,
    Tm_sigma = 1,
    lambda_neggs = 1,
    lambda_cs = 10
  )
}

## Bundle data
eggs_data <- list(
  n_clutches = n_clutches$egg_masses,
  clutch_size = clutch_size$clutch_size,
  N.records = nrow(n_clutches),
  N.snails = nrow(snails),
  N.experiments = nrow(experiments),
  N.populations = nrow(populations),
  N.obs_cs = nrow(clutch_size),
  temp = experiments$temperature,
  exp_id = snails$experiment_id,
  snail_id = n_clutches$snail_id,
  snail_id_cs = clutch_size$snail_id,
  pop_id = experiments$population_id
)

## Run the model
eggs_fit <- jags(
  data = eggs_data,
  inits = egg_inits,
  model.file = "trait_fits/egg_laying.txt",
  parameters.to.save = c("q_mu", "T0_mu", "Tm_mu"),
  n.thin = nthin,
  n.chains = nchains,
  n.burnin = nburnin,
  n.iter = niter
)


##########################
######## Lifespan ########
##########################

## Experiments by temperature
temperatures <- unique(survival_data[, .(temperature, experiment_index, population_index)])

## Set initial values
survival_inits <- function(){list(
  lf_max_mu = 5,
  Tm_mu = 35,
  T0_mu = 5,
  lf_max_sigma = 0.5,
  T0_sigma = 1,
  Tm_sigma = 1,
  temp_sigma = 1,
  weeks_survived = ifelse(is.na(survival_data$weeks_survived) & (survival_data$source == ""), survival_data$exp_length+survival_data$survived, NA)
)}

## Bundle data
survival_jags_data <-list(
    survived = survival_data$survived,
    weeks_survived = survival_data$weeks_survived,
    exp_index = survival_data$experiment_index,
    exp_length = survival_data$exp_length,
    sample_size = survival_data$sample_size,
    N_obs_raw = nrow(survival_data[is.na(sample_size)]), 
    N_obs_total = nrow(survival_data),
    N_groups = nrow(temperatures),
    N_populations = max(survival_data$population_index),
    temp = temperatures$temperature,
    pop_index = temperatures$population_index
)

## Run the model
survival_fit <- jags(
  data = survival_jags_data,
  inits = survival_inits,
  parameters.to.save = c("q_mu", "T0_mu", "Tm_mu"),  # parameters to save are init parameters, plus predictions
  model.file = "trait_fits/survival.txt",
  n.thin = nthin,
  n.chains = nchains,
  n.burnin = nburnin,
  n.iter = niter
)


##########################
#### Maturation time #####
##########################

# intial values
mattime_inits <- function(){list(
  T0 = 10.1,
  GDD = 150,
  sigma = rlnorm(1))}

# Priors
mattime_priors <- t(matrix(data = c(4,12,30,200), nrow = 2, ncol = 2))

# Observations
maturing_time_data <- traits_data[parameter == "maturation_time"]

##### Bundle Data
mattime_jags_data <-list(
  trait = maturing_time_data$trait,
  temp = maturing_time_data$temperature, 
  N_obs = nrow(maturing_time_data), 
  priors = mattime_priors
)

maturation_fit <- jags(
    data = mattime_jags_data,
    inits = mattime_inits,
    parameters.to.save = c("T0", "GDD"),
    model.file = "trait_fits/gdd.txt",
    n.thin = nthin,
    n.chains = nchains,
    n.burnin = nburnin,
    n.iter = niter
)


##########################
#### Hatching success ####
##########################

# intial values
hatching_success_inits <-function(){list(
  q = 0.5,
  sigma = 1)}

# priors
hatching_success_priors <- matrix(data = c(0,1, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)

# Data
hatching_success_data <- traits_data[parameter == "hatching_success"]

##### Bundle Data
hatching_success_jag_data <-list(
  trait = hatching_success_data$trait,
  N_obs = nrow(hatching_success_data), 
  priors = hatching_success_priors
)

# Run jags!!
hatching_success_fit <- jags(
    data = hatching_success_jag_data,
    inits = hatching_success_inits,
    parameters.to.save = c("q"),
    model.file = "trait_fits/constant.txt",
    n.thin = nthin,
    n.chains = nchains,
    n.burnin = nburnin,
    n.iter = niter
)

##########################
##### Hatching time ######
##########################

# intial values
hatching_time_inits <- function(){list(
  T0 = 8.1,
  GDD = 150,
  sigma = rlnorm(1))}

# Priors
hatching_time_priors <- t(matrix(data = c(4,12,50,300), nrow = 2, ncol = 2))

# Data
hatching_time_data <- traits_data[parameter == "hatching_time"]

##### Bundle Data
hatching_time_jag_data <-list(
  trait = hatching_time_data$trait,
  temp = hatching_time_data$temperature,
  N_obs = nrow(hatching_time_data), 
  priors = hatching_time_priors
)

# Run jags!!
hatching_time_fit <- jags(
    data = hatching_time_jag_data,
    inits = hatching_time_inits,
    parameters.to.save = c("T0", "GDD"),
    model.file = "trait_fits/gdd.txt",
    n.thin = nthin,
    n.chains = nchains,
    n.burnin = nburnin,
    n.iter = niter
)

###### EXPORT RESULTS
models <- list(hatching_time_fit, eggs_fit, maturation_fit, hatching_success_fit, survival_fit)
filenames = c("hatchingtime", "egglaying", "maturation", "hatchingsuccess", "lifespan")

# Export mean estimates
params_not_to_save <- c("sigma", "pred", "deviance")

models_summary <- data.table(trait = character(), parameter = character(), mean = numeric(), sd = numeric(), "2.5%" = numeric(), "97.5%" = numeric())

for (i in 1:5){
  model <- models[[i]]

  params <- model$parameters.to.save[!sapply(model$parameters.to.save, function(x) x %in% params_not_to_save)]

  filename <- paste0("trait_fits/fits/", filenames[i],".csv")

  model_summary <- data.table(model$BUGSoutput$summary[params,c("mean", "sd", "2.5%", "97.5%"), drop = FALSE])
  model_summary$trait <- filenames[[i]]
  model_summary$parameter <- params

  models_summary <- rbind(models_summary, model_summary)

  write.csv(model$BUGSoutput$summary[params,c("mean", "sd", "2.5%", "97.5%"), drop = FALSE], filename)
}

write.csv(models_summary, "trait_fits/fits/param_estimates.csv", row.names = FALSE)

# Export 5000 draws for each parameter
iterations <- 1000
draws_dt <- data.table(index = c(1:(iterations*nchains))) # get 5 times iterations samples, because there are 5 chains

for (i in 1:5){
  model <- models[[i]]
  name <- filenames[i]

  parameters <- model$parameters.to.save[!sapply(model$parameters.to.save, function(x) x %in% params_not_to_save)]

  samples <- jags.samples(model$model, variable.names = parameters, n.iter = iterations)

  for (param in parameters){
    param_name <- paste(name, gsub("_mu", "", param), sep = "_")
    draws_dt[, param_name] <- samples[[param]]
  }
}

write.csv(draws_dt[, -"index"], "trait_fits/fits/draws_from_posterior.csv", row.names = FALSE)

#######################
#### Plotting data ####
#######################

# Here, calculate some naive means without Bayesian analysis to plot in the figure

### Egg-laying
mean_n_clutches <- n_clutches[, .(n_clutches = mean(egg_masses)), by = c("temperature", "country")]
mean_clutch_size <- clutch_size[, .(clutch_size = mean(clutch_size)), by = c("temperature", "country")]

# Merge number of clutches and clutch size
egg_means <- n_clutches[, .(n_clutches = mean(egg_masses), N_n_clutches = .N), by = c("temperature", "country")][
  clutch_size[, .(clutch_size = mean(clutch_size), N_clutch_size = .N), by = c("temperature", "country")], on = c("temperature", "country")
]

# calculate number of eggs and effective sample size (geometric mean)
egg_means[, `:=`(trait =  n_clutches * clutch_size, sample_size = sqrt(N_n_clutches * N_clutch_size))]
egg_means_to_plot = egg_means[, .(temperature, country, trait, sample_size, parameter = "eggs")]

### Lifespan
# Get the last week where snails from a group were alive and find the death rate assuming exponential decay
surv_means <- survival_data[
  source == "", .(
    len = max(ifelse(is.na(weeks_survived), exp_length+as.integer(1), weeks_survived)), 
    surv = sum(survived),
    sample_size = .N), 
  by = c("temperature", "country")]

surv_means[, trait := 1/(-log(surv/sample_size)/len)]
surv_means_to_plot = surv_means[, .(temperature, country, trait, sample_size, parameter = "lifespan")]

surv_means_other_studies = survival_data[source != "", .(temperature, country = population_index, trait = 1/(-log(survived/sample_size)/exp_length), sample_size, parameter = "lifespan")]

data_points_plot <- rbind(egg_means_to_plot, surv_means_to_plot, surv_means_other_studies, maturing_time_data, hatching_success_data, hatching_time_data,
  fill = TRUE)

write.csv(data_points_plot, "trait_fits/data/data_points_plot.csv", row.names = FALSE, na = "")

### Plotting itself is done in Julia

