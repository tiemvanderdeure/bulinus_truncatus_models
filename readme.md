### Code repo for van der Deure, Maes, Huyse & Stensgaard (2024)
This repository contains all code used to run models and create figures for "Climate change could fuel urinary schistosomiasis transmission in Africa and Europe" by van der Deure, T., Maes, T., Huyse, T., and Stensgaard,. A-S.

In this publication, we use correlative and mechanistic modelling to predict the current and future distribution of schistosomiasis intermediate host _Bulinus truncatus_.

The code used for correlative and mechansitic modelling are found in separate folders.

To reproduce the correlative modelling, first install the necessary packages using the script provided. Note that we used biomod2 version 4.2.4. Then download environmental layers by running `03_layers.R`, fit the models using `04_models.R` and run the projections using `06_projections.R`. All model settings used are specified in `04_models.R`. The other scripts in this folder apply spatial thinning to the occurrence points, build background layers, and investigate correlations between environmental variables.

Most mechanistic modelling was performed in Julia. The packages used including their versions are provided in `Manifest.toml`. Use `Pkg.instantiate()` to install all packages. Run `climate_data.jl` to download the necessary climate layers. `run_model.jl` generates model runs, calling the other scripts as required.

Most of the life history trait data used in the mechanistic models stem from an experiment that has been separately published as a [preprint](https://www.biorxiv.org/content/10.1101/2024.01.02.573866v1).

Anyone is welcome to use, adapt, or redistribute (parts of) the code provided here. When using this code as part of a scientific publication, please cite the forthcoming publication.

The corresponding author for the publication is Tiem van der Deure, email tvd@sund.ku.dk.


