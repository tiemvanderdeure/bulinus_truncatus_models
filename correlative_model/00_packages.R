# Before running

install.packages("devtools")
require(devtools)
# install biomod2 version 4.2.4 specifically. Version 4.2.5 deprecated functions that later scripts use.
install_version("biomod2", version = "4.2.4", repos = "http://cran.r-project.org") 
install.packages("terra")
install.packages("doParallel")
install.packages("dplyr")
install.packages("spThin")
install.packages("ggplot")
install.packages("rnaturalearth")