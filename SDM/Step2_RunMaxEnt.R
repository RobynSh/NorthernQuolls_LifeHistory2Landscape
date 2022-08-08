#############################################
#                   STEP 2                  #  
#           SDM: MaxEnt Modelling           #
#         Data set: Dasyurus hallucatus     #
#            Author: Robyn Shaw             #
#             Date: 11/05/2021              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

#############################
# LOAD IN REQUIRED PACKAGES #
#############################




# Install:
list.of.packages <- c("SDMtune",
                      "raster",
                      "zeallot",
                      "ENMeval",
                      "wesanderson",
                      "cowplot",
                      "pals",
                      "ggplot2",
                      "ggsn")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(SDMtune)
library(raster)
library(zeallot)
library(ENMeval)
library(wesanderson)
library(cowplot)
library(pals)
library(ggplot2)
library(ggsn)


# NOTE: Using Maxent in R
# The file maxent.jar can be downloaded here:
# https://biodiversityinformatics.amnh.org/open_source/maxent/
# Need MaxEnt version >= 3.4.1 (Phillips et al. 2017)). This file must be copied into the right folder to be available for the dismo package (Hijmans et al. 2017). Copy the file maxent.jar into the folder named java that is located inside the folder returned by the following command:
system.file(package="dismo")

# The function checkMaxentInstallation checks that Java JDK and rJava are installed, and that the file maxent.jar is in the correct folder.
checkMaxentInstallation()

# If everything is correctly configured for dismo, the following command will return the MaxEnt version (make sure that the version is >= 3.4.1):
dismo::maxent()


###################################
#  Load environmental Variables   #
###################################

# Load up rasters
# Note: rasters have been prepared/manipulated prior to this step
# Some rasters were derived from other layers, and they all have been aggregated to the same resolution based on what was deemed most biologically appropriate (mean, median, min, max, etc.), and cropped/masked to the same extent.
Rasters <- stack(list.files("Rasters_Shapefiles/", pattern = ".asc", full.names = TRUE))

# Remove inappropriate layers (i.e. rasters used to test different hypotheses later on)
Rasters <- Rasters[[-which(names(Rasters) %in% c("Dist2Water.1km", "Dist2Water.5km", "Fortescue", "MajorRivers", "RiverBasins", "VRM.5", "VRM.10", "VRM.20"))]]

# Create vectors with all names, and formatted names for labels
Rasters.AllNames <- names(Rasters)
AllNames <- c("ADI", "ADM", "ADX", "B03", "B04", "B05", "B06", "B07", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17", "B29", "B30", "B31", "B32", "B33", "BPC1", "BPC2", "BPC3", "BPC4", "Clay", "Coarse Fragments", "DEM", "Distance to Water", "Fire Frequency (00-08)", "Persist. Forest Cover", "Rocky Soils", "Sand", "Silt", "Soils", "Spinifex", "VRM", "Weathering Index")


# Some of these variables may not be correlated, but still probably shouldn't be in the same model
# Make a set for each unique combination of the climate, ruggedness and soil layer types (6 sets in total)
# For use in Maxent further down

ClimatePC <- names(Rasters)[grepl(pattern = "BPC", names(Rasters))]
ClimateBIO <- c(names(Rasters)[grepl(pattern = "AD", names(Rasters))], 
                names(Rasters)[grepl(pattern = "B0", names(Rasters))], 
                names(Rasters)[grepl(pattern = "B1", names(Rasters))], 
                names(Rasters)[grepl(pattern = "B2", names(Rasters))], 
                names(Rasters)[grepl(pattern = "B3", names(Rasters))])
ContSoil <- c("Clay", "CoarseFrags", "Sand", "Silt")


# Set 1: 
# climate represented by: Bioclim variables
# Soil represented by: categorical soil layer
Rasters.Set1 <- Rasters[[-which(names(Rasters) %in% c(ClimatePC, ContSoil, "RockySoils"))]]

# Set 2: 
# climate represented by: PCA of bioclim layers
# Soil represented by: categorical soil layer
Rasters.Set2 <- Rasters[[-which(names(Rasters) %in% c(ClimateBIO, ContSoil, "RockySoils"))]]

# Set 3: 
# climate represented by: Bioclim variables
# Soil represented by: continuous soil layers
Rasters.Set3 <- Rasters[[-which(names(Rasters) %in% c(ClimatePC, "Soils_Preferred", "RockySoils"))]]

# Set 4: 
# climate represented by: PCA of bioclim layers
# Soil represented by: continuous soil layers
Rasters.Set4 <- Rasters[[-which(names(Rasters) %in% c(ClimateBIO, "Soils_Preferred", "RockySoils"))]]

# Set 5: 
# climate represented by: Bioclim variables
# Soil represented by: binary rocky layer
Rasters.Set5 <- Rasters[[-which(names(Rasters) %in% c(ClimatePC, ContSoil, "Soils_Preferred"))]]

# Set 6: 
# climate represented by: PCA of bioclim layers
# Soil represented by: binary rocky layer
Rasters.Set6 <- Rasters[[-which(names(Rasters) %in% c(ClimateBIO, ContSoil, "Soils_Preferred"))]]


#############################################
# Prepare presence and background locations #
#############################################

# Background sites

# We don't have absences so will use background sites which are points used to represent the environment across the landscape of interest. Background sites can fall on known presence and absence sites-this is legitimate. Background points will be generated using the bias layer produced in step 1. 

# First, load in bias layer
bias <- raster("SDM/Data/CWR_Bias_Layer.asc")
crs(bias) <- "+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs"
par(mfrow= c(1,1))
plot(bias)


# Extract points
# Set seed to make sure it's reproducible
set.seed(1)
# Studies suggest using a large number of 'absences' (10,000)
# Add an extra 500 to account for any that get removed because they're outside of the study area
n.bg <- 10500

# Change NA values to zeros otherwise the density can't be calculated
values(bias)[is.na(values(bias))] <- 0

# Generate background points (one from each cell), based on bias layer
bg <- xyFromCell(bias, sample(ncell(bias), n.bg, prob=values(bias)))
plot(bg)

# convert to spatial points
bg <- as.data.frame(bg)
coordinates(bg) <- ~x + y
crs(bg) <- crs(bias)
plot(bg)

# Check by plotting (use WII)
WIIpal <- colorRampPalette(c("black", "brown", "yellow"))
Biaspal <- colorRampPalette(c("white", "blue"))
plot(Rasters[[38]], legend=FALSE, col = WIIpal(10))
plot(bias, alpha=0.7, add=T, legend=FALSE, col = Biaspal(50))
plot(bg, add=T, col="black", legend=FALSE, pch=16, cex = 0.1)

# Match background points with environment at each location
Enviro.bg <- cbind(bg@coords, as.data.frame(raster::extract(Rasters, cbind(bg@coords))))
nrow(Enviro.bg)

# Remove any records that might have NA's as environmental data (i.e., fall into the water).
# Remove outlier records
if (any(is.na(rowSums(Enviro.bg)))) {
  Enviro.bg <- Enviro.bg[-which(is.na(rowSums(Enviro.bg))), ]
}

nrow(Enviro.bg)
# 10467

# Barely any were NAs, so randomly sample down to 10,000 again
Enviro.bg <- Enviro.bg[sample(10000, replace = FALSE), ]
nrow(Enviro.bg)
# 10000

# Check by plotting
plot(Enviro.bg$x, Enviro.bg$y)
plot(bg@coords)

# Save background points
# write.csv(Enviro.bg[, 1:2], "BackgroundPoints.RasterChoice/Dh_FinalBackgroundCoords.csv", row.names = FALSE)
# Important note!! Unfortunately I forgot to set the seed when I ran this originally. This means the background points used in the study are slightly different (just by chance). Therefore, I have provided the data for the background points used to generate the SDM presented in the paper, which I load in below. The process used to obtain these points was identical to that documented above.


#################################
# Load in presence/absence data #
#################################

# Presence coordinates
p_coords <- read.csv("SDM/Data/Dh_FinalPresenceCoords.csv")

# Absence coordinates
bg_coords <- read.csv("SDM/Data/Dh_FinalBackgroundCoords.csv")


##############
# Run Maxent #
##############

# Loop through the different raster sets for the the full SDM process documented here: https://consbiol-unibern.github.io/SDMtune/

for (i in 2:2) {
  
  # Clear memory
  gc()
  removeTmpFiles()
  
  # Create an SWD object
  # Before training a model, prepare the data in the correct format. An SWD object stores the species name, the coordinates of the species at presence and absence/background locations and the value of the environmental variables at the locations. The argument categorical indicates which (if any) environmental variables are categorical. The value of the environmental variables for each location is then extracted (while locations that have NA value for at least one environmental variable are excluded).
  
  # Specify raster set to run
  Set <- i
  
  CatVars <- c("Soils_Preferred", "RockySoils")
  
  # Create SWD object
  Dh.SWD <- prepareSWD(species = "Dasyurus hallucatus", p = p_coords, a = bg_coords, env = get(paste0("Rasters.Set", Set)), categorical = names(get(paste0("Rasters.Set", Set)))[names(get(paste0("Rasters.Set", Set))) %in% CatVars])
  
  # Name as object
  Rasters <- get(paste0("Rasters.Set", Set))
  
  # Create/Set the output directory
  Results.dir <- paste0("SDM/SDM_outputs/Results_Set", Set)
  if(!dir.exists(Results.dir)){
    dir.create(Results.dir)
  }
  

  ###################################
  # Model built using AUC: Overview #
  ###################################
  
  # This method uses AUC as a test statistic to evaluate model performance with cross validation
  # Variables are dropped and model hyperparameters are tuned based on AUC
  
  # Occurence records are randomly split into three parts: 
  # training, validation and testing datasets
  # Hold back 20% of the observations to use as a testing dataset to evaluate the final model
  # Use 20% of the data to drive the hyperparameter tuning
  # Use 60% of the observations to train and evaluate the tuned model
  # Train each model using the 10,000 background locations
  # The function trainValTest() allows you to split the data in three folds containing the provided percentage of data.
  c(train, tune, test) %<-% trainValTest(Dh.SWD, test = 0.2, val = 0.2, only_presence = TRUE, seed = 25)
  
  ########################################################################
  # Model built using cross validation and AUC as the performance metric #
  ########################################################################
  
  # Area under the receiver-operator characteristic curve (AUC) is a Threshold-independent metric (TIM). 
  # TIMs generally evaluate the ability of the raw predictions (ranging between 0 to 1) to reflect the propensity of the species to be present or absent at sites.
  
  # AUC has a range from 0 to 1:
  # 0 to 0.5: Model performs worse than random
  # 0.5: Model performs no better or worse than random
  # 0.5 to 1: Model performs better than random
  # In general, AUC values of 0.5-0.7 are considered low and represent poor model performance, values of 0.7-0.9 are considered moderate, and values above 0.9 represent excellent model performance.
  
  
  ####################
  # Cross validation #
  ####################
  
  # If test sites are close to training sites then they probably don't really represent very independent test sites. As a result, model performance will be falsely elevated. Cross validation involves splitting test sites from training sites to use geographic g-folds which divide sites spatially into sections. This helps increase the independence between them and increases the reliability of the evaluation metric.
  
  # The Checkerboard2 method partitions the data into k=4 bins. This is done by aggregating the input raster at two scales. Presence and background points are assigned to a bin with respect to where they fall in checkerboards of both scales.
  
  # Split the training dataset into four random folds (based on a checkerboard pattern) to perform cross validation
  folds.cb <- get.checkerboard2(occs = train@coords[train@pa == 1, ], 
                                envs = Rasters, 
                                bg = train@coords[train@pa == 0, ], 
                                aggregation.factor = c(10, 10))
  
  # Create colour palette
  palette(wes_palette("Rushmore1")[c(2, 4, 3, 5)])
  
  # Check by plotting
  plot(Rasters[[1]], col= "grey", legend=FALSE, axes = FALSE)
  points(train@coords[train@pa == 0, ], pch = 21, bg=folds.cb$bg.grp, cex = 0.2, col = folds.cb$bg.grp)
  points(train@coords[train@pa == 1, ], pch = 24, bg=folds.cb$occ.grp, col='white', cex= 1, lwd = 0.9)
  
  # Save Plot as pdf
  pdf(paste0(Results.dir, "/Checkerboard.CV.pdf"), width = 20/cm(1), height = 15/cm(1), useDingbats = FALSE)
  print(plot(Rasters[[1]], col= "grey", legend=FALSE, axes = FALSE))
  print(points(train@coords[train@pa == 0, ], pch = 21, bg=folds.cb$bg.grp, cex = 0.2, col = folds.cb$bg.grp))
  print(points(train@coords[train@pa == 1, ], pch = 24, bg=folds.cb$occ.grp, col='white', cex= 1, lwd = 0.6))
  dev.off()
  
  
  ###########################
  # Train the default Model #
  ###########################
  
  # Train model with hyperparameters set to default values (will tune later)
  Default.Model.CV.auc <- train(method = "Maxent", data = train, folds = folds.cb)
  
  # Plot ROC curves for each CV fold
  DefaultROC.1 <- plotROC(Default.Model.CV.auc@models[[1]], test = test)
  DefaultROC.2 <- plotROC(Default.Model.CV.auc@models[[2]], test = test)
  DefaultROC.3 <- plotROC(Default.Model.CV.auc@models[[3]], test = test)
  DefaultROC.4 <- plotROC(Default.Model.CV.auc@models[[4]], test = test)
  
  # Put all on same plot
  DefaultROC.AUC <- plot_grid(DefaultROC.1, DefaultROC.2, DefaultROC.3, DefaultROC.4, labels = c("CV Fold 1", "CV Fold 2", "CV Fold 3", "CV Fold 4"), label_size = 10)

  # Save Plot as pdf
  pdf(paste0(Results.dir, "/Default.CV.ROC.pdf"), width = 20/cm(1), height = 16/cm(1), useDingbats = FALSE, pointsize = 6)
  print(DefaultROC.AUC)
  dev.off()
  
  ###############################
  # Remove correlated variables #
  ###############################
  
  # Remove highly correlated variables by repeating the following steps:
  # 1. Rank the variables according to the percent contribution.
  # 2. Check if the variable ranked as most important is highly correlated with other variables, according to a Spearman correlation coefficient > 0.7. If the algorithm finds correlated variables it moves to the next step, otherwise it checks the other variables in the rank.
  # 3. Next, perform a leave one out Jackknife test among the correlated variables;
  # 4. Remove the variable that decreases the model performance the least according to AUC
  # 5. The process is repeated until the remaining variables have a correlation coefficient lower than 0.7
  
  # The function varSel removes the correlated variable with the lowest importance (based on the default model)
  # Variable importance is calculated using the percent contribution computed by Maxent software
  Default.Model.CV.varSel.auc <- varSel(Default.Model.CV.auc, metric = "auc", test = test, bg4cor = Dh.SWD, method = "spearman", cor_th = 0.7, permut = 10, use_pc = TRUE)

  # Plot ROC curves for each CV fold
  Df.Varsel.ROC.1 <- plotROC(Default.Model.CV.varSel.auc@models[[1]], test = test)
  Df.Varsel.ROC.2 <- plotROC(Default.Model.CV.varSel.auc@models[[2]], test = test)
  Df.Varsel.ROC.3 <- plotROC(Default.Model.CV.varSel.auc@models[[3]], test = test)
  Df.Varsel.ROC.4 <- plotROC(Default.Model.CV.varSel.auc@models[[4]], test = test)
  
  # Combine onto the same plot
  Df.Varsel.ROC.AUC <- plot_grid(Df.Varsel.ROC.1, Df.Varsel.ROC.2, Df.Varsel.ROC.3, Df.Varsel.ROC.4, labels = c("CV Fold 1", "CV Fold 2", "CV Fold 3", "CV Fold 4"), label_size = 10)
  
  # Save Plot as pdf
  pdf(paste0(Results.dir, "/Default.CV.VarSel.ROC.pdf"), width = 20/cm(1), height = 16/cm(1), useDingbats = FALSE, pointsize = 6)
  print(Df.Varsel.ROC.AUC)
  dev.off()
  
  ##############################
  # Tune model hyperparameters #
  ##############################
  
  # Maxent includes its own "tuning" procedures to help alleviate over-fitting to data and account for the "correct" amount of complexity in the modeled response. While generally robust, these procedures were refined using a particular (though large) data set. So the model tuning done behind the scenes by Maxent may fit most circumstances fairly well but aren't guaranteed to work well for any particular species. 
  
  # When you tune the model hyperparameters you iteratively adjust the hyperparameters while monitoring the changes in the evaluation metric computed using the testing dataset. In this process, the information contained in the testing dataset leaks in the model and therefore, at the end of the process, the testing dataset no longer represents an independent set to evaluate the model (Muller and Guido 2016). 
  
  # This is why I've split the records into training, validation and testing datasets. The training dataset is used to train the model, the validation datasets to drive the hyperparameter tuning and the testing dataset to evaluate the tuned model.
  
  
  # Fine tune the reduced variable model's hyperparameters using the "tune" dataset 
  # The optimizeModel function checks for the increase in the mean validation AUC across the four cross validation folds. 
  # Search for the best set of hyperparameters among the following values: 
  # 1) feature classes combinations: lq, lh, lqp, lqh, lph, lqph, with l representing linear, q quadratic, p product and h hinge
  # 2) regularization multiplier ranging from 0.2 to 5 with increments of 0.2
  # 3) maximum number of iterations: 300, 500 or 700
  
  # List arguments:
  h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lp", "lqp", "lqph"), iter = c(300, 500, 700))
  
  # Optimise the model
  Model.CV.varSel.Tune.auc <- optimizeModel(Default.Model.CV.varSel.auc, hypers = h, metric = "auc", test = tune, pop = 5, gen = 20, seed = 798)
  
  # save results
  write.csv(Model.CV.varSel.Tune.auc@results, paste0(Results.dir, "/OptimizeMod.Results.AUC.csv"), row.names = FALSE)
  
  # Save best model
  Model.CV.varSel.Tune.best.auc <- Model.CV.varSel.Tune.auc@models[[1]]  # Best model
  
  # Plot ROC curves for each CV fold
  Df.Varsel.Tune.ROC.1 <- plotROC(Model.CV.varSel.Tune.best.auc@models[[1]], test = test)
  Df.Varsel.Tune.ROC.2 <- plotROC(Model.CV.varSel.Tune.best.auc@models[[2]], test = test)
  Df.Varsel.Tune.ROC.3 <- plotROC(Model.CV.varSel.Tune.best.auc@models[[3]], test = test)
  Df.Varsel.Tune.ROC.4 <- plotROC(Model.CV.varSel.Tune.best.auc@models[[4]], test = test)
  
  # Put all on one plot
  Df.Varsel.Tune.ROC <- plot_grid(Df.Varsel.Tune.ROC.1, Df.Varsel.Tune.ROC.2, Df.Varsel.Tune.ROC.3, Df.Varsel.Tune.ROC.4, labels = c("CV Fold 1", "CV Fold 2", "CV Fold 3", "CV Fold 4"), label_size = 10)
  
  # Save Plot as pdf
  pdf(paste0(Results.dir, "/Tuned.CV.VarSel.ROC.pdf"), width = 20/cm(1), height = 16/cm(1), useDingbats = FALSE, pointsize = 6)
  print(Df.Varsel.Tune.ROC)
  dev.off()
  
  
  ########################################
  # Remove variables with low importance #
  ########################################
  
  # Using the tuned model, remove variables with a low permutation importance (<2%) 
  # The reduceVar function uses the Jackknife approach to control for the decrease in the mean validation AUC across the four cross validation folds
  
  # Reduce low importance variables
  Model.CV.varSel.Tune.Best.RedVar.auc <- reduceVar(Model.CV.varSel.Tune.best.auc, th = 2, metric = "auc", test = test, permut = 10, use_jk = TRUE, use_pc = TRUE)
  
  # Plot ROC curves for each CV fold
  Varsel_reduced.Tune.ROC.1 <- plotROC(Model.CV.varSel.Tune.Best.RedVar.auc@models[[1]], test = test)
  Varsel_reduced.Tune.ROC.2 <- plotROC(Model.CV.varSel.Tune.Best.RedVar.auc@models[[2]], test = test)
  Varsel_reduced.Tune.ROC.3 <- plotROC(Model.CV.varSel.Tune.Best.RedVar.auc@models[[3]], test = test)
  Varsel_reduced.Tune.ROC.4 <- plotROC(Model.CV.varSel.Tune.Best.RedVar.auc@models[[4]], test = test)
  
  # Put all on the same grid: 
  Varsel_reduced.Tune.ROC <- plot_grid(Varsel_reduced.Tune.ROC.1, Varsel_reduced.Tune.ROC.2, Varsel_reduced.Tune.ROC.3, Varsel_reduced.Tune.ROC.4, labels = c("CV Fold 1", "CV Fold 2", "CV Fold 3", "CV Fold 4"), label_size = 10)

  # Save Plot as pdf
  pdf(paste0(Results.dir, "/Tuned.CV.VarSel_reduced.ROC.pdf"), width = 20/cm(1), height = 16/cm(1), useDingbats = FALSE, pointsize = 6)
  print(Varsel_reduced.Tune.ROC)
  dev.off()
  
  
  ########################
  # Evaluate final model #
  ########################
  
  # Use the best tuned model, with the optimal set of reduced variables to create/evaluate the final model. 
  # Before evaluating this model, merge the training (and the cross validation folds within) and the tuning datasets together to increase the number of locations.
  # Train a new model with the merged observations and the tuned configuration. 
  # Variables have been removed, so cannot directly merge the original datasets containing all environmental variables. 
  # Instead, extract the train dataset with the selected variables from the output of the best tuned model and merge it with the validation dataset using the function mergeSWD():
  new_train.AUC <- Model.CV.varSel.Tune.Best.RedVar.auc@data
  merged_data.AUC <- mergeSWD(new_train.AUC, tune, only_presence = TRUE) 
  # NOTE: Merge presence data only
  
  # Train final model
  final_model.AUC <- train("Maxent", data = merged_data.AUC, reg = Model.CV.varSel.Tune.best.auc@models[[1]]@model@reg, fc = Model.CV.varSel.Tune.best.auc@models[[1]]@model@fc, iter = Model.CV.varSel.Tune.best.auc@models[[1]]@model@iter)

  # Save Plot as pdf
  pdf(paste0(Results.dir, "/FinalModel.AUC.ROC.pdf"), width = 16/cm(1), height = 16/cm(1), useDingbats = FALSE, pointsize = 6)
  print(plotROC(final_model.AUC, test = test))
  dev.off()
  
  ###################
  # Save all models #
  ###################
  
  save(object = Default.Model.CV.auc, file = paste0(Results.dir, "/Default.Model.CV.auc.rdata"))
  save(object = Default.Model.CV.varSel.auc, file = paste0(Results.dir, "/Default.Model.CV.varSel.auc.rdata"))
  save(object = Model.CV.varSel.Tune.auc, file = paste0(Results.dir, "/Model.CV.varSel.Tune.auc.rdata"))
  save(object = Model.CV.varSel.Tune.best.auc, file = paste0(Results.dir, "/Model.CV.varSel.Tune.best.auc.rdata"))
  save(object = Model.CV.varSel.Tune.Best.RedVar.auc, file = paste0(Results.dir, "/Model.CV.varSel.Tune.Best.RedVar.auc.rdata"))
  save(object = final_model.AUC, file = paste0(Results.dir, "/final_model.AUC.rdata"))
  
  
  #######################
  # Variable Importance #
  #######################
  
  # Maxent models provide the variable importance values in the output 
  # These values are stored in the model object and can be displayed using the maxentVarImp() command
  # this function extracts the variable importance values and formats them in a more human readable way
  VI.AUC <- maxentVarImp(final_model.AUC)

  # Save Plot as pdf
  pdf(paste0(Results.dir, "/FinalModel.VarImp.AUC.pdf"), width = 16/cm(1), height = 16/cm(1), useDingbats = FALSE, pointsize = 6)
  print(plotVarImp(VI.AUC))
  dev.off()
  
  
  #############################
  # Create a distribution map #
  #############################
  
  # Plot the final model by projecting to the full extent of the study area, applying the cloglog transformation (Phillips et al., 2017) to the raw output of the model
  map.AUC <- predict(final_model.AUC, data = Rasters, type = "cloglog")
  
  # Make a nice plot for publication
  # For ggplot, have to convert raster to data frame
  map.AUC_df <- as.data.frame(map.AUC, xy = TRUE) 
  
  # Remove scientific notation from UTMs in plot
  options(scipen = 999)
  
  # Create a colour palette
  HabPal <- c("black", parula(n = 20))
  
  Map.AUC.plot <- ggplot() +
    geom_raster(data = map.AUC_df, 
                aes(x = x, y = y, fill = layer)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(210000, 1030000), ylim = c(7305000, 7792000)) +
    scale_fill_gradientn(colours = HabPal, na.value = "white") +
    labs(x = "UTMX", y = "UTMY", fill = "Habitat\nSuitability\n") +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", 
                                      fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 8, color = "black"),
          axis.title.y = element_text(size = 9, 
                                      margin = margin(0,10,0,0), 
                                      face ="bold"),
          axis.title.x = element_text(size = 9, 
                                      margin = margin(10,0,0,0), 
                                      face ="bold"),
          axis.ticks = element_line(size = 0.2),
          legend.text = element_text(size = 8, color = "black"),
          legend.title = element_text(size = 9, 
                                      margin = margin(10,0,0,0), 
                                      face ="bold"),
          legend.position = "right") + 
    scalebar(location = "bottomleft", dist_unit = "km", 
             transform = FALSE,
             x.min =  250000, x.max = 1000000,
             y.min = 7330000, y.max = 7782000,
             dist = 100, st.dist = 0.025, st.size = 2.5,
             height = 0.02, border.size = 0.15) +
    north(location = "topleft", scale = 0.05, symbol = 12,
          x.min = 220000, x.max = 1000000, 
          y.min = 7330000, y.max = 7782000)
  
 
  # Save as pdf
  ggsave(filename = "FinalModel_MapPredict.AUC.pdf",
         plot = Map.AUC.plot,
         device = "pdf",
         path = Results.dir,
         scale = 1,
         width = 22,
         height = 12,
         units = "cm",
         dpi = 300)
  
  # Save as raster
  writeRaster(map.AUC, 
              paste0(Results.dir, "/final_model.AUC_SDM.tif"), 
              format='GTiff', 
              datatype= dataType(map.AUC))
  
  # Create another map that also has presence points
  Map.AUC.samps.plot <- Map.AUC.plot + geom_point(data = p_coords, aes(x = UTMX, y = UTMY), 
                                                  size = 1.5, stroke = 0.2, shape = 24, 
                                                  fill = "white", colour = "black")
  
  
  # Save as pdf
  ggsave(filename = "FinalModel_MapPredict.AUC_Samps.pdf",
         plot = Map.AUC.samps.plot,
         device = "pdf",
         path = Results.dir,
         scale = 1,
         width = 22,
         height = 12,
         units = "cm",
         dpi = 300)
  
  
  
  ################################
  # Plot species response curves #
  ################################
  
  # With the function plotResponse() it is possible to plot the cloglog marginal and the univariate response curve.
  # On top is displayed the rug of the presence locations and on bottom the rug of the background locations. 
  
  # In the case of an SDMmodelCV() the response curve shows the averaged value of the prediction together with one Standard Deviation error interval. We can also use the cross validation model (where we have standard error across folds)
  
  # Use lapply to loop through each variable and plot response curves for both the best, tuned cross validation model and the final model
  # Note that the plotResponse function creates a ggplot object, so can edit this using ggplot themes etc.
  
  # Define variables
  vars.AUC <- Rasters.AllNames[Rasters.AllNames %in% colnames(final_model.AUC@data@data)]
  var.AUC.names <- AllNames[Rasters.AllNames %in% colnames(final_model.AUC@data@data)]
  
  #### MARGINAL RESPONSE CURVES ####
  
  # Cross validation model
  CV.resp.marg = lapply(1:length(var.AUC.names), function(X)
    plotResponse(Model.CV.varSel.Tune.Best.RedVar.auc, var = vars.AUC[X], type = "cloglog", only_presence = TRUE, marginal = TRUE, fun = mean, rug = TRUE, color = "red") + 
      labs(x = var.AUC.names[X]) +
      scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.25)) +
      theme(panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(linetype = "solid", colour = "black", 
                                        fill = "transparent", size = 0.2),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            axis.ticks = element_line(size = 0.2),
            axis.text = element_text(size = 8, color = "black"),
            axis.title = element_text(size = 9, 
                                      face ="bold")))
  
  # Final model
  # Note, I'm going to overlay this onto the plots above, so I've made most of the plot transparent (e.g. borders, ticks, text etc.) except for the curve
  AUC.resp.marg = lapply(1:length(var.AUC.names), function(X)
    plotResponse(final_model.AUC, var = vars.AUC[X], type = "cloglog", only_presence = TRUE, marginal = TRUE, fun = mean, rug = FALSE, color = "transparent") + 
      labs(x = "", y = "") + 
      geom_line(linetype = "dashed", size = 0.3, color = "darkblue") +
      scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.25)) +
      theme(panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(linetype = "solid", colour = "transparent", 
                                        fill = "transparent", size = 0.2),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            axis.ticks = element_line(size = 0.2, color = "transparent"),
            axis.text = element_text(size = 8, color = "transparent"),
            axis.title = element_text(size = 9, color = "transparent",
                                      face ="bold")))
  
  # Align the cross validation response curve and the final model curve for each variable
  resp.plots.marg = lapply(1:length(var.AUC.names), function(X)
    ggdraw(align_plots(AUC.resp.marg[[X]], CV.resp.marg[[X]], align="v")[[1]]) + 
      draw_plot(align_plots(AUC.resp.marg[[X]], CV.resp.marg[[X]], align="v")[[2]]))
  
  # Put all on the same plot
  resp.grid.marg.AUC <- plot_grid(plotlist = resp.plots.marg, nrow = ceiling(length(var.AUC.names)/2), ncol = 2)
  
  
  # Save as pdf
  ggsave(filename = "Response_AUC_Marg.final.bestCV.pdf",
         plot = resp.grid.marg.AUC,
         device = "pdf",
         path = Results.dir,
         scale = 1,
         width = 20,
         height = 30,
         units = "cm",
         dpi = 300)
  
  
  #### UNIVARIATE RESPONSE CURVES ####
  
  # The univariate option reruns the maxent model using only the variable of interest
  
  # Cross validation model
  CV.resp.uni = lapply(1:length(var.AUC.names), function(X)
    plotResponse(Model.CV.varSel.Tune.Best.RedVar.auc, var = vars.AUC[X], type = "cloglog", only_presence = TRUE, marginal = FALSE, fun = mean, rug = TRUE, color = "red") + 
      labs(x = var.AUC.names[X]) +
      scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.25)) +
      theme(panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(linetype = "solid", colour = "black", 
                                        fill = "transparent", size = 0.2),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            axis.ticks = element_line(size = 0.2),
            axis.text = element_text(size = 8, color = "black"),
            axis.title = element_text(size = 9, 
                                      face ="bold")))
  
  # Final model
  # Note, I'm going to overlay this onto the plots above, so I've made most of the plot transparent (e.g. borders, ticks, text etc.) except for the curve
  AUC.resp.uni = lapply(1:length(var.AUC.names), function(X)
    plotResponse(final_model.AUC, var = vars.AUC[X], type = "cloglog", only_presence = TRUE, marginal = FALSE, fun = mean, rug = FALSE, color = "transparent") + 
      labs(x = "", y = "") + 
      geom_line(linetype = "dashed", size = 0.3, color = "darkblue") +
      scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.25)) +
      theme(panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(linetype = "solid", colour = "transparent", 
                                        fill = "transparent", size = 0.2),
            panel.grid.minor = element_blank(), 
            panel.grid.major = element_blank(),
            axis.ticks = element_line(size = 0.2, color = "transparent"),
            axis.text = element_text(size = 8, color = "transparent"),
            axis.title = element_text(size = 9, color = "transparent",
                                      face ="bold")))
  
  # Align the cross validation response curve and the final model curve for each variable
  resp.plots.uni = lapply(1:length(var.AUC.names), function(X)
    ggdraw(align_plots(AUC.resp.uni[[X]], CV.resp.uni[[X]], align="v")[[1]]) + 
      draw_plot(align_plots(AUC.resp.uni[[X]], CV.resp.uni[[X]], align="v")[[2]]))
  
  # Put all on the same plot
  resp.grid.uni.AUC <- plot_grid(plotlist = resp.plots.uni, nrow = ceiling(length(var.AUC.names)/2), ncol = 2)
  
  
  # Save as pdf
  ggsave(filename = "Response_AUC_Uni.final.bestCV.pdf",
         plot = resp.grid.uni.AUC,
         device = "pdf",
         path = Results.dir,
         scale = 1,
         width = 20,
         height = 30,
         units = "cm",
         dpi = 300)

}
