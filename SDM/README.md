<p align="center">
  <img src="https://zenodo.org/badge/520347690.svg">
</p>

<p align="center">
<b><i>Species Distribution Modelling: Generating a habitat suitability map for D. hallucatus</i></b>
</p>
<div align="center">
    <img src="../Figs/Paper_Figure3.png" width="700px"</img> 
</div>


# Step 1: Record Cleaning 

The R code *Step1_RecordCleaning.R* documents the process of cleaning Atlas of Living Australia and NatureMap records for:
* *D. hallucatus* - for modelling its distribution 
* All Western Australian mammals - for generating a survey bias layer using those in the Pilbara within the Critical Weight Range (CWR)

### Input files:
* *D. hallucatus* raw occurrence records - *Data/Dh.ALA_NM_Raw.Pilbara.Records.csv*
* All Western Australian mammal raw occurrence records - *Data/Mammal.spp_NM_Raw.Records.csv*
* Information on CWR mammals for creating bias layer - *Data/MammalSpecies_Info.csv*
* Shape file of study region (*PilbaraIBRA7_UTM.shp*) and rasters representing the Pilbara environment and landscape characteristics in the *Rasters_Shapefiles* folder in the home directory (*.asc* files, specific files documented in code).

### Output files:
This code generates the following files (output to the Data folder):
* *Dh_FinalPresenceCoords.csv* - the final cleaned record coordinates, i.e. presences for MaxEnt modelling in step 2.
* *CWR_Bias_Layer.asc* - a probability function reflecting survey effort represented spatially as a raster. This is used in step 2 for generating background points (pseudo-absences).
* *RasterCorr_Total1km_0.7Only.pdf* - plot displaying raster correlations (only those with at least one pairwise comparison greater than/less than 0.7/-0.7 are shown).

Prior to this step, we merged ALA and NatureMap data sets and removed unnecessary columns and duplicated records (which left only one ALA record, all of the others were already represented in the NatureMap data base). Records with missing coordinates were filtered out prior to this script. Rasters were also derived/prepared/manipulated prior to this step. For example, all rasters have been aggregated to the same resolution based on the most biologically appropriate metric (mean, median, min, max, etc.), and cropped/masked to the same extent (the Pilbara IBRA region, with some coverage into surrounding regions). For details on raster manipulation/creation see the Supplementary material for the paper (Appendix S2-S3).



# Step 2: Species Distribution Modelling

The R code *Step2_RunMaxEnt.R* documents the process of:
- Loading in rasters, defining categorical/continuous variables, and splitting them into different sets, i.e. different representations of the same environmental variables (for example, categorical vs. continuous soil variables).
- Generating background points
- Running MaxEnt with the package SDMtune (looping through for each raster set), by first splitting data for training (with cross validation), tuning and testing. Next, we run the default model, drop correlated variables, tune hyperparameters, drop variables with low importance and run the final model using optimised hyperparameters and variables from the previous steps. 

### Input files:
* Shape file of study region (*PilbaraIBRA7_UTM.shp*) and rasters representing the Pilbara environment and landscape characteristics in the *Rasters_Shapefiles* folder in the home directory (*.asc* files, specific files documented in code).
* Bias layer created in step 1 to generate background points, where more points are taken from heavily surveyed areas - *Data/CWR_Bias_Layer.asc*
* In this script, background points are loaded back in. Unfortunately, I didn't set the seed when I ran this script originally. Therefore, to make the study repeatable, I have provided the data for the background points used to generate the SDM presented in the paper. The process used to obtain these points was identical to that documented in the code - *Data/Dh_FinalBackgroundCoords.csv* 

### Output files:
Throughout MaxEnt modelling, figures/tables are output into the folder SDM_outputs, in the relevant subfolders for each ratser set:
 - Results_Set1
 - Results_Set2
 - Results_Set3
 - Results_Set4
 - Results_Set5
 - Results_Set6
 
Figures/tables/data outputs: 
  - *Checkerboard.CV.pdf* - plot displaying geographically separated cross validation folds
  - *Default.CV.ROC.pdf* - ROC Curves for each CV fold for the default model
  - *Default.CV.VarSel.ROC.pdf* - ROC Curves for each CV fold for the model after correlated variables dropped
  - *Default.Model.CV.auc.rdata* - Default MaxEnt model R object
  - *Default.Model.CV.varSel.auc.rdata* - MaxEnt model R object after correlated variables dropped
  - *final_model.AUC_SDM.tif* - Habitat suitability raster generated from the final model
  - *final_model.AUC.rdata* - The final model R object
  - *FinalModel_MapPredict.AUC_Samps.pdf* - Figure showing the habitat suitability predictions from the final model (with presence locations)
  - *FinalModel_MapPredict.AUC.pdf*  - Figure showing the habitat suitability predictions from the final model (without presence locations)
  - *FinalModel.AUC.ROC.pdf* - ROC Curve for the final model
  - *FinalModel.VarImp.AUC.pdf* - Variable importance bar chart for the final model
  - *Model.CV.varSel.Tune.auc.rdata* - R object for hyperparameter tuning
  - *Model.CV.varSel.Tune.best.auc.rdata*- R object for the model after correlated variables dropped, with best/tuned hyperparameters
  - *Model.CV.varSel.Tune.Best.RedVar.auc.rdata*- R object for the model after correlated variables dropped, hyperparameter tuning, and low importance variables dropped
  - *OptimizeMod.Results.AUC.csv* - Hyperparameter tuning results saved in table format
  - *Response_AUC_Marg.final.bestCV.pdf* - Marginal response curves for the final model
  - *Response_AUC_Uni.final.bestCV.pdf* - Univariate response curves for the final model
  - *Tuned.CV.VarSel_reduced.ROC.pdf* - ROC Curves for each CV fold for the model after correlated variables dropped, hyperparameters tuned and low importance variables dropped 
  - *Tuned.CV.VarSel.ROC.pdf* - ROC Curves for each CV fold for the model after correlated variables dropped and hyperparameters tuned


**Note that the best performing model presented in the paper is from raster set 3.**


&nbsp;

&nbsp;

<div align="center">
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
</div>
