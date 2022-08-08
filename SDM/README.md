<div align="center">
    <img src="Data/Quoll.png" width="250px"</img> 
</div>
<p align="center">
<b><i>Dasyurus hallucatus</i></b>
</p>

# Step 1: Record Cleaning 

The R code *Step1_RecordCleaning.R* documents the process of cleaning Atlas of Living Australia and NatureMap records for:
* *D. hallucatus* - for modelling its distribution (input file, *Dh.ALA_NM_Raw.Pilbara.Records.csv*, is found in the Data folder)
* All Critical Weight Range mammals - for generating a survey bias layer (input files, *Mammal.spp_NM_Raw.Records.csv* and *MammalSpecies_Info.csv*, is found in the Data folder)

This code generates the following files (output to the Data folder):
* *Dh_FinalPresenceCoords.csv* - the final cleaned record coordinates, i.e. presences for MaxEnt modelling in step 2
* *CWR_Bias_Layer.asc* - a raster 

Prior to this step, we merged ALA and NatureMap data sets and removed irrelevant columns and duplicated records (which left only one ALA record, all of the others were already represented in the NatureMap data base). Missing coordinate info was filtered out prior to this script (as this was used to find duplicate records across data bases). Rasters were also derived/prepared/manipulated prior to this step. For example, all rasters have been aggregated to the same resolution based on the most biologically appropriate metric (mean, median, min, max, etc.), and cropped/masked to the same extent (the Pilbara IBRA region, with some coverage into surrounding regions). For details on raster manipulation/creation see the Supplementary material for the paper (Appendix S2-S3).





# Step 2: Species Distribution Modelling




Data:
CWR_Bias_Layer.asc
Dh_FinalBackgroundCoords.csv
Dh_FinalPresenceCoords.csv
Dh.ALA_NM_Clean.Pilbara.Records.csv
Dh.ALA_NM_Raw.Pilbara.Records.csv
Mammal.spp_NM_Raw.Records.csv
MammalSpecies_Info.csv"  

SDM_outputs:
RasterCorr_Total1km_0.7Only.pdf
Results_Set1
Results_Set2
Results_Set3
Results_Set4
Results_Set5
Results_Set6

Checkerboard.CV.pdf
Default.CV.ROC.pdf
Default.CV.VarSel.ROC.pdf
Default.Model.CV.auc.rdata
Default.Model.CV.varSel.auc.rdata
final_model.AUC_SDM.tif
final_model.AUC.rdata
FinalModel_MapPredict.AUC_Samps.pdf
FinalModel_MapPredict.AUC.pdf
FinalModel.AUC.ROC.pdf
FinalModel.VarImp.AUC.pdf
Model.CV.varSel.Tune.auc.rdata
Model.CV.varSel.Tune.best.auc.rdata
Model.CV.varSel.Tune.Best.RedVar.auc.rdata
OptimizeMod.Results.AUC.csv
Response_AUC_Marg.final.bestCV.pdf
Response_AUC_Uni.final.bestCV.pdf
Tuned.CV.VarSel_reduced.ROC.pdf
Tuned.CV.VarSel.ROC.pdf


&nbsp;
<div align="center">
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
</div>
