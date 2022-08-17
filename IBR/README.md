<p align="center">
<b><i>Isolation-By-Environment: Generalised Dissimilarity Modelling for D. hallucatus</i></b>
</p>
<div align="center">
    <img src="IBE_outputs/Paper_Figure6.png" width="700px"</img> 
</div>


# Step 8: Isolation-By-Environment 

The R code *Step8_IBE.R* documents the process of:
* Calculating mean pairwise Bray Curtis genetic distance within a 1km area (to match the pixel size of the rasters used in this analysis)
* Creating the environmental predictor table from the relevant rasters
* Running univariate GDMs to aid in variable selection - i.e. any that explain less than 5% of the variance are dropped at this stage
* Dropping correlated variables using Spearman's correlation and the Variance Inflation Factor (by keeping the one that explained the most variance)
* Running a multivariate GDM on the uncorrelated variable set
* Calculating variable importance and using backwards elimination to remove any variables that are not significant in the model (to avoid over fitting)
* Plotting spatially interpolated genetic gradients using a PCA of the GDM-transformed environmental predictors, assigning PCs to an RGB palette
* Carrying out variation partitioning, by running uni and multvariate GDMs with all combinations of variables to determine the unique and combined contributions of each variable to the final model


### Input files:
* *D. hallucatus* genlight object (*Cleaned.Unrelated.gl.Dh.rdata*) after SNP filtering and sample cleaning and associated metadata (*Cleaned.Unrelated.Ind.metadata.Dh.csv*), imported from the *SNP_Filtering/SampleClean_outputs* folder.
* Shape file of study region (all Islands removed except for Dolphin Island: *IBRA7_Mainland.and.DolphinIslandOnly.shp*) and rasters (1 km resolution) to test IBE hypotheses in the *Rasters_Shapefiles* folder in the home directory (*.asc* files, specific files documented in code).


### Output files:
This code generates the following files (output to the *IBE_outputs* folder):
* *Dh_GDM_Ispline_1000bs.pdf* - Figure: Ispline functions representing partial contributions of individual variables to variation in genetic distance in the final model, with error bars representing 1000 bootstraps
* *Dh_GDM_Ispline.pdf* - Figure: Relationship between environmental dissimilarity and genetic dissimilarity from the final model
* *Dh_GDM_VarImp.pdf* - Figure: Relative contribution of geographic distance and environmental variables to the final multivariate GDM
* *GDM.PCA.Map.pdf* - Figure: Spatial interpolation of genetic gradients using RGB colour composites derived from a principal components analysis of GDM-transformed environmental predictors
* *Venn.PercVar.GDM.pdf* - Figure: the unique and combined contributions of each variable to the explained model deviance (variation partitioning). Note that B15 was not included due to the difficulty of interpreting a four-way Venn diagram.

*Note that the final GDM figure in the paper was combined/edited in Adobe Illustrator*


&nbsp;

&nbsp;

<div align="center">
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
</div>
