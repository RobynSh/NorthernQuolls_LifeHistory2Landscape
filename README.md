<p align="center">
<b><i>Linking life history to landscape for threatened species conservation in a multi-use region</i></b>
</p>
<div align="center">
    <img src="Figs/Paper_Figure1.png" width="700px"</img> 
</div>


## **Complete code for Conservation Biology paper**
### *doi: 10.1111/cobi.13989*

This repository contains R code, data and outputs for all analyses presented in the paper. The analyses are split into the following folders to reflect the main sections/subheadings in the paper (presented in order of analysis - note that outputs generated are sometimes used as inputs in the next step/s, so the order of analysis matters):
* **SDM:**
  - *Step1_RecordCleaning.R*
  - *Step2_RunMaxent.R*
* **SNP_Filtering:**
  - *Step3_DArTSNPFilt.R*
  - *Step4_SampleCleaning.R*
* **IBB:**
  - *Step5_ClusteringAnalyses.R*
  - *Step6_Genalex.GENHET.R*
* **IBR:**
  - *Step7_IBR.R*
* **IBE:**
  - *Step8_IBE.R*

There are README files in each of the analysis sub-folders listed above, which describe the analysis, the data (inputs) and the resulting output files. 

The remaining sub folders contain data used throughout the different analyses, as described below:
* **Figs:** This folder contains figures from the paper that are presented in the README files
* **Rasters_Shapefiles:** This folder contains rasters and shapefiles covering the Pilbara region of Western Australia. All rasters in this folder have been aggregated to a 1 km resolution. For details on spatial layers, see the Supplementary Information for this paper (Appendix S2), and the figure below.
* **Rasters_5km:** This folder contains a subset of the rasters described above, that have been further aggregated to a 5 km resolution (for Isolation-By-Resistance modelling). More detials are in the README file under **IBR**.


<p align="center">
<b><i>Spatial layers explored in this study with associated acronyms, the hypothesised mechanism driving spatial-environmental associations, and how these link to the life cycle/demography of the northern quoll (numbers match above figure).</i></b>
</p>
<div align="center">
    <img src="Figs/Paper_Figure2.png" width="700px"</img> 
</div>

##### The specific hypotheses include HS= habitat suitability, tested using Species Distribution Modelling; IBB= isolation-by-barrier, environmental associations tested using Generalised Dissimilarity Modelling (GDM); IBR= isolation-by-resistance, tested using linear mixed effects models with MLPE; IBE= isolation-by-environment, tested using GDM.



&nbsp;

&nbsp;

<div align="center">
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
</div>
