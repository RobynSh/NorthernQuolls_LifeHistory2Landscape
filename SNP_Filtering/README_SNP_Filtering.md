## README

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

Repository containing R code (*DArTSNPFilt.R*), tutorial, data, and outputs for DArTSeq (https://www.diversityarrays.com/) SNP visualisation and filtering. 

This repository steps through the process of filtering DArTSeq SNPs (although the steps are relevant to any reduced representation SNP data) using the following files in the "Data" folder as the starting input:  

### Data supplied by DArT:
* 1-row SNP csv (*Report_DDasy19-4717_1Row_NameEdit.csv*)  
* 2-row SNP csv (*Report_DDasy19-4717_2Row_NameEdit.csv*)  
* Read count csv (*ReadCounts_DDasy19-4717_NameEdit.csv*)  

### Study specific data:  
* Sample meta-data csv (including ID, pop, lon/lat, other info; *Dasyurus_hallucatus_ind.metadata.csv*)  
* Shape file of study region (*PilbaraIBRA.shp*)  
* Optional: Tasmanian devil chromosome info (*TasDevil7.1_Scaffold_Info.csv*)  

Filtering outputs (plots, filtered data) can be found in the "Filtering_outputs" folder. The full process is outlined in detail in the html tutorial (*DArTSNPFilt.html*).