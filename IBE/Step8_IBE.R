#############################################
#                 STEP 8                    #  
#       Isolation-By-Environment (IBE)       #
#         Data set: Dasyurus hallucatus     #
#            Author: Robyn Shaw             #
#             Date: 23/02/2022              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

#############################
# LOAD IN REQUIRED PACKAGES #
#############################

library(raster)
library(ecodist)
library(reshape2)
library(dplyr)
library(splitstackshape)
library(gdm)
library(GGally)
library(legendary)
library(usdm)
library(rgdal)
library(RStoolbox)
library(ggsn)
library(gtools)
library(tidyr)
library(stringr)
library(Vennerable)


###################################################
# Create Bray Curtis genetic dissimilarity matrix #
###################################################

# Note that because the raster pixel size is 1km, I want to get the mean genetic dissimilarity of individuals within this 1km area
# I am using Bray Curtis distance as my measure of dissimilarity

# Load in genetic data
load("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.gl.Dh.rdata")
gl <- Dh.gl

# Load in meta-data
samps <- read.csv("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.Ind.metadata.Dh.csv")

# Add meta-data to genlight
gl@other$ind.metrics <- samps[match(gl@ind.names, samps$id), ]

# Read in one of the processed rasters so that I can aggregate the genetic data based on the pixel size of the rasters 
Raster <- raster("Rasters_Shapefiles/ADI.asc")

# Get xy coordinates of raster pixels that overlap with quoll samples
cell <- data.frame(id = gl@ind.names, cell = raster::extract(x = Raster, y = gl@other$ind.metrics[, c("UTMX", "UTMY")], cellnumbers=T)[, 1], Dat = raster::extract(x = Raster, y = gl@other$ind.metrics[, c("UTMX", "UTMY")], cellnumbers=T)[, 2])

# Find any samples outside of raster extent (to remove later)
rmDh <- cell$id[is.na(cell$Dat)]
xy <- cbind(cell, xyFromCell(Raster, cell[,2]))

# How many unique points?
length(unique(xy$cell)) #126

# Create a new column called "pixel" to identify duplicated coords easily (i.e. all coordinates that fall within the same 1km pixel)
# Create a new data frame with only the unique coordinates
pixel.group <- data.frame(cell = unique(xy$cell), pixel = c(paste0("P00", 1:9), paste0("P0", 10:99), paste0("P", 100:length(unique(xy$cell)))))

# Rename individuals in genlight with pixel id
gl@ind.names <- as.character(pixel.group[match(xy$cell, pixel.group$cell), "pixel"])

# Remove samples outside of raster extent
gl <- gl[!gl@other$ind.metrics$id %in% rmDh, ]

# Order by ind.name
gl <- gl[order(gl@ind.names), ]

# Create a distance matrix
# Some of the genetic distance methods can't deal with missing data
# I'll interpolate missing data by using the mean
gl.mat <- round(tab(gl, NA.method="mean"), 0)

# Calculate Bray Curtis distance between genetic samples
BCGenDist <- distance(gl.mat, method = "bray-curtis")

# Convert lower triangular to table
BCGenDist.mat <- as.matrix(BCGenDist)
BCGenDist.tab <- melt(BCGenDist.mat)[melt(lower.tri(BCGenDist.mat))$value, ]
names(BCGenDist.tab) <- c("Ind2", "Ind1", "BCDist")

# Remove comparisons between same pixel group
BCGenDist.tab <- BCGenDist.tab[!BCGenDist.tab$Ind1 == BCGenDist.tab$Ind2,]

# Create a new column to group by (and take the mean of)
BC.pixdf <- data.frame(Ind2 = as.numeric(sub(".*P", "", BCGenDist.tab$Ind2)), Ind1 = as.numeric(sub(".*P", "", BCGenDist.tab$Ind1)))

# Make a new column with both pixel groups so I can get the mean across duplicates
BCGenDist.tab$Pair <- paste(ifelse(apply(BC.pixdf, 1, which.min) == 1, as.character(BCGenDist.tab$Ind2), as.character(BCGenDist.tab$Ind1)), ifelse(apply(BC.pixdf, 1, which.max) == 1, as.character(BCGenDist.tab$Ind2), as.character(BCGenDist.tab$Ind1)), sep = "_")

# Group by pair and then take the mean
BCGenDist.Av <- BCGenDist.tab %>%
  group_by(Pair) %>%
  summarise(BCDistMean = mean(BCDist))

# Check order is correct
BCGenDist.Av$Pair

# Check correct number of values
(length(unique(gl@ind.names))*length(unique(gl@ind.names))-length(unique(gl@ind.names)))/2 == nrow(BCGenDist.Av)

# Add back in same pair comparisons (for the diagonal zero value in the matrix)
BCGenDist.Av <- rbind(BCGenDist.Av, data.frame(Pair = paste0(unique(gl@ind.names), "_", unique(gl@ind.names)), BCDistMean = 0))

# Split back into ind1, ind2 columns
BCGenDistMean.tab <- cbind(cSplit(BCGenDist.Av[,1], 1, sep="_"), BCGenDist.Av$BCDistMean)

# Order
BCGenDistMean.tab <- BCGenDistMean.tab[order(BCGenDistMean.tab$Pair_1, BCGenDistMean.tab$Pair_2), ]

# Create the final genetic distance matrix
GenDistMat.BC <- as.dist(xtabs(V2~., BCGenDistMean.tab) + t(xtabs(V2~., BCGenDistMean.tab)), upper = TRUE, diag = TRUE)

# Write to file
GenDistMat.BC <- as.matrix(GenDistMat.BC)
GenDistMat.BC <- cbind(IDNo = rownames(GenDistMat.BC), GenDistMat.BC)
rownames(GenDist.BC) <- 1:nrow(GenDistMat.BC)

########################################
# Create environmental variable matrix #
########################################

# Get the coordinates for the pixel the sample groups fall in
xy.pix <- xy[!xy$id %in% rmDh, c("id", "x", "y")]
xy.pix <- cbind(gl@ind.names[match(xy.pix$id, gl@other$ind.metrics$id)], xy.pix[, c("x", "y")])
colnames(xy.pix)[1] <- "IDNo"
xy.pix <- unique(xy.pix)

# Save as spatial points object
Dh.sp <- SpatialPointsDataFrame(xy.pix[, c("x", "y")], proj4string = CRS("+proj=utm +zone=50 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"), data = xy.pix)
plot(Dh.sp)


##############################
# Generate raster predictors #
##############################

# Load in rasters
All.rasters <- stack(list.files("Rasters_Shapefiles/", pattern = ".asc", full.names = TRUE))

# Remove irrelevant layers
All.rasters <- All.rasters[[-which(names(All.rasters) %in% c("Clay", "CoarseFrags", "DEM", "Dist2Water", "PersistForestCover", "Sand", "Silt", "Soils_Preferred", "Soils_Resistance", "Spinifex", "WeatheringIndex"))]]



######### HYPOTHESES #########
# Does climate drive genetic dissimilarity/turnover?
## Test aridity indices, bioclim variables and bioclim PC axes

# Are rocky/rugged versus savanna quolls genetically different
## Test Rocky soils (rock vs. not rock), VRM (rugged vs. not rugged, and continuous layer)

# Test whether quolls are different north vs. south of the fortescue (and if other major rivers act as barriers)
## Test Fortescue and MajorRivers layers

# Are quolls near water (i.e. in productive areas) genetically different from those further from water?
## Test Dist2Water layers

# Are quolls in frequently burnt habitat genetically different?
# Test fire frequency layer


#################################
# Create the "Predictors" table #
#################################

# Match individual coordinates with environment at each location
# i.e. sample the raster pixel for each individual's coordinate
Dh.Env <- cbind(data.frame(UTMX = Dh.sp@coords[,1], UTMY = Dh.sp@coords[,2], IDNo = xy.pix$IDNo), as.data.frame(raster::extract(x = All.rasters, y = cbind(Dh.sp@coords))))

# Remove NA samples
any(is.na(rowSums(Dh.Env[4:ncol(Dh.Env)])))

#######################
# RUN UNIVARIATE GDMS #
#######################

# GDM is a nonlinear extension of permutational matrix regression that uses flexible splines and a GLM to accommodate two types of nonlinearity common in ecological datasets: (1) variation in the rate of compositional turnover (non-stationarity) along environmental gradients, and (2) the curvilinear relationship between biological distance and environmental and geographical distance.

# This loop: 
# 1. Combines the biological and environmental data into site-pair format. 
# 2. Fits generalized dissimilarity model

# Exclude variable if explain <5% variation (when geo = F) (Supple paper uses 5%)

CatRast <- c("Fortescue", "MajorRivers", "RiverBasins", "RockySoils", "Dist2Water.1km", "Dist2Water.5km", "VRM.10", "VRM.20", "VRM.5")
ContRast <- names(All.rasters)[-which(names(All.rasters) %in% CatRast)]

Include <- vector()
VarExp <- vector()

for (i in 1:length(names(All.rasters))) {
  envTab = Dh.Env[, c(1:3, i+3)]
  if (colnames(envTab)[4] %in% ContRast) {
  gdm.dis <- formatsitepair(GenDist.BC, 
                            bioFormat=3, 
                            XColumn="UTMX", 
                            YColumn="UTMY",
                            predData=envTab, 
                            siteColumn="IDNo")
  gdm <- gdm(gdm.dis, geo = FALSE)
  
  # Need to format differently if categorical variable
  } else if (colnames(envTab)[4] %in% CatRast) {
    gdm.dis <- formatsitepair(GenDist.BC, 
                              bioFormat=3, 
                              XColumn="UTMX", 
                              YColumn="UTMY",
                              predData=envTab, 
                              siteColumn="IDNo")
    
    gdm.dis$s1 <- 0
    gdm.dis$s2 <- ifelse(gdm.dis[,7] == gdm.dis[,8], 0, 1)
    
    Cols <- colnames(gdm.dis)[7:8]
    
    gdm.dis <- gdm.dis[, -c(7:8)]
    colnames(gdm.dis)[7:8] <- Cols
    
    gdm <- gdm(gdm.dis, geo = FALSE)
  }
  
  if (is.null(gdm)) {
  } else if (gdm$explained >= 5) {
    Include <- c(Include, colnames(envTab)[ncol(envTab)])
    VarExp <- c(VarExp, gdm$explained)
    print(summary(gdm))
  } else if (gdm$explained >= 0.001) {
    print(summary(gdm))
  }
}

###############################################
# DECIDE ON REDUCED/UNCORRELATED VARIABLE SET #
###############################################

# Spearman's rho to check correlations of remaining variables
# i.e. only those that explain > 5% variance
envTab.5Perc.Var = Dh.Env[, c("UTMX", "UTMY", "IDNo", Include)]
# Create correlation plots to help pick variables
ggcorr(envTab.5Perc.Var[, 4:ncol(envTab.5Perc.Var)], method=c("pairwise","spearman"), label= TRUE)

# Create spoke plot to help pick variables containing only those with correlation > or < 0.7/-0.7
correl <- cor(envTab.5Perc.Var[, 4:ncol(envTab.5Perc.Var)], method="spearman")
pos <- correl >= 0.7
neg <- correl <= -0.7
spoke(pos, neg)

# Check contributions (choose variable with greatest contribution in correlated sets)
names(VarExp) <- Include
VarExp[order(VarExp)]

# Use correlation plot to make obvious choices - when variables are correlated > |0.7|, take the one that explains the most variance
rmCor <- c("ADI", "BPC2", "B11", "Fortescue", "B03", "B05", "B31", "B16", "BPC1", "BPC4", "B04", "B32", "MajorRivers", "B07", "ADX", "B13")

# Remove unwanted variables and check no longer correlated
correl <- cor(envTab.5Perc.Var[, -which(colnames(envTab.5Perc.Var) %in% c("UTMX", "UTMY", "IDNo", rmCor))], method="spearman")
pos <- correl >= 0.7
neg <- correl <= -0.7
spoke(pos, neg)

# Reduce to uncorrelated variables in site-pair table
envTab.5Perc.Var = Dh.Env[, colnames(envTab.5Perc.Var)[-which(colnames(envTab.5Perc.Var) %in% rmCor)]]

# Check Variance Inflation Factor
vif(envTab.5Perc.Var[, -which(colnames(envTab.5Perc.Var) %in% c("UTMX", "UTMY", "IDNo"))])
VarExp[names(VarExp) %in% colnames(envTab.5Perc.Var)]
# Start by removing those that explain the least variance
vif(envTab.5Perc.Var[, -which(colnames(envTab.5Perc.Var) %in% c("UTMX", "UTMY", "IDNo", "ADM"))])
vif(envTab.5Perc.Var[, -which(colnames(envTab.5Perc.Var) %in% c("UTMX", "UTMY", "IDNo", "ADM", "RiverBasins"))])
vif(envTab.5Perc.Var[, -which(colnames(envTab.5Perc.Var) %in% c("UTMX", "UTMY", "IDNo", "ADM", "RiverBasins", "BPC3"))])
vif(envTab.5Perc.Var[, -which(colnames(envTab.5Perc.Var) %in% c("UTMX", "UTMY", "IDNo", "ADM", "RiverBasins", "BPC3", "B12"))])
# Now all below 5 - good!
envTab.5Perc.Var = Dh.Env[, c("UTMX", "UTMY", "IDNo", "B06", "B10", "B15", "B29")]


########################
# RUN MULTIVARIATE GDM #
########################

# Rerun GDM with all variables that explain >5% of the variation and that have vif < 5 variables (e.g. multi-variate model)
# Format site-pair table
gdm.dis.5Perc.Var <- formatsitepair(GenDist.BC, bioFormat=3, XColumn="UTMX", YColumn="UTMY", predData=envTab.5Perc.Var, siteColumn="IDNo")

# Reformat any categorical variables (in this case there are none)
for (i in 1:ncol(gdm.dis.5Perc.Var)) {
  if (colnames(gdm.dis.5Perc.Var)[i] %in% paste0("s1.", CatRast)) {
    var <- CatRast[which(colnames(gdm.dis.5Perc.Var)[i] == paste0("s1.", CatRast))]
    s1 <- 0
    s2 <- ifelse(gdm.dis.5Perc.Var[,paste0("s1.", var)] == gdm.dis.5Perc.Var[,paste0("s2.", var)], 0, 1)
    gdm.dis.5Perc.Var[,paste0("s1.", var)] <- s1
    gdm.dis.5Perc.Var[,paste0("s2.", var)] <- s2
  }
} 

# This time, run the model with geographic distance included
gdm.5Perc.Var <- gdm(gdm.dis.5Perc.Var, geo=T)
summary(gdm.5Perc.Var)
plot(gdm.5Perc.Var)


###########################
# VARIABLE IMPORTANCE AND #
#  BACKWARDS ELIMINATION  #
###########################

# Quantify model significance and variable importance/significance in gdm using matrix permutation
# For initial model building, I'm only using 100 permutations
VarImp <- gdm.varImp(gdm.dis.5Perc.Var, geo = T, fullModelOnly = FALSE, nPerm = 100)

# Plot variable importance for full model
par(mfrow = c(1,1))
barplot(sort(VarImp[[2]][,1], decreasing=T), las=2)

# Look at percent deviance explained
# Variable percent deviance explained 
# Variable p values
# Backwards variable elimination
VarImp
# B29 is not significant in fullModel
# Percent deviance explained barely changes between fullModel and fullModel-1, but drops by 1 after fullModel-1
# Therefore, I'm going to drop B29 and go for fullModel-1


################################
# RUN FINAL MULTIVARIATE MODEL #
################################

# Rerun GDM with final variables (based on backwards elimination, i.e. drop B29)
# Reformat site pair table with final variable set
envTab.Fin = envTab.5Perc.Var[, -which(colnames(envTab.5Perc.Var) %in% "B29")]
gdm.diss.Fin.Vars <- formatsitepair(GenDist.BC, bioFormat=3, XColumn="UTMX", YColumn="UTMY", predData=envTab.Fin, siteColumn="IDNo")

# Reformat any Categorical variables as binary response (none here)
for (i in 1:ncol(gdm.diss.Fin.Vars)) {
  if (colnames(gdm.diss.Fin.Vars)[i] %in% paste0("s1.", CatRast)) {
    var <- CatRast[which(colnames(gdm.diss.Fin.Vars)[i] == paste0("s1.", CatRast))]
    s1 <- 0
    s2 <- ifelse(gdm.diss.Fin.Vars[,paste0("s1.", var)] == gdm.diss.Fin.Vars[,paste0("s2.", var)], 0, 1)
    gdm.diss.Fin.Vars[,paste0("s1.", var)] <- s1
    gdm.diss.Fin.Vars[,paste0("s2.", var)] <- s2
  }
} 

# Rerun gdm with final variable set
gdm.diss.Final <- gdm(gdm.diss.Fin.Vars, geo=T)
# Get model summary
summary(gdm.diss.Final)
# Plot I-splines
plot(gdm.diss.Final, plot.layout = c(1, 2), xlim = c(0.4, 0.9), ylim = c(0.4, 0.9))

# Output plots
pdf("IBE/IBE_outputs/Dh_GDM_Ispline.pdf", width = 1.85, height = 1.85, useDingbats = FALSE, pointsize = 2)
par(mfrow = c(1,1))
plot(gdm.diss.Final$predicted, gdm.diss.Final$observed, xlab = "Predicted Environmental Dissimilarity", ylab = "Observed Genetic Dissimilarity", pch = 16, cex = 0.3, lwd = 0.3, col = "grey")
abline(lm(gdm.diss.Final$observed ~ gdm.diss.Final$predicted, data = mtcars), col = "black")
dev.off()

# Output Isplines with error bars/uncertainty based on 1000 bootstrap iterations
options(scipen = 999)
pdf("IBE/IBE_outputs/Dh_GDM_Ispline_1000bs.pdf", width = 3, height = 3, useDingbats = FALSE, pointsize = 5)
plotUncertainty(gdm.diss.Fin.Vars, sampleSites = 0.7, bsIters = 1000, geo = TRUE, splines = NULL, knots = NULL, splineCol="blue", errCol="grey80", plot.linewidth=2.0, plot.layout=c(2,2), parallel=FALSE, cores=2)
dev.off()

# Re-calculate variable importance - only major variables from model testing above (include backwards elimination step so that I can get the percent deviance explained by the geo only model), 500 iterations
VarImp.Final <- gdm.varImp(gdm.diss.Fin.Vars, geo = TRUE, fullModelOnly = FALSE, nPerm = 500)
VarImp.Final

# Plot variable importance
pdf("IBE/IBE_outputs/Dh_GDM_VarImp.pdf", width = 3, height = 4, useDingbats = FALSE, pointsize = 7)
par(mfrow = c(1,1))
barplot(sort(VarImp.Final[[2]][,1], decreasing=T), las=2)
dev.off()

##############
# PLOT MODEL #
##############

# Get raster subset that includes variables in final model
SubsRasts <- All.rasters[[which(names(All.rasters) %in% c("B06", "B10", "B15"))]]
names(SubsRasts)

# Rename enviro variables, then transform using spline functions derived from the final GDM
envTrans <- envTab.Fin
tabTrans <- gdm.transform(gdm.diss.Final, envTrans)
rastTrans <- gdm.transform(gdm.diss.Final, SubsRasts)

# Reduce to principal components and map
rastDat <- na.omit(getValues(rastTrans))
pcaSamp <- prcomp(rastDat)
pcaRast <- predict(rastTrans, pcaSamp, index=1:3)

# Scale rasters (should scale to 255, but too bright - so scaling to 240 instead)
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*240
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*240
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*240

# Check by plotting
par(mfrow= c(1,1))
plotRGB(pcaRast, r=3, g=1, b=2)

# Read in IBRA subregions to add to map
IBRA <- readOGR("Rasters_Shapefiles/IBRA7_Mainland.and.DolphinIslandOnly.shp")
IBRA.Pilb <- IBRA[IBRA$REG_NAME_7 == "Pilbara",]

# Crop and mask to IBRA
pcaRast.c <- crop(pcaRast, IBRA.Pilb)
pcaRast.m <- mask(pcaRast.c, IBRA.Pilb)

# Prepare for plotting
IBRA.Pilb@data$id <- rownames(IBRA.Pilb@data)
IBRA.Pilb.df = fortify(IBRA.Pilb, region = "id")

# Create nice plot for publication
Map.GDM.PCA.plot <- ggRGB(pcaRast.m, r=3, g=1, b=2) +
  geom_polygon(data = IBRA.Pilb.df, 
               aes(x = long, y = lat, group = group), 
               fill = "transparent", col = "lightgrey") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(278000, 1000000), ylim = c(7350000, 7800000)) +
  labs(x = "UTMX", y = "UTMY") +
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
           x.min =  310000, x.max = 1000000,
           y.min = 7380000, y.max = 7782000,
           dist = 100, st.dist = 0.025, st.size = 2.5,
           height = 0.02, border.size = 0.15) +
  north(location = "topleft", scale = 0.05, symbol = 12,
        x.min = 290000, x.max = 1000000, 
        y.min = 7330000, y.max = 7782000)

# View
Map.GDM.PCA.plot

# Save as pdf
ggsave(filename = "GDM.PCA.Map.pdf",
       plot = Map.GDM.PCA.plot,
       device = "pdf",
       path = "IBE/IBE_outputs/",
       scale = 1,
       width = 16.2,
       height = 10.8,
       units = "cm",
       dpi = 300)


##########################
# VARIATION PARTITIONING #
##########################

# The following code runs all combinations of uni/multivariate GDMs to determine the amount of variance explained by each model variable. The final output is a Venn diagram displaying the amount of variance explained (and shared variance across variables). This gives us an idea of how much information each  variable adds to the final model. 

# Note that B15 was not included due to the small amount of variance explained by this variable alone (and it made the Venn diagram super complicated to show overlaps etc. across 4 variables)
# Be prepared for some really hacky code below to get this into the right format for creating the 3-way Venn diagram!!

# Partition variance across all combinations of variables
FinVars <- c("B06", "B10", "Geo")

# Create a loop to get every unique combination of variables
VarComb <- vector()

for (i in 1:length(FinVars)) {
  Var.Comb.Mat <- as.data.frame(combinations(length(FinVars), i, FinVars, repeats.allowed = FALSE))
  VarComb <- c(VarComb, unite(data = Var.Comb.Mat, col = "Comb", sep = "_")$Comb)
}

# Create a data frame with a row each variable combo and a column for the percent deviance explained (to be filled using a loop)
FinalMod.Dat <- data.frame(Variables = VarComb, Percent.Dev = as.numeric(rep(NA, length(VarComb))))

# For each set of variables, run a gdm and add percent deviance to empty dataframe
for (i in 1:length(FinalMod.Dat$Variables)) {
  
  # Get each variable name
  gdmVars <- unlist(str_split(as.character(FinalMod.Dat$Variables[i]), "_"))
  
  # If geo only model, fill df using model from previous final model run
  # The Geo only model is included in the backwards elimination
  if (length(gdmVars) == 1 & gdmVars[1] == "Geo") {
    FinalMod.Dat[i, "Percent.Dev"] <- VarImp.Final[[1]][2, 4]
  } else {
    # Otherwise run the GDM for the specified combination of variables
    # Prepare data
    envTab.all.mods = envTab.Fin[, c("UTMX", "UTMY", "IDNo", gdmVars[!gdmVars %in% "Geo"])]
    # Format site-pair table
    gdm.all.mods.tab <- formatsitepair(GenDist.BC, bioFormat=3, XColumn="UTMX", YColumn="UTMY", 
                                       predData=envTab.all.mods, siteColumn="IDNo")
    
    # If Geo is included in variable list, use "geo=T" in the gdm
    if (any(gdmVars %in% "Geo")) {
      # Run with geographic distance
      gdm.all.mods.G <- gdm(gdm.all.mods.tab, geo=T)
      FinalMod.Dat[i, "Percent.Dev"] <- gdm.all.mods.G$explained
    } else {
      # Otherwise run without geographic distance
      gdm.all.mods <- gdm(gdm.all.mods.tab, geo=F)
      FinalMod.Dat[i, "Percent.Dev"] <- gdm.all.mods$explained
    }
  }
}

# Condense temperature variables so that I have a three-way Venn (can't do proportional Venn with 4 variables)
FinalMod.Dat.Condensed <- FinalMod.Dat
FinalMod.Dat.Condensed$Variables <- str_replace(FinalMod.Dat.Condensed$Variables, "Geo", "Geog. dist")

# Now calculate each segment of deviance explained
FinalMod.Dat.Condensed$Segs <- rep(NA, nrow(FinalMod.Dat.Condensed))
Temp.single <- vector()
Temp.single.overlap <- vector()

# For each set of variables combined, calculate overlap
# Start with single variables
for (i in 1:length(FinalMod.Dat.Condensed$Variables)) {
  # Get each variable name
  gdmVars <- unlist(str_split(as.character(FinalMod.Dat.Condensed$Variables[i]), "_"))
  
  if (length(gdmVars) == 1) {
    # Calculate deviance explained by single variable only
    single <- FinalMod.Dat.Condensed$Percent.Dev[which(sapply(str_split(as.character(FinalMod.Dat.Condensed$Variables), "_"), FUN = length) == 3)] - FinalMod.Dat.Condensed$Percent.Dev[which(!grepl(gdmVars, FinalMod.Dat.Condensed$Variables) & sapply(str_split(as.character(FinalMod.Dat.Condensed$Variables), "_"), FUN = length) == 2)]
    # Add to data frame
    FinalMod.Dat.Condensed[i, "Segs"] <- single
    # Add to vector to make other calculations easier
    Temp.single <- c(Temp.single, single)
    
    # Create a vector with the overlaps within each single var
    Temp.single.overlap <- c(Temp.single.overlap, FinalMod.Dat.Condensed$Percent.Dev[i] - single)
    names(Temp.single.overlap)[length(Temp.single.overlap)] <- gdmVars
  }
}

# Next, calculate 2-way totals
Temp2way <- vector()

for (i in 1:length(FinalMod.Dat.Condensed$Variables)) {
  # Get each variable name
  gdmVars <- unlist(str_split(as.character(FinalMod.Dat.Condensed$Variables[i]), "_"))
  
  if (length(gdmVars) == 2) {
    # Calculate overlap for combined variables
    Temp2way <- c(Temp2way, sum(FinalMod.Dat.Condensed$Percent.Dev[FinalMod.Dat.Condensed$Variables %in% gdmVars]) - FinalMod.Dat.Condensed$Percent.Dev[i])
    names(Temp2way)[length(Temp2way)] <- FinalMod.Dat.Condensed$Variables[i]
  }
}

# Calculate 3-way overlap
Middle <- vector()

for (i in 1:length(Temp.single.overlap)) {
  Middle <- c(Middle, sum(Temp2way[grepl(pattern = names(Temp.single.overlap)[i], x = names(Temp2way))]) - Temp.single.overlap[i])
}

if (length(unique(round(Middle, 6))) == 1) {
  FinalMod.Dat.Condensed$Segs[which(sapply(str_split(as.character(FinalMod.Dat.Condensed$Variables), "_"), FUN = length) == 3)] <- unique(round(Middle, 6))
} else {
  print("Something went wrong!")
  stop()
}

FinalMod.Dat.Condensed$Segs[match(names(Temp2way), FinalMod.Dat.Condensed$Variables)] <- unname(Temp2way - unique(round(Middle, 6)))
round(sum(FinalMod.Dat.Condensed$Segs),4) == round(FinalMod.Dat.Condensed$Percent.Dev[which(sapply(str_split(as.character(FinalMod.Dat.Condensed$Variables), "_"), FUN = length) == 3)],4)

# Add unexplained
FinalMod.Dat.Condensed <- rbind(data.frame(Variables = "Unexplained", Percent.Dev = 100, Segs = (100 - sum(FinalMod.Dat.Condensed$Segs))), FinalMod.Dat.Condensed)

# Create venn object for plot
FinVars.cond <- c("B06", "B10", "Geog. dist")
VennGDM <- Venn(Weight = FinalMod.Dat.Condensed$Segs, SetNames = FinVars.cond)
VennGDM@IndicatorWeight
VennOrder <- c("Unexplained", unname(sapply(apply(VennGDM@IndicatorWeight[,1:length(FinVars.cond)], 1, function(x) names(which(x==1))), paste, collapse="_"))[-1])
VennGDM <- Venn(Weight = FinalMod.Dat.Condensed$Segs[match(VennOrder, as.character(FinalMod.Dat.Condensed$Variables))], SetNames = FinVars.cond)
VennGDM@IndicatorWeight
VennGDM@IndicatorWeight[,ncol(VennGDM@IndicatorWeight)] <- round(VennGDM@IndicatorWeight[,ncol(VennGDM@IndicatorWeight)] , 1)
plot(VennGDM, type = "squares", doWeights = TRUE)

# Save (editing done in illustrator)
pdf("IBE/IBE_outputs/Venn.PercVar.GDM.pdf", width = 5, height = 5, pointsize = 4)
plot(VennGDM, type = "squares", doWeights = TRUE)
dev.off()

# For ease of visualisation, B15 not presented, however it is included in final model
# (sig and explained >2% variance)
# Add note in figure to say unexplained also includes precip
# The other 0.981% explained by precip contained within other variables
VarImp.Final[[1]][2,1] - FinalMod.Dat.Condensed$Percent.Dev[FinalMod.Dat.Condensed$Variables == "B06_B10_Geog. dist"]


