#######################################
#           Filter on HWE             #
#   Output Genalex formatted file     #
# Calculate individual heterozygosity #
# Interpolate across the landscape    #
#######################################

library(dartR)
library(dplyr)
library(splitstackshape)
library(stringr)
library(gstat)
library(sp)
library(raster)
library(rgdal)
library(ggplot2)
library(ggsn)
library(viridis)

setwd("/Users/RobynShaw/Google Drive/Work/Pilbara_Small_Mammals/Manuscripts/Quoll_Paper/Full_Analysis/IBB/")

############################
# Filter SNPs based on HWE #
############################

# Read in sample names and concensus clusters (from genetic clustering analyses)
IndsClusts <- read.csv("Dh.Structure.Tess.LEA.DAPC.csv", stringsAsFactors = FALSE)

# Load in genlight (post-filtering, related individuals removed) and metadata
load("../Final_Datasets/Cleaned.Unrelated.gl.Dh.rdata")
Dh.MD <- read.csv("../Final_Datasets/Cleaned.Unrelated.Ind.metadata.Dh_New.csv")

Dh.gl$other$ind.metrics <- Dh.MD[match(Dh.gl@ind.names, Dh.MD$id), ]

# Remove dup locations to match IBB analysis
Dh.gl.rmDup <- Dh.gl[Dh.gl@ind.names %in% IndsClusts$sampleID, ]

# Check to see if any individuals were removed
length(intersect(Dh.gl.rmDup@ind.names, IndsClusts$sampleID))
length(setdiff(Dh.gl.rmDup@ind.names, IndsClusts$sampleID))
length(setdiff(Dh.gl.rmDup@ind.names, IndsClusts$sampleID))
# Same order?
sum(!Dh.gl.rmDup@ind.names == IndsClusts$sampleID)

# Add cluster to genlight
Dh.gl.rmDup@pop <- as.factor(IndsClusts$GenClust.Tess[match(Dh.gl.rmDup@ind.names, IndsClusts$sampleID)])

# Remove admixed individuals
Dh.gl.rmAd <- Dh.gl.rmDup[!IndsClusts$Prop.Tess[match(Dh.gl.rmDup@ind.names, IndsClusts$sampleID)] == "Admixed", ]


# Calculate HWE (with Bonferroni Corrected significance)
HWE_Clusts <- gl.report.hwe(Dh.gl.rmAd, bonf = TRUE)
HWE_Clusts$Locus <- as.character(HWE_Clusts$Locus)

# Replace "ns" (non-significant) with a zero, and everything else (i.e. significant loci "*") with a one. Note that I am using the bonferroni significance, rather than just the 0.05 significance which is why there are ns results.
HWE_Clusts$BonSig_1.0 <- ifelse(HWE_Clusts$BonSig == "ns", 0, 1)

# Calculate the proportion of populations where that locus is out of HWE
HWE.sig.pop.count <- HWE_Clusts %>% 
  group_by(Locus) %>% 
  summarise(PropPopsSig.HWE = sum(BonSig_1.0)/length(unique(Dh.gl.rmAd@pop)))

# Create a new data frame with all loci (rather than just the subset out of HWE)
Locus_HWE <- data.frame(Locus = locNames(Dh.gl.rmAd))
Locus_HWE$Locus <- as.character(Locus_HWE$Locus)
Locus_HWE <- left_join(Locus_HWE, HWE.sig.pop.count, by = "Locus")

# Replace NA (i.e. those that were not out of HWE) with a 0
Locus_HWE$PropPopsSig.HWE[is.na(Locus_HWE$PropPopsSig.HWE)] <- 0

# Visualise
hist(Locus_HWE$PropPopsSig.HWE)
sum(Locus_HWE$PropPopsSig.HWE > 0.3) # Loci just out in one pop 
sum(Locus_HWE$PropPopsSig.HWE > 0.6) # Loci out in two pops

# In this case I'm going to remove everything because there are only two (real) pops
gl.HWE <- Dh.gl.rmAd[, !(Dh.gl.rmAd@loc.names %in% Locus_HWE$Locus[Locus_HWE$PropPopsSig.HWE > 0.3])]


############################
# Export in Genalex format #
############################

# Save genetic data as a matrix
GenDat <- as.matrix(gl.HWE[,])

# Change SNP coding so that it is in Genalex format
GenDat[GenDat == 0] <- 1.1
GenDat[GenDat == 1] <- 1.2
GenDat[GenDat == 2] <- 2.2
GenDat[is.na(GenDat)] <- 0
GenDat[1:10,1:10] #Check

#Split so that each locus = two columns
Genalexout <- cSplit(GenDat, c(1:ncol(GenDat)),sep=".") # Split cols
Genalexout[is.na(Genalexout)] <- 0 # Replace NAs in second col with 0

# Check that order matches so I can bring over the metadata
sum(!gl.HWE@other$ind.metrics$id == rownames(GenDat)) # if both dfs match = 0

# Add sample ID and pop (genetic cluster) columns
Genalexout_IDs <- cbind(gl.HWE@ind.names, gl.HWE@pop, Genalexout) 
colnames(Genalexout_IDs)[1:2] <- c("sampleID", "pop")
Genalexout_IDs[1:10,1:10] # Check

# Add a spacer column before/after coordinate info (remember to delete column name in excel)
NA_Col <- rep(NA, nrow(Genalexout_IDs))
Genalex.Final <- cbind(Genalexout_IDs, NA_Col, gl.HWE@other$ind.metrics[, c("UTMX", "UTMY")], NA_Col, gl.HWE@other$ind.metrics[, c("id", "pop", "lat", "lon", "sex", "species", "year", "IBRA", "original_id", "sample_location", "Buffer.15km")])

# Write genalex formatted file
write.csv(x = Genalex.Final, file = "Dh.HWE.Filt.genalex.csv", row.names = FALSE, na = "")

##############################################################
# Calculate internal heterozygosity/homozygosity/relatedness #
##############################################################

# The function GENHET is available here: http://www.aureliecoulon.net/research/ac-computer-programs.html
# This function calculates the five most important individual heterozygosity estimates all at once without any size limitation (loci, alleles or individuals). 

# First create the function (copied over from script provided from the website above)
"GENHET"<- function(dat,estimfreq,locname,alfuser){
    
    nbloc=(ncol(dat)-1)/2
    nbind=nrow(dat)
    
    #estimation of allele frequencies (only if estimfreq=T)
    if(estimfreq=="T")
      
    {
      #creation of the list of alleles
      datv=vector(length=nbind*nbloc*2)
      for (i in 2:ncol(dat)) datv[(nrow(dat)*(i-2)+1):(nrow(dat)*(i-1))]=dat[,i]
      al=sort(na.omit(unique(datv)))
      
      #count of the number of times each allele appears + nb of missing data
      alcount=matrix(nrow=(length(al)+1),ncol=(nbloc+1))
      alcount[,1]=c(al,NA)
      for(j in 1:(nrow(alcount)-1))
        for(k in 1:(ncol(alcount)-1))
          alcount[j,(k+1)]=sum(dat[,(k*2):(k*2+1)]==alcount[j,1],na.rm=T)
      for(l in 2:ncol(alcount))
        alcount[nrow(alcount),l]=(2*nbind-sum(alcount[1:(nrow(alcount)-1),l]))
      
      
      #creation of the table of allele frequencies
      alfreq=matrix(nrow=length(al),ncol=(nbloc+1))
      colnames(alfreq)=c("Allele",locname)
      alfreq[,1]=al
      for(m in (1:nrow(alfreq)))
        for (n in 2:ncol(alfreq)) alfreq[m,n]=alcount[m,n]/(nbind*2-alcount[nrow(alcount),n])
      
    }
    
    else alfreq=alfuser
    
    dat=as.data.frame(dat)
    library(gtools)
    res=matrix(nrow=nrow(dat),ncol=6)
    colnames(res)=c("sampleid","PHt","Hs_obs","Hs_exp","IR","HL")
    res[,1]=as.character(dat[,1])
    
    #estimation of E per locus (for HL and Hs_exp)
    E=vector(length=nbloc)
    alfreq2=alfreq[,2:ncol(alfreq)]*alfreq[,2:ncol(alfreq)]
    for(k in 1:ncol(alfreq2)) E[k]=1-sum(alfreq2[,k],na.rm=T)
    
    #estimation of the mean heterozygosity per locus
    mHtl=vector(length=nbloc)
    ctNAl=0
    ctHtl=0
    for(l in 1:ncol(dat))
    {if (even(l)==T)
    {
      for (m in 1:nrow(dat))
      { if (is.na(dat[m,l])==T) ctNAl=(ctNAl+1)
      else if (is.na(dat[m,(l+1)])==T) ctNAl=(ctNAl+1)
      else if (dat[m,l]!=dat[m,(l+1)]) ctHtl=(ctHtl+1)
      }
      mHtl[l/2]=ctHtl/(nrow(dat)-ctNAl)
      ctNAl=0
      ctHtl=0
    }
    }
    
    #the program in itself
    ctHt=0
    ctNA=0
    ctHm=0
    smHtl=0
    mmHtl=0
    sE=0
    mE=0
    sfl=0
    sEh=0
    sEj=0
    
    for(i in 1:nrow(dat))
    { for (j in 2:(nbloc*2))
    { if (even(j)==T)
    {
      if (is.na(dat[i,j])==T) ctNA=(ctNA+1)
      else if (is.na(dat[i,(j+1)])==T) ctNA=(ctNA+1)
      else {
        if (dat[i,j]!=dat[i,(j+1)])
        {
          ctHt=(ctHt+1)
          sEj=sEj+E[j/2]
        }
        else sEh=sEh+E[j/2]
        smHtl=smHtl+mHtl[j/2]
        sE=sE+E[j/2]
        sfl=sfl+alfreq[alfreq[,1]==as.numeric(dat[i,j]),(j/2+1)]+alfreq[alfreq[,1]==as.numeric(dat[i,j+1]),(j/2+1)]
      }
    }
    }
      res[i,2]=ctHt/(nbloc-ctNA)
      mmHtl=smHtl/(nbloc-ctNA)
      res[i,3]=(ctHt/(nbloc-ctNA))/mmHtl
      mE=sE/(nbloc-ctNA)
      res[i,4]=(ctHt/(nbloc-ctNA))/mE
      ctHm=nbloc-ctHt-ctNA
      res[i,5]=(2*ctHm-sfl)/(2*(nbloc-ctNA)-sfl)
      res[i,6]=sEh/(sEh+sEj)
      ctHt=0
      ctNA=0
      ctHm=0
      smHtl=0
      mmHtl=0
      sE=0
      mE=0
      sfl=0
      sEh=0
      sEj=0
    }
    return(res)
  }


# There are several options:
# homozygosity by locus (HL) calculated in the R package 'genhet' (Coulon, 2010). HL was used to represent genetic diversity in landscape genetic analyses (below) as the HL index weights each locus by its allelic variability to determine its overall contribution to the homozygosity index (Aparicio, Ortego, & Cordero, 2006).
# 
#proportion of heterozygous loci (PHt), standardized heterozygosity relative to the mean expected heterozygosity (Hs_exp, Coltman 1999),
#standardized heterozygosity relative to the mean observed heterozygosity (Hs_obs), IR (Amos 2001), HL (Aparicio 2006)


# Reformat Genalex data frame for GENHET
# Include total dataset but remove loci out of HWE

# In this case I'm going to remove everything because there are only two (real) pops
gl.genhet <- Dh.gl[, !(Dh.gl@loc.names %in% Locus_HWE$Locus[Locus_HWE$PropPopsSig.HWE > 0.3])]

# Save genetic data as a matrix
GenDat.genhet <- as.matrix(gl.genhet[,])

# Change SNP coding so that it is in GENHET format
GenDat.genhet[GenDat.genhet == 0] <- 1.1
GenDat.genhet[GenDat.genhet == 1] <- 1.2
GenDat.genhet[GenDat.genhet == 2] <- 2.2
GenDat.genhet[1:10,1:10] #Check

#Split so that each locus = two columns
genhetout <- cSplit(GenDat.genhet, c(1:ncol(GenDat.genhet)),sep=".") # Split cols
genhetout <- cbind(rownames(GenDat.genhet), genhetout)
genhetout[1:10,1:10] # Check

GenHet <- as.data.frame(genhetout)

# Loci should be in two column format, with each named as a and b
colnames(GenHet) <- str_replace(colnames(GenHet), "_1", "a")
colnames(GenHet) <- str_replace(colnames(GenHet), "_2", "b")
# Column names need to start with a character so R doesn't put an "X" in front
colnames(GenHet) <- c("sampleid", paste0("SNP", colnames(GenHet)[2:ncol(GenHet)]))

# Create data frame with locus names (i.e. without a/b) 
LocNames <- colnames(GenHet)[seq(2, ncol(GenHet), 2)]
LocNames <- str_replace(LocNames, "a", "")

# Run GENHET function
IndHet <- GENHET(dat = GenHet, estimfreq = "T", locname = LocNames)

# Output GENHET results as table
write.table(IndHet,"GENHET_IndHeterozygosity.Dh.txt",sep="\t",row.names=F,quote=F)
# IndHet <- read.table("GENHET_IndHeterozygosity.Dh.txt", header = TRUE, stringsAsFactors = FALSE)


#################################################
# Interpolate internal metrics across landscape #
#################################################

# Convert to data frame
IndHet.df <- data.frame(IndHet, stringsAsFactors = FALSE)
IndHet.df <- cbind(IndHet.df, Dh.gl@other$ind.metrics[match(IndHet.df$sampleid, Dh.gl@other$ind.metrics$id), c("UTMX", "UTMY")])

Buffers <- gl.genhet@other$ind.metrics[, c("id", "pop", "Buffer.1km", "Buffer.2km", "Buffer.5km", "Buffer.10km", "Buffer.15km")]
colnames(Buffers)[1] <- "sampleid"
IndHet.df <- left_join(IndHet.df, Buffers, by= "sampleid")
IndHet.df$Buffer.15km[IndHet.df$pop == "Dolphin Island"] <- "B15.00"

# Group by pair and then take the mean
IndHet.df.mean <- IndHet.df[, c("Buffer.15km", "UTMX", "UTMY", "PHt", "Hs_obs", "Hs_exp", "IR", "HL")] %>%
  group_by(Buffer.15km) %>%
  summarise_all(mean)

# Get SD for means
IndHet.df.mean$PHt.HotCold <- (IndHet.df.mean$PHt - mean(IndHet.df.mean$PHt))/sd(IndHet.df.mean$PHt)
IndHet.df.mean$Hs_obs.HotCold <- (IndHet.df.mean$Hs_obs - mean(IndHet.df.mean$Hs_obs))/sd(IndHet.df.mean$Hs_obs)
IndHet.df.mean$Hs_exp.HotCold <- (IndHet.df.mean$Hs_exp - mean(IndHet.df.mean$Hs_exp))/sd(IndHet.df.mean$Hs_exp)
IndHet.df.mean$IR.HotCold <- (IndHet.df.mean$IR - mean(IndHet.df.mean$IR))/sd(IndHet.df.mean$IR)
IndHet.df.mean$HL.HotCold <- (IndHet.df.mean$HL - mean(IndHet.df.mean$HL))/sd(IndHet.df.mean$HL)

# Get SD for Inds
IndHet.df$PHt.HotCold <- (IndHet.df$PHt - mean(IndHet.df$PHt))/sd(IndHet.df$PHt)
IndHet.df$Hs_obs.HotCold <- (IndHet.df$Hs_obs - mean(IndHet.df$Hs_obs))/sd(IndHet.df$Hs_obs)
IndHet.df$Hs_exp.HotCold <- (IndHet.df$Hs_exp - mean(IndHet.df$Hs_exp))/sd(IndHet.df$Hs_exp)
IndHet.df$IR.HotCold <- (IndHet.df$IR - mean(IndHet.df$IR))/sd(IndHet.df$IR)
IndHet.df$HL.HotCold <- (IndHet.df$HL - mean(IndHet.df$HL))/sd(IndHet.df$HL)


# Write
write.csv(IndHet.df, "IndHet.df15km.csv", row.names = FALSE)


# Write
write.csv(IndHet.df.mean, "IndHet.df.15km_mean.csv", row.names = FALSE)

# Create new SpatialPointsDataFrame object
IndHet.mean.sp <- SpatialPointsDataFrame(coords = IndHet.df.mean[, c("UTMX", "UTMY")], proj4string = CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"), data = IndHet.df.mean)

# Interpolation: Inverse Distance Weighthing (IDW)

# Load in a raster to use as a template
# I'll use the centroids of the pixels as the points for interpolation
# This way the extent/resolution will line up and I'll be able to mask the interpolated surface to the Pilbara using template
raster.template <- raster("../../../../Spatial_layers/Processed/UTM50S_Agg1km.aligned/NaturalPerrenial_EucDistWater_UTM50S.1km.Align.bl.tif")

# Get pixel coordinates
Rast.coords <- as.data.frame(coordinates(raster.template))

# Convert to spatial points
coordinates(Rast.coords) <- c("x", "y")

# Convert to SpatialPixel object
gridded(Rast.coords) <- TRUE

# Convert to SpatialGrid object
fullgrid(Rast.coords) <- TRUE

# Add projection information to the empty grid
proj4string(Rast.coords) <- proj4string(raster.template)

# Load in map (polygon) of the Pilbra IBRA region so I can use this to mask out islands in interpolated maps
# Note: this shape file has been cleaned to remove islands (except for Dolphin Island)
IBRA <- readOGR("../../../../GIS/Australian_shp/IBRA7_subregions_states_Fixed_Geom/IBRA7_subregions_states_fixed_Pilb.clean.shp")

# Run interpolation
Metric <- colnames(IndHet.df.mean[,-c(1:3)])

for (i in 1:length(Metric)) {
  # Interpolate
  idw <- gstat::idw(formula = get(Metric[i]) ~ 1, locations = IndHet.mean.sp, newdata = Rast.coords, idp = 2, nmin = 12)
  
  # Convert to raster and mask to buffer
  r <- raster(idw)
  r <- mask(r, raster.template)
  
  # Mask out islands
  r <- mask(r, IBRA)
  
  # Save in workspace
  assign(paste0("idw.", Metric[i]), idw)
  assign(paste0("Raster.idw.", Metric[i]), r)
  
  rm(r, idw)
}

IBRA.Pilb <- IBRA[IBRA$REG_NAME_7 == "Pilbara", ]
Raster.idw.HL.HotCold.crop <- crop(Raster.idw.HL.HotCold, IBRA.Pilb)
Raster.idw.HL.HotCold.crop <- mask(Raster.idw.HL.HotCold.crop, IBRA.Pilb)

# Make a nice plot for publication
# For ggplot, have to convert raster to data frame
idwHL_df <- as.data.frame(Raster.idw.HL.HotCold.crop, xy = TRUE) 

# Create a colour palette
Pal <- viridis(n = 20, direction = -1, option = "A", end = 0.9)

# Remove scientific notation from UTMs in plot
options(scipen = 999)

HL.plot <- ggplot() +
  geom_raster(data = idwHL_df, 
              aes(x = x, y = y, fill = var1.pred)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(210000, 1070000), ylim = c(7270000, 7830000)) +
  scale_fill_gradientn(colours = Pal, na.value = "white") +
  labs(x = "UTMX", y = "UTMY", fill = "HL") +
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
           x.min =  228000, x.max = 1000000,
           y.min = 7300000, y.max = 7782000,
           dist = 100, st.dist = 0.025, st.size = 2.5,
           height = 0.02, border.size = 0.15) +
  north(location = "topleft", scale = 0.05, symbol = 12,
        x.min = 228000, x.max = 1000000, 
        y.min = 7330000, y.max = 7820000)


HL.plot





# Make a nice plot for publication
# For ggplot, have to convert raster to data frame
idw.PHt_df <- as.data.frame(Raster.idw.PHt, xy = TRUE) 

# Remove scientific notation from UTMs in plot
options(scipen = 999)

PHt.plot <- ggplot() +
  geom_raster(data = idw.PHt_df, 
              aes(x = x, y = y, fill = var1.pred)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(210000, 1070000), ylim = c(7270000, 7830000)) +
  scale_fill_gradientn(colours = rev(Pal), na.value = "white") +
  labs(x = "UTMX", y = "UTMY", fill = "PHt") +
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
           x.min =  228000, x.max = 1000000,
           y.min = 7300000, y.max = 7782000,
           dist = 100, st.dist = 0.025, st.size = 2.5,
           height = 0.02, border.size = 0.15) +
  north(location = "topleft", scale = 0.05, symbol = 12,
        x.min = 228000, x.max = 1000000, 
        y.min = 7330000, y.max = 7820000)


PHt.plot





###########################################################
# Make a nice plot for Mantel test carried out in Genalex #
###########################################################


mantel <- read.csv("Mantel_Dat.csv")
colnames(mantel)[1] <- "Population"
#mantel$Geographic_Distance <- mantel$Geographic_Distance/1000


poplab <- c("Cluster 1", "Cluster 2")
names(poplab) <- c("Cluster1", "Cluster2")

options(scipen = 999)

mantel.plot <- ggplot(mantel, aes(x = Geographic_Distance, y = Genetic_Distance)) +
  geom_point(size = 1, colour = "#4f4f4f")  +
  geom_hline(yintercept = 0, size = 0.2, colour = "black", linetype = "dotted") +
  facet_wrap(~ Population, scales = "free_x", labeller = labeller(Population = poplab)) +
  labs(x = "Geographic Distance (km)", y = "Genetic Distance") +
  coord_cartesian(ylim = c(1000, 3000)) +
  geom_smooth(method='lm', formula= y~x, se = FALSE, size = 0.5) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_text(colour="black", size=9,face="bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(colour="black", size=9, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text.x = element_text(colour="black", size=7,face="plain"),
        axis.text.y = element_text(colour="black", size=7,face="plain"),
        axis.ticks = element_line(size = 0.2),
        strip.text.x = element_text(size=9, face = "bold"),
        strip.background = element_rect(fill = "lightgrey", colour = "black", size = 0.2))
mantel.plot


pdf("Figures/Mantel.Pop.Plots.pdf", width = 8, height = 5, useDingbats=FALSE)
mantel.plot
dev.off()

jpeg("Figures/Mantel.Pop.Plots.jpg", width = 8, height = 5, units = "in", res = 300)
mantel.plot
dev.off()
