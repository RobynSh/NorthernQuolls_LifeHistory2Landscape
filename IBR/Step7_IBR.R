#############################################
#                 STEP 7                    #  
#       Isolation-By-Resistance (IBR)       #
#         Data set: Dasyurus hallucatus     #
#            Author: Robyn Shaw             #
#             Date: 22/02/2022              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

#############################
# LOAD IN REQUIRED PACKAGES #
#############################

library(raster)
library(ecodist)
library(reshape2)
library(dplyr)
library(stringr)
library(ResistanceGA)
library(parallel)
library(doParallel)
library(rgdal)
library(pointdexter)
library(viridis)
library(ggplot2)
library(ggnewscale)
library(ggsn)

####################
# SET UP VARIABLES #
####################

# Set species abreviation
Sp <- "Dh"

# Set path to a directory for saving data/outputs
Out.path <- "IBR/IBR_outputs/"


###################################################
# Import layers and convert to appropriate format #
###################################################

### SAMPLE COORDINATES ###
Samps <- read.csv("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.Ind.metadata.Dh.csv", stringsAsFactors = FALSE)

# Remove island samples (I'm not interested in showing that the ocean is/isn't a barrier - I just want to look at what facilitates gene flow on the Pilbara mainland)
Samps <- Samps[!Samps$pop == "Dolphin Island", ]

####################################
# Get pixel coordinate from raster #
####################################

# Read in one of the processed rasters (make sure its one with a 5km pixel size) so that I can aggregate the genetic data based on the pixel size of the rasters 
Raster <- raster("Rasters_5km/Clay_Res5km.asc")

# Get xy coordinates of raster pixels that overlap with quoll samples
cell <- data.frame(id = Samps[, "id"], cell = extract(Raster, Samps[, c("UTMX", "UTMY")], cellnumbers=T)[, 1], Dat = extract(Raster, Samps[, c("UTMX", "UTMY")], cellnumbers=T)[, 2])
# Find any samples outside of raster extent (to remove later)
rmDh <- cell$id[is.na(cell$Dat)] # one sample
xy <- cbind(cell, xyFromCell(Raster, cell[,2]))

# Create a new column called "pixel" to identify duplicated coords easily (i.e. all coordinates that fall within the same 1km pixel)
# Create a new data frame with only the unique coordinates
pixel.group <- data.frame(cell = unique(xy$cell), pixel = c(paste0("P00", 1:9), paste0("P0", 10:length(unique(xy$cell)))))

# Match these to the main data set
Samps <- cbind(Samps, pixel.group[match(xy$cell, pixel.group$cell), "pixel"])
colnames(Samps)[ncol(Samps)] <- "pixel.grp"

# Remove samples outside of raster extent
Samps <- Samps[!Samps$id == rmDh, ]

###################################################
# READ IN FILTERED GENETC DATA - GENLIGHT FORMAT  #
###################################################

# Load in genlight processed in the SNP filtering/Sample cleaning scripts
load("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.gl.Dh.rdata")

# Make sure genlight matches Samps df
gl <- Dh.gl[Dh.gl@ind.names %in% Samps$id, ]

# Make sure ind metrics match ids (for later on)
sum(!gl@ind.names == gl@other$ind.metrics$id) # Should be 0

# Rename Ids as the Coordinate group so it's easier to get mean
indNames(gl) <- Samps$pixel.grp[match(gl@ind.names, Samps$id)]

# Some of the genetic distance methods can't deal with missing data
# I'll interpolate missing data by using the mean
gl.mat <- round(tab(gl, NA.method="mean"), 0)



##############################
# Calculate genetic metrics  #
##############################

# According to Shirk et al. 2017
# https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12684?casa_token=kDgBnlvK8aMAAAAA:cAGe_TzoPPF-u4Qf0p-3WXiB0jhTLo0vUcOV6suW05lCUJifsq8w-nvWECqKeNuigFtiPssskhxAxcuI
# The best metrics for landscape genetic analysis are PC eigenvectors (as long as more than one axis is included), Euclidean distance, Bray-Curtis distance and Proportion of shared alleles (Bray-Curtis and Prop shared alleles were perfectly correlated).
# I ran a preliminary test on all three, and there was very little difference.
# For this reason, I chose to use Euclidean genetic distance and take the average across individuals if there were multiple samples per pixel


#### EUCLIDEAN GENETIC DISTANCE #### 

# Calculate Euclidean genetic distance between samples
EucGenDist <- distance(gl.mat, method = "euclidean")

# Convert lower triangular to table
# Have already tested to make sure Euc dist order is the same as ResistanceGA lower convenience function
# I just need to make sure I sort by Individual 1, then individual 2 in the code below
EucGenDist.mat <- as.matrix(EucGenDist)
EucGenDist.tab <- melt(EucGenDist.mat)[melt(lower.tri(EucGenDist.mat))$value, ]
names(EucGenDist.tab) <- c("Ind2", "Ind1", "EucDist")

# Remove comparisons between same pixel group
EucGenDist.tab <- EucGenDist.tab[!EucGenDist.tab$Ind1 == EucGenDist.tab$Ind2,]

# Create a new column to group by (and take the mean of)
Euc.pixdf <- data.frame(Ind2 = as.numeric(sub(".*P", "", EucGenDist.tab$Ind2)), Ind1 = as.numeric(sub(".*P", "", EucGenDist.tab$Ind1)))

# Make a new column with both pixel groups so I can get the mean across duplicates
EucGenDist.tab$Pair <- paste(ifelse(apply(Euc.pixdf, 1, which.min) == 1, as.character(EucGenDist.tab$Ind2), as.character(EucGenDist.tab$Ind1)), ifelse(apply(Euc.pixdf, 1, which.max) == 1, as.character(EucGenDist.tab$Ind2), as.character(EucGenDist.tab$Ind1)), sep = "_")

# Group by pair and then take the mean
EucGenDist.Av <- EucGenDist.tab %>%
  group_by(Pair) %>%
  summarise(EucDistMean = mean(EucDist))

# Check order is correct
EucGenDist.Av$Pair

# Check correct number of values
(length(unique(Samps$pixel.grp))*length(unique(Samps$pixel.grp))-length(unique(Samps$pixel.grp)))/2 == nrow(EucGenDist.Av)

# Create final vector to run in ResistanceGA
EucDist.Vect <- EucGenDist.Av$EucDistMean


########################################
# WRITE DATA FOR INPUT TO RESISTANCEGA #
########################################

# Write genetic data
genetic.df <- data.frame(EucDist = EucDist.Vect)
write.csv(genetic.df, paste0(Out.path, "GeneticMetrics_", Sp, ".csv"), row.names = FALSE)

# Write coordinates
Coords <- Samps[!duplicated(Samps$pixel.grp), c("UTMX", "UTMY", "pixel.grp")]
Coords <- Coords[order(Coords$pixel.grp), c("UTMX", "UTMY")]

# Check numbers match up 
nrow(genetic.df) == ((nrow(Coords)*nrow(Coords))-nrow(Coords))/2
write.csv(Coords, paste0(Out.path, "Coords_", Sp, ".csv"), row.names = FALSE)



####################
# Run ResistanceGA #
####################

# Set data directory
data.dir <- "IBR/IBR_outputs/"

# Set path/filename to Coordinate file
Coords <- "IBR/IBR_outputs/Coords_Dh.csv"

# Set path/filename to genetic file 
GenMet.df.name <- "IBR/IBR_outputs/GeneticMetrics_Dh.csv"

# Raster folder
Rasts <- "Rasters_5km/"

# Define Genetic Metric (i.e. column name)
GenMet <- "EucDist"

# Set seed for single surface
seed <- 1234

# Replicate (should replicate 2x)
rep <- 1

# Set number of cores to use when running parallel (i.e. number of VCPUs available to instance)
Para <- 16

###############
# Import data #
###############

# Load in sample coordinates
Samps <- as.matrix(read.csv(Coords, stringsAsFactors = FALSE))

# Load in genetic distance (in vector format)
# Already ordered so matches coords
GenMet.df <- read.csv(GenMet.df.name)
Gen_vect <- GenMet.df[, GenMet]


#####################################################
# Run Single Surface Optimisation with ResistanceGA #
#####################################################

# I'm doing an initial single surface optimisation of all rasters to find out which raster performs best when layers are correlated (I'll drop the others for multisurface optimisation so that I'm only including uncorrelated variables)
# Create raster stack (only use those with a 5km resolution))
Rstack <- stack(list.files(Rasts, pattern = "_Res5km.asc", full.names = TRUE))

# Create a vector for the names
Rstack.names <- str_remove(names(Rstack), "_Res5km")
names(Rstack) <- Rstack.names

# Run a loop to analyse each surface separately
for (i in 1:length(Rstack.names)) {
  
  # Create/Set results directory:
  Results.dir <- paste0("IBR/IBR_outputs/SS_Optim/", Rstack.names[i])
  ifelse(!dir.exists(file.path(Results.dir)), dir.create(Results.dir, recursive = TRUE), FALSE)
  
  # Run single surfaces
  GA.inputs <- GA.prep(method = "LL",
                       ASCII.dir = Rstack[[i]],
                       Results.dir = paste0(Results.dir, "/"),
                       seed = seed + rep,
                       parallel = Para)
  
  gdist.inputs <- gdist.prep(n.Pops = nrow(Samps),
                             samples = Samps,
                             response = Gen_vect,
                             method = 'commuteDistance')
  
  # First run all single surfaces, Multi-surface is response variable
  SS_RESULTS <- SS_optim(gdist.inputs = gdist.inputs, 
                         GA.inputs = GA.inputs, 
                         dist_mod = TRUE, 
                         null_mod = TRUE)
  
}


####################
# Run all combined #
####################

# After single surface optimisation, I have decided which (uncorrelated) rasters to include in the final analysis
# I will run the all_comb function for multi-surface optimisation
# Need to be in the results directory or function thros an error
setwd("IBR/IBR_outputs/")

# Directory containing final rasters
Rasts.final <- "../../Rasters_5km/Final_IBR_Set/"


GA.inputs <- GA.prep(method = "LL",
                     ASCII.dir = Rasts.final,
                     Results.dir = "all.comb",
                     seed = 1234 + rep,
                     parallel = Para)

gdist.inputs <- gdist.prep(n.Pops = nrow(Samps),
                           samples = Samps,
                           response = Gen_vect,
                           method = 'commuteDistance')

AC_RESULTS <- all_comb(gdist.inputs = gdist.inputs,
                       GA.inputs = GA.inputs, 
                       results.dir = "all.comb",
                       max.combination = 4,
                       iters = 1000,
                       replicate = 1,
                       sample.prop = 0.75,
                       dist_mod = TRUE,
                       null_mod = TRUE)


########################
# PLOT TRANSFORMATIONS # 
########################

# Plot transformations for layers included in the top-ranked resistance surface.
# Parameter values for transformations are provided in the MLPE output for multisurface optimisations (Multisurface_Optim_Summary.txt)

options(scipen = 999) # remove scientific notation so that axis is in 1000's

# Distance to water
pdf("Dist2WaterTrans.pdf", width = 8, height = 8)
Plot.trans(transformation = 1.50704, PARM = c(6.715469, 1621.32), Resistance = "../../Rasters_5km/Final_IBR_Set/EucDistWater_Res5km.asc")
dev.off()

# Silt
pdf("SiltTrans.pdf", width = 8, height = 8)
Plot.trans(transformation = 1.640336, PARM = c(0.592705, 2352.323), Resistance = "../../Rasters_5km/Final_IBR_Set/Silt_Res5km.asc")
dev.off()


#########################################
# LOAD IN TOP-RANKED RESISTANCE SURFACE #
#    AND PREPARE FOR CIRCUITSCAPE       #
#########################################

# The top-ranked, optimised resistance surface from ResistanceGA was Distance to water and silt
Res <- raster("all.combrep_1/EucDistWater_Total.Silt_Total/EucDistWater_Total.Silt_Total.asc")
plot(Res)

# Load in study region polygon (IBRA borders) 
IBRA <- readOGR("../../Rasters_Shapefiles/IBRA7_Mainland.and.DolphinIslandOnly.shp")

# Subset to the Pilbara
IBRA.Pilb <- IBRA[IBRA$REG_NAME_7 == "Pilbara",]
plot(IBRA.Pilb)

# Prepare nodes for generating current map
# Koen et al. (2014) doi: 10.1111/2041-210X.12197 suggest that using a buffer that is 20% the size of the study width and randomly placing nodes around the perimeter of the study boundary is the best way to reduce bias associated with node placement in a current density map.
# I'm going to crop my current map down to the Pilbara IBRA boundary, so I'll add a 20% buffer around this.
# Note that I'm using the Y axis to determine the 'width' of the study region. I think I'm interpreting the paper recommendations correctly, in that the width (shorter length) is still okay for a rectangular study area (don't really have an option because 20% of the X axis length will put the buffer outside of my raster extent).
IBRA.Pilb.buff <- buffer(x = IBRA.Pilb, (extent(IBRA.Pilb)[4] - extent(IBRA.Pilb)[3])*0.20, dissolve = TRUE)

# Check it all looks okay
plot(IBRA.Pilb.buff)
plot(IBRA.Pilb, add=TRUE)
plot(Res, add=TRUE)

# Add points around the perimeter of the buffer polygon
# I'm sure there's a better way to do this - but I can't figure it out
# Have to convert polygon to raster, then back again!

# Create a raster template
rast <- raster()
extent(rast) <- extent(IBRA.Pilb.buff)
res(rast) <- 5000
values(rast) <- -Inf

# rasterise polygon using template
p <- rasterize(IBRA.Pilb.buff, rast)

# convert back to polygon
pp <- rasterToPolygons(p, dissolve=TRUE)

# Get points along boundary
pp.pointbound <- GetPolygonBoundaries(pp)
# Randomly subsample the points to 100
# Note that Koen et al. used 15-20, but suggest this number will depend on the size of the study area
# They also suggest that apart from computation time, there's no penalty to using more nodes
# Set seed so repeatable
set.seed(126)
pp.points <- pp.pointbound[sample(1:nrow(pp.pointbound), size = 100, replace = FALSE), ]
pp.points.sp <- SpatialPoints(pp.points)

# Check by plotting
plot(pp.points.sp, col= "white")
plot(Res, legend = FALSE, add= TRUE)
plot(pp.points.sp, add = TRUE)

####################
# RUN CIRCUITSCAPE #
####################

# Create circuitscape function
# Point to directory where circuitscape exe is (i.e. software must be installed locally)
# Function by Bernd Gruber
# May turn up in the dartr package at some point - keep an eye out
# https://github.com/green-striped-gecko/dartR

gl.runCS <- function(landscape, locs,outpath=tempdir(), plot=FALSE, CS_exe = 'C:/"Program Files"/Circuitscape/cs_run.exe' )
{
  # Cost surface
  cost <- landscape
  # Locs
  sites <- locs
  
  # Plot it
  if (plot)
  {
    plot(cost)
    points(sites,pch=19,col=2)
  }
  # Rasterize points using the cost extent
  
  sites <- rasterize(x = sites,y = cost)
  # Write rasters to your working directory
  
  writeRaster(sites,file.path(outpath, "sites_rast.asc"),overwrite=TRUE)
  writeRaster(cost,file.path(outpath,"cost_rast.asc"),overwrite=TRUE)
  
  # Make an .ini file
  CS_ini <- c("[circuitscape options]",            
              "data_type = raster",
              "scenario = pairwise",
              "write_cur_maps = 1",
              paste(c("point_file =",
                      "habitat_file =",
                      "output_file ="),
                    paste(outpath,c("sites_rast.asc",
                                    "cost_rast.asc",
                                    "CS.out"),
                          sep="/")))
  
  # Write it to your working directory
  writeLines(CS_ini,file.path(outpath,"myini.ini"))
  
  # Make the CS run cmd
  CS_run <- paste(CS_exe, file.path(outpath,"myini.ini"))
  
  # Run the command
  system(CS_run)
  
  # Import the effective resistance
  CSdis <- as.dist(read.csv(file.path(outpath,"CS_resistances.out"),sep=" ",row.names=1,header=1))
  
  CS_currentmap <- raster(file.path(outpath, "CS_cum_curmap.asc"))
  
  return(list(CSdis=CSdis, map=CS_currentmap))
}

# Run function using the resistance surface, with the border points as the nodes
CS <- gl.runCS(landscape = Res, outpath = tempdir(), locs = pp.points)

# Take the log of the current map (as suggested in Circuitscape)
CSmap <- log(CS$map)

# I'm going to create a slightly smaller buffer to crop my current map with
# This is because my resolution is quite coarse, and the IBRA polygon is very detailed
# Essentially, I just want to smooth out the detail by adding a small buffer
# I'll probably clip to the actual polygon in illustrator
IBRA.Pilb.buff.small <- buffer(x = IBRA.Pilb, (extent(IBRA.Pilb)[4] - extent(IBRA.Pilb)[3])*0.05, dissolve = TRUE)

# Plot to check
plot(IBRA.Pilb.buff)
plot(IBRA.Pilb.buff.small, add = TRUE)
plot(IBRA.Pilb, add = TRUE)

# Crop and mask to remove nodes (i.e. to remove original 20% buffer)
CSmap <- crop(CSmap, IBRA.Pilb.buff.small)
CSmap <- mask(CSmap, IBRA.Pilb.buff.small)


###################################
# CREATE FINAL COMBINED PLOT WITH #
#    INDIVIDUAL HETEROZYGOSITY    #
###################################

# I'm going to colour samples by heterozygosity by locus
# Read in Ind het metrics
IndHet <- read.csv("../../IBB/IBB_outputs/IndHet.df.15km_mean.csv")

# For ggplot, have to convert raster to data frame
map.CS_df <- as.data.frame(CSmap, xy = TRUE) 

# Remove scientific notation from UTMs in plot
options(scipen = 999)

# Create a colour palette
# There are a lot more negative values in this current map, so need to tweak palette so that brighter colours only correspond to positive values (i.e. higher dispersal)
CurrentPal <- c(viridis(n = 20, begin = 0, end = 0.1, option = "A"), viridis(n = 10, begin = 0.1, option = "A"))

# Format IBRA polygon for plotting
IBRA.Pilb@data$id <- rownames(IBRA.Pilb@data)
IBRA.Pilb.df = fortify(IBRA.Pilb, region = "id")

# Create map
Map.CS.plot <- ggplot() +
  geom_raster(data = map.CS_df, 
              aes(x = x, y = y, fill = layer)) +
  geom_polygon(data = IBRA.Pilb.df, 
               aes(x = long, y = lat, group = group), 
               fill = "transparent", col = "lightgrey") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(278000, 1000000), ylim = c(7350000, 7800000)) +
  scale_fill_gradientn(colours = CurrentPal, na.value = "white") +
  labs(x = "UTMX", y = "UTMY", fill = "Current\nDensity\n") +
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
                                    hjust = 0.5,
                                    margin = margin(10,0,0,0), 
                                    face ="bold"),
        legend.position = "right") + 
  new_scale("fill") +
  geom_point(data = IndHet, aes(x = UTMX, 
                                y = UTMY, 
                                fill = PHt.HotCold,
                                size = n), 
             stroke = 0.2, pch = 21, col = "black") +
  scale_fill_viridis(option = "D", begin = 0.2) +
  scale_size(breaks = c(1, seq(10, 30, 10), max(IndHet$n)), range = c(3, 8)) +
  labs(fill = "Genetic\nDiversity\n(Std. Dev.)\n") +
  theme(legend.key=element_blank()) +
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
Map.CS.plot

# Check colours are okay for colourblind
colorblindr::cvd_grid(Map.CS.plot)

# Save as pdf
ggsave(filename = "ResGA.CS.CurrentMap_IndHet.pdf",
       plot = Map.CS.plot, 
       useDingbats=FALSE,
       device = "pdf",
       path = ".",
       scale = 1,
       width = 26,
       height = 15.6,
       units = "cm",
       dpi = 300)


