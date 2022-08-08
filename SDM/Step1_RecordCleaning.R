#############################################
#                 STEP 1                    #  
#             Record cleaning               #
#         Data set: Dasyurus hallucatus     #
#            Author: Robyn Shaw             #
#             Date: 11/05/2021              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

#############################
# LOAD IN REQUIRED PACKAGES #
#############################

# Install:
list.of.packages <- c("dplyr",
                      "rgdal",
                      "raster",
                      "corrplot",
                      "enmSdm",
                      "MASS")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Further instructions for installing enmSdm are here:
# https://github.com/adamlilith/enmSdm
# May need to install the package from the zipTarFile 
# and then install the following dependencies: 
# devtools::install_github('adamlilith/omnibus')
# devtools::install_github('adamlilith/statisfactory')
# devtools::install_github('adamlilith/legendary')

library(dplyr)
library(rgdal)
library(raster)
library(corrplot)
library(enmSdm)
library(MASS)


############################
#     Record Cleaning      #
############################

# Read in raw records
RawRecords <- read.csv("SDM/Data/Dh.ALA_NM_Raw.Pilbara.Records.csv", stringsAsFactors = FALSE)

################################
# Remove duplicate coordinates #
################################

# How many records are we starting with?
nrow(RawRecords)
# 5378

# Coordinates
# Are there any records with missing coordinates?
sum(is.na(RawRecords$Lat))
sum(is.na(RawRecords$Long))
sum(RawRecords$Lat == "")
sum(RawRecords$Long == "")
# No missing coordinate info (already removed prior to this script)
# Will look for outliers through mapping later on

# Remove duplicates with exact same coordinates and important meta-data
# Even if these aren't true duplicates, I want to thin to one record per pixel later on, 
# so may as well remove these now
RawRecords.rmDup <- RawRecords[!duplicated(RawRecords[, c("Species", "Long", "Lat", "BasisOfRecord", "Accuracy.m", "Year", "MethodOfCapture", "Certainty")]), ]

# How many records have we removed?
nrow(RawRecords) - nrow(RawRecords.rmDup)
# 2405

# How many records remain?
nrow(RawRecords.rmDup)
# 2973

# Find duplicate coordinates, that have different entries in important meta-data columns
# Then create decision rules to decide which is the best record to keep out of these
# (probably making this unnecessarily complicated, but at least I'll know I've kept the most robust record)
DupCoords <- RawRecords.rmDup[duplicated(RawRecords.rmDup[ , c("Long", "Lat")], fromLast = FALSE) | duplicated(RawRecords.rmDup[ , c("Long", "Lat")], fromLast = TRUE), ]

nrow(DupCoords) # 490 duplicated records
nrow(unique(DupCoords[ , c("Long", "Lat")])) # of 212 unique coordinates

# Create a new Column called "Dup.Group" to identify duplicated coords easily
Dup.Group <- data.frame(lon_lat = unique(paste0(DupCoords$Long, "_", DupCoords$Lat)), Coord.Group = c(paste0("Dup.Group_000", 1:9), paste0("Dup.Group_00", 10:99), paste0("Dup.Group_0", 100:212)))

# Match these to the dup data set
DupCoords$lon_lat <- paste0(DupCoords$Long, "_", DupCoords$Lat)
DupCoords <- left_join(DupCoords, Dup.Group, by = "lon_lat")

# Order by Coord group
DupCoords <- DupCoords[order(DupCoords$Coord.Group), ]

# Set rules for keeping record
DupCoords$KeepDup <- ifelse(is.na(DupCoords$Year), "N", ifelse(DupCoords$Year < 1990, "N", ifelse(is.na(DupCoords$Accuracy.m), "N", ifelse(DupCoords$Accuracy.m > 1000, "N", "Y"))))

# summaries so that 1= all records in group break some/all of the rules vs.
# 2 = one (or more) record is better than the others and should be preferentially kept
Keep.Y.N <- DupCoords %>%
  group_by(lon_lat) %>%
  summarise(KeepDup.count = length(unique(KeepDup)))
# Add to DF
DupCoords <- left_join(DupCoords, Keep.Y.N, by = "lon_lat")

# Split up by keep count and either randomly choose one duplicate to keep or keep the best record
DupCoords.rand <- DupCoords[DupCoords$KeepDup.count == 1, ]
DupCoords.rand <- DupCoords.rand[!duplicated(DupCoords.rand$lon_lat), ]
DupCoords.Choose <- DupCoords[DupCoords$KeepDup.count == 2 & DupCoords$KeepDup == "Y", ]
DupCoords.Choose <- DupCoords.Choose[!duplicated(DupCoords.Choose$lon_lat), ]

KeepRecords <- c(DupCoords.rand$RecordID, DupCoords.Choose$RecordID)

# Remove the rest of the dups
RawRecords.rmDup.2 <- RawRecords.rmDup[!(duplicated(RawRecords.rmDup[ , c("Long", "Lat")], fromLast = FALSE) | duplicated(RawRecords.rmDup[ , c("Long", "Lat")], fromLast = TRUE)), ]
RawRecords.rmDup.2 <- rbind(RawRecords.rmDup.2, DupCoords.rand[, 1:20])
RawRecords.rmDup.2 <- rbind(RawRecords.rmDup.2, DupCoords.Choose[DupCoords.Choose$RecordID %in% KeepRecords, 1:20])

# Any duplicates left?
sum(duplicated(RawRecords.rmDup.2[, c("Long", "Lat")])) # No duplicates

# How many records are left?
nrow(RawRecords.rmDup.2)
# 2695

# How many did we remove?
nrow(RawRecords.rmDup) - nrow(RawRecords.rmDup.2)
# 278


##################################
# Remove missing/erroneous years #
##################################

# Record year
sum(is.na(RawRecords.rmDup.2$Year) | RawRecords.rmDup.2$Year %in% "") 
# 63 with missing year info

# Check these out:
RawRecords.rmDup.2[(is.na(RawRecords.rmDup.2$Year) | RawRecords.rmDup.2$Year %in% ""), ]
# I could possibly look into when the specific surveys were done, but I think it's probably safest to just remove these, since there aren't too many
Records_rmYear <- RawRecords.rmDup.2[!(is.na(RawRecords.rmDup.2$Year)), ]

# Which years are records from?
table(Records_rmYear$Year)
# There is one clearly incorrect year - take a look
Records_rmYear[Records_rmYear$Year == 2099,] 
# Return id is from 2009, I'm going to assume this is the correct year
Records_rmYear$Year[Records_rmYear$Year == 2099] <- 2009
# Leave at this for the moment. Will take a subset later on that represents "contemporary" records

# How many records left?
nrow(Records_rmYear)
# 2632

# How many removed?
nrow(RawRecords.rmDup.2) - nrow(Records_rmYear)
# 63

#############################
# Check coordinate accuracy #
#############################

# Accuracy: 
# Remove or validate records with no accuracy info
# When were records with missing data sampled?
table(Records_rmYear[is.na(Records_rmYear$Accuracy.m), "Year"])
Records_rmYear[is.na(Records_rmYear$Accuracy.m), c("Long", "Lat")]
# All coordinates have 4dps which is ~10 m accuracy (suggesting a GPS was used)
# However, 4 records are from 1990, where it is possible that gps wasn't available
# Will remove all records prior to 2000 to be sure
Records_rmAcc <- Records_rmYear[!(is.na(Records_rmYear$Accuracy.m) & Records_rmYear$Year < 2000), ]
# Assume/replace post 2000 samples are at 1000 m accuracy
Records_rmAcc$Accuracy.m <- ifelse(is.na(Records_rmAcc$Accuracy.m), 1000, Records_rmAcc$Accuracy.m)

# Remove records with accuracy > 5 km for the moment
# I will likely use a higher res, but this is a first pass clean
table(Records_rmAcc$Accuracy.m)
Records_rmAcc[Records_rmAcc$Accuracy.m > 5000, ]
Records_rmAcc <- Records_rmAcc[!(Records_rmAcc$Accuracy.m > 5000), ]
# How many records left?
nrow(Records_rmAcc)
# 2558
# How many removed?
nrow(Records_rmYear) - nrow(Records_rmAcc)
# 74

#######################################
# Check WAM museum IDs for duplicates #
#######################################

# Check for duplicate WAM.IDs from the museum: 
sum(duplicated(Records_rmAcc$WAM.ID[!(is.na(Records_rmAcc$WAM.ID) | Records_rmAcc$WAM.ID %in% c("", " ", "N/A", "n/a", "Scats Only"))]))
DupWAM <- Records_rmAcc$WAM.ID[!(is.na(Records_rmAcc$WAM.ID) | Records_rmAcc$WAM.ID %in% c("", " ", "N/A", "n/a", "Scats Only"))]
DupWAM.Vect <- DupWAM[duplicated(DupWAM)]
DupWAM.df <- Records_rmAcc[Records_rmAcc$WAM.ID %in% DupWAM.Vect, ]
DupWAM.df <- DupWAM.df[order(DupWAM.df$WAM.ID), ]
# Check dups to see if same animal has just been caught multiple times
# Otherwise, if true duplicate, remove
DupWAM.df
# Only two records seem dubious- remove
rmWAM.Dup <- c("Record_5343", "Record_5317")
Records_rmWAM <- Records_rmAcc[!Records_rmAcc$RecordID %in% rmWAM.Dup, ]


# How many were removed?
nrow(Records_rmAcc) - nrow(Records_rmWAM)
# 2
# How many remain?
nrow(Records_rmWAM)
# 2556

#################################
# Check certainty of species ID #
#################################

# Quality of record
# In the regulation 17 returns, the certainty field means:
# 1 Not sure (Identified to family level or lower)
# 2 Moderately certain (Identified to genus level)
# 3 Certain (Identified to species level)
# 4 WAM Vouchered 
table(Records_rmWAM$Certainty)
sum(is.na(Records_rmWAM$Certainty))

# I will delete the records where certainty is listed as 1, 2, "Don't know"/uncertain, etc. and look at the ones with no info. I'll keep moderately certain, as quolls are pretty hard to confuse with other animals if you get a good look
Records_rmCert <- Records_rmWAM[!(Records_rmWAM$Certainty %in% c("1", "2", "Don't know", "Not sure", "possibly Pseudantechinus", "Uncertain")), ]
Records_rmCert[(Records_rmCert$Certainty %in% c("") | is.na(Records_rmCert$Certainty)), ]

unique(Records_rmCert[(Records_rmCert$Certainty %in% c("") | is.na(Records_rmCert$Certainty)), "MethodOfCapture"])
unique(Records_rmCert[(Records_rmCert$Certainty %in% c("") | is.na(Records_rmCert$Certainty)), "BasisOfRecord"])

# Remove records with no certainty info, if no useful info in method of capture or basis of record (that describes the quality of the ID - i.e. did they get a look at it in the hand?)
Records_rmCert <- Records_rmCert[!((Records_rmCert$Certainty %in% c("") | is.na(Records_rmCert$Certainty)) & (Records_rmCert$BasisOfRecord %in% c("NA", "Sighting", "sighting", "", "sighting only") | is.na(Records_rmCert$BasisOfRecord) & (Records_rmCert$MethodOfCapture %in% c("NA", "Sighting", "") | is.na(Records_rmCert$MethodOfCapture)))), ]

# How many were removed?
nrow(Records_rmWAM) - nrow(Records_rmCert)
# 111
# How many remain?
nrow(Records_rmCert)
# 2445


############################################
# Check Method of Capture/ Basis of record #
############################################

table(Records_rmCert$BasisOfRecord)
# Rename common variables
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("cage trap", "cage trapping", "trapping - cage", "Trapping - cage", "trapping-cage", "trapping - cage; searches", "trapping", "trapping - cages and elliots", "Trapping elliott", "Trapping sheffield", "trapping - sheffield")] <- "Cage/Elliot/Unspecified Trap"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("Carcass found on road.", "carcass found on side of road", "Found dead on road")] <- "Road Kill/Carcass"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("Immediate release at capture point", "Released", "Release", "RELEASED")] <- "Capture Release"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("", "N/A")] <- NA
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("observation only", "opportunistic", "Opportunistic", "Opportunistic sighting")] <- "Opportunistic Sighting/Spotlight Survey"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("remote camera", "Remote camera", "Remote camera survey", "video observation only")] <- "Camera Trap or Video"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("Preserved Specimen", "PreservedSpecimen")] <- "Preserved Specimen"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("scat search", "Scat search", "scat searches", "scat searches and trapping")] <- "Sign/Scats"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("searches", "Searches", "Survey")] <- "Unspecified Survey/Search"
Records_rmCert$BasisOfRecord[Records_rmCert$BasisOfRecord %in% c("Radio tracking - VHF, daytime tracking to dens")] <- "Radio Tracking"

table(Records_rmCert$MethodOfCapture)
# Rename common variables
Records_rmCert$MethodOfCapture[Records_rmCert$MethodOfCapture %in% c("1 capture (1 individual)", "2 captures (2 individuals)", "4 captures (4 individuals)", "7 captures (2 individuals)", "7 captures (4 individuals)", "8 captures (8 individuals)", "capture", "Capture", "cage", "Cage trap", "Cage Trap", "Caught or trapped", "Elliot Trap (large)", "ELLIOTT TRAP", "Individual (alive)", "WET PITFALL TRAP")] <- "Cage/Elliot/Unspecified Trap"
Records_rmCert$MethodOfCapture[Records_rmCert$MethodOfCapture %in% c("Dead", "Remains", "ROAD KILL")] <- "Road Kill/Carcass"
Records_rmCert$MethodOfCapture[Records_rmCert$MethodOfCapture %in% c("", "N/A")] <- NA
Records_rmCert$MethodOfCapture[Records_rmCert$MethodOfCapture %in% c("Day sighting", "Fauna: Spotlighting / Nocturnal Survey", "Night sighting", "OBSERVED", "Opportunistic", "Sighting")] <- "Opportunistic Sighting/Spotlight Survey"
Records_rmCert$MethodOfCapture[Records_rmCert$MethodOfCapture %in% c("Camera", "Camera record", "Camera trap", "Fauna: Motion Camera", "Motion camera", "Motion Camera", "MOTION CAMERA", "Motion camera observation", "remote camera ", "Remote Camera")] <- "Camera Trap or Video"
Records_rmCert$MethodOfCapture[Records_rmCert$MethodOfCapture %in% c("Fauna: Tracks / Diggings / Scratching", "Burrow", "Scat", "Scat ", "Scat (new)", "Scat (old)", "Secondary sign", "Track")] <- "Sign/Scats"
Records_rmCert$MethodOfCapture[Records_rmCert$MethodOfCapture %in% c("Fauna: Monitoring Impact Site", "Fauna: Monitoring Reference Site", "IMAGE")] <- "Unspecified Survey/Search"

table(Records_rmCert$BasisOfRecord)
table(Records_rmCert$MethodOfCapture)
table(Records_rmCert$BasisOfRecord, Records_rmCert$MethodOfCapture)
sum(is.na(Records_rmCert$BasisOfRecord) & is.na(Records_rmCert$MethodOfCapture))

# Split into separate dfs for filtering
# df1 - no info on how record was obtained:
df1 <- Records_rmCert[is.na(Records_rmCert$BasisOfRecord) & is.na(Records_rmCert$MethodOfCapture), ]
# df2 - any of the lower quality IDs (i.e. sightings or sign/scats)
df2 <- Records_rmCert[Records_rmCert$BasisOfRecord %in% c("Opportunistic Sighting/Spotlight Survey", "Sign/Scats", "Unspecified Survey/Search") | Records_rmCert$MethodOfCapture %in% c("Opportunistic Sighting/Spotlight Survey", "Sign/Scats", "Unspecified Survey/Search"), ]
# df3 - high quality IDs (captures, camera, specimens)
df3 <- Records_rmCert[!(Records_rmCert$BasisOfRecord %in% c("Opportunistic Sighting/Spotlight Survey", "Sign/Scats", "Unspecified Survey/Search") | Records_rmCert$MethodOfCapture %in% c("Opportunistic Sighting/Spotlight Survey", "Sign/Scats", "Unspecified Survey/Search") | (is.na(Records_rmCert$MethodOfCapture) & is.na(Records_rmCert$BasisOfRecord))), ]

# Look at survey info/collector for df1
table(df1$DataProvider)
table(df1$Database)
unique(df1$SiteName)
unique(df1$Collector)
unique(df1$ApproximateLocation)
table(df1$Certainty) 
sum(is.na(df1$DataProvider))
sum(is.na(df1$Database))
sum(is.na(df1$SiteName))
sum(is.na(df1$Collector))
sum(is.na(df1$ApproximateLocation))
sum(is.na(df1$Certainty))

# No missing info in other fields describing location/collector/etc. Many of the collectors/providers are consultancies and/or familiar names of ecologists/people in dept. All listed as certain or WAM vouchered. I think it's okay to keep all of the records in df1.

# Look at survey info/collector for df2
table(df2$DataProvider)
table(df2$Database)
unique(df2$SiteName)
unique(df2$Collector)
unique(df2$ApproximateLocation)
table(df2$Certainty) 
sum(is.na(df2$DataProvider))
sum(is.na(df2$Database))
sum(is.na(df2$SiteName))
sum(is.na(df2$Collector))
sum(is.na(df2$ApproximateLocation))
sum(is.na(df2$Certainty))

# Again, almost no missing info in other fields describing location/collector/etc and collectors/providers are consultancies and/or familiar names. Many listed as certain/ very certain. I'll just remove those with "Moderately certain".
df2 <- df2[!df2$Certainty == "Moderately certain", ]

# All of df3 can make it through (high quality)
table(df3$BasisOfRecord, df3$MethodOfCapture)

# Merge back together
Records_rmBasisMeth <- rbind(df1, df2, df3)

# How many were removed?
nrow(Records_rmCert) - nrow(Records_rmBasisMeth)
# 7
# How many remain?
nrow(Records_rmBasisMeth)
# 2438

#########################
# LOOK FOR OUTLIERS     #
# Plot records on a map #
#########################

# Load in map (polygon) of Australia's IBRA regions
Aus <- readOGR("Rasters_Shapefiles/IBRA7_subregions_states.shp")

# Convert records to spatial points
Dhal_Records.sp <- SpatialPointsDataFrame(coords = Records_rmBasisMeth[,2:3], data = Records_rmBasisMeth, proj4string = crs(Aus))

# Crop Aus shapefile to records
Aus.cropped <- crop(Aus, extent(Dhal_Records.sp))

# Plot Map
plot(Aus.cropped, col='gray80')
points(Dhal_Records.sp, bg='mediumseagreen', pch=21)

# Most look fine - I'll remove no data values later using rasters (i.e. samples outside of raster extent)


################################################
# Fine-tune cleaning for this project based on #
#            years and resolution              #
################################################

records <- Records_rmBasisMeth
nrow(records) # 2438

# Dates:
# I want this SDM to represent the contemporary distribution in the Pilbara
# Therefore, I think it's probably best to use post 2000 records only

# Look at distribution:
hist(records$Year,
     xlab='Collection year',
     ylab='Number of records',
     main='Collection year (all years)'
)

# Only retain records collected post 2000
# Luckily most records are collected in this timeframe
records <- records[records$Year >= 2000, ]
hist(records$Year,
     xlab='Collection year',
     ylab='Number of records',
     main='Collection year (1970-2000)'
)
nrow(records)
# 2410 

# Check coordinate accuracy
hist(records$Accuracy.m,
     xlab='Accuracy (m)',
     ylab='Number of records',
     main='Coordinate Accuracy'
)

# I'm going to create an SDM with a 1km resolution
# Remove any samples with accuracy > 1km
records <- records[records$Accuracy.m <= 1000, ]
hist(records$Accuracy.m,
     xlab='Accuracy (m)',
     ylab='Number of records',
     main='Coordinate Accuracy'
)
nrow(records)
# 2386 

# Plot records on a map

# Load in map (polygon) of the Pilbra IBRA region:
# UTM CRS, cleaned to remove islands (except for Dolphin Island)
Pilb <- readOGR("Rasters_Shapefiles/PilbaraIBRA7_UTM.shp")

# For the records, first coerce the points into a spatial object
recordsSpatial <- SpatialPointsDataFrame(
  coords=cbind(records$Long, records$Lat), 
  data=records,
  proj4string=CRS('+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs')
)
colnames(recordsSpatial@coords) <- c("Longitude", "Latitude")

# Covert records to UTM50S to match enviro rasters later on
recordsSpatial <- spTransform(recordsSpatial,  CRS("+proj=utm +zone=50 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Add northing and easting to original dataframe
records$UTMX <- recordsSpatial@coords[, 1]
records$UTMY <- recordsSpatial@coords[, 2]

# Change colnames
colnames(recordsSpatial@coords) <- c("UTMX", "UTMY")

# Manually create a bounding box around the points then use that box to zoom in.
plot(recordsSpatial, col='white')
plot(Pilb, add=TRUE, col='gray80')
points(recordsSpatial, pch=21, bg='mediumseagreen')


###################################
#  Load environmental Variables   #
###################################

# Load up rasters
# Note: rasters have been prepared/manipulated prior to this step
# Some rasters were derived from other layers, and they all have been aggregated to the same resolution based on what was deemed most biologically appropriate (mean, median, min, max, etc.), and cropped/masked to the same extent.
Rasters <- stack(list.files("Rasters_Shapefiles/", pattern = ".asc", full.names = TRUE))

# Remove inappropriate layers (i.e. rasters used to test different hypotheses later on)
Rasters <- Rasters[[-which(names(Rasters) %in% c("Dist2Water.1km", "Dist2Water.5km", "Fortescue", "MajorRivers", "RiverBasins", "VRM.5", "VRM.10", "VRM.20"))]]

# Calculate correlations between rasters
# Randomly sample 10,000 points from raster (and get xy coordinates)
points <- xyFromCell(Rasters, sample(ncell(Rasters), 10000))

# Convert to data frame, then spatial points object
points <- as.data.frame(points)
coordinates(points) <- ~x + y
crs(points) <- crs(Rasters[[1]])

# Plot to check
par(mfrow=c(1,1))
plot(Rasters[[1]])
plot(points, add = T)

# Sample enviro variables in raster stack using points
points.corr <- cbind(points@coords, as.data.frame(raster::extract(Rasters, cbind(points@coords))))

# Remove NA points (i.e. in masked out area)
if (any(is.na(rowSums(points.corr)))) {
  points.corr <- points.corr[-which(is.na(rowSums(points.corr))), ]
}

# Plot to check
plot(points.corr[,1:2])

# Test correlation
correl <- cor(points.corr[, -c(1:2)], method='spearman')
diag(correl) <- NA
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)

# Only include variables with a correlation above/below 0.7/-0.7 (otherwise gets too busy to visualise)
CorrNames <- names(as.data.frame(correl))[colMin(as.data.frame(correl)) <= -0.7 | colMax(as.data.frame(correl)) >= 0.7]
correl <- cor(points.corr[, CorrNames], method='spearman')

# Create pairwise correlation plots (jpeg and pdf)
jpeg("SDM/SDM_outputs/RasterCorr_Total1km_0.7Only.jpg", width = 25, height = 25, units = "cm", pointsize = 12, res = 300)
corrplot.mixed(corr = correl, upper = "circle", lower = "number", order = "AOE", lower.col = "black", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = .5, upper.col = c(viridis::viridis(n = 10, option = "D", begin = 0.2, end = 0.9), viridis::viridis(n = 10, option = "A", direction = -1, begin = 0.6)), outline = TRUE)
dev.off()

pdf("SDM/SDM_outputs/RasterCorr_Total1km_0.7Only.pdf", width = 25/cm(1), height = 25/cm(1), pointsize = 12)
corrplot.mixed(corr = correl, upper = "circle", lower = "number", order = "AOE", lower.col = "black", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = .5, upper.col = c(viridis::viridis(n = 10, option = "D", begin = 0.2, end = 0.9), viridis::viridis(n = 10, option = "A", direction = -1, begin = 0.6)), outline = TRUE)
dev.off()


###############################################################
# Filter to one record per cell, match Records to Environment #
#                 and scrutinise for outliers                 #
###############################################################

# Remove all but one record per cell
# use one raster as a "template" for finding cell duplicates
recordsNoDups <- elimCellDups(x=recordsSpatial, r=Rasters[[1]], longLat=c('UTMX', 'UTMY'))
nrow(recordsSpatial)
nrow(recordsNoDups)
# Went from 2386 records to 624 records after spatially thinning to one per pixel

# plot
par(mfrow=c(1,1))
plot(recordsNoDups, col = "white")
plot(Pilb, add=TRUE, col='gray70')
points(recordsNoDups, pch=21, bg='mediumseagreen')

# match species' records with environment at each location
Enviro.Records <- cbind(recordsNoDups@coords, RecordID = recordsNoDups@data$RecordID, as.data.frame(raster::extract(Rasters, cbind(recordsNoDups@coords))))

# Remove any records that might have NA's as environmental data (i.e., over the ocean).
# Remove outlier records
if (any(is.na(rowSums(Enviro.Records[ , 4:ncol(Enviro.Records)])))) {
  Enviro.Records <- Enviro.Records[-which(is.na(rowSums(Enviro.Records[ , 4:ncol(Enviro.Records)]))), ]
}

nrow(Enviro.Records)
# 609

# Originally, another spatial layer was included in the data set which removed 11 additional records due to missing data in some locations (Grazing intensity- see supplementary info). This layer was subsequently dropped, however, the record cleaning preceded this step and so the SDM in the paper doesn't include these records. I have viewed them on a map, and while in unique pixels, they all occur in areas well sampled (so unlikely to be unique environments) - thus I will remove them rather than redoing the SDM.
ExtraNA <- c("Record_1138",
             "Record_2041",
             "Record_2439",
             "Record_3171",
             "Record_5088",
             "Record_2441",
             "Record_2619",
             "Record_2780",
             "Record_2883",
             "Record_3073",
             "Record_3670")

Enviro.Records <- Enviro.Records[!(Enviro.Records$RecordID %in% ExtraNA), ]
nrow(Enviro.Records)
# 598

# Save final presence records
write.csv(Enviro.Records[, 1:2], "SDM/Data/Dh_FinalPresenceCoords.csv", row.names = FALSE)


########################################
#   Clean records for Critical Weight  #
#  Range mammals to create bias layer  #
########################################

# To correct for bias in sampling, we generated background points by sampling 'pseudo-absence' points, i.e. points from more restricted area than 'background'. By creating an absence probability layer based on similar species: Phillips et al. (2009) propose using target background sites, which are presences of other species that could have been recorded given the methods used to observe the focal species. 

# To do this, we downloaded all mammal records from NatureMap and cleaned them using the script below following a similar process to that used for quolls.

# Read in all mammal records from NatureMap
spp.raw <- read.csv("SDM/Data/Mammal.spp_NM_Raw.Records.csv", stringsAsFactors = FALSE)

#########################################
# Remove all samples not in the Pilbara #
#########################################

# Convert data to spatial points object:
spp.raw.sp <- SpatialPointsDataFrame(coords = spp.raw[, c("GDA_LONG", "GDA_LAT")], data = spp.raw, proj4string = CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
plot(spp.raw.sp)

# Remove samples outside of Pilbara
spp.raw.Pilb <- crop(spp.raw.sp, crop(Aus, extent(c(113.6231, 122.7897, -24.90922, -19.53748))))
# Plot
plot(spp.raw.Pilb)

# Save as df
spp.raw.Pilb.df <- spp.raw.Pilb@data


#######################################
# Remove missing coordinates or years #
#######################################

# How many records are we starting with?
nrow(spp.raw.Pilb.df)
# 69412

# Coordinates
# Are there any records with missing coordinates?
sum(is.na(spp.raw.Pilb.df$GDA_LAT))
sum(is.na(spp.raw.Pilb.df$GDA_LONG))
sum(spp.raw.Pilb.df$GDA_LAT == "")
sum(spp.raw.Pilb.df$GDA_LONG == "")
# No missing coordinate info

# Record year
sum(is.na(spp.raw.Pilb.df$COLDATEY) | spp.raw.Pilb.df$COLDATEY %in% "") 
# 1264 with missing year info

# I could possibly look into when the specific surveys were done, but I think this is overkill for such a large dataset. I'll just remove them.
Records_rmYear <- spp.raw.Pilb.df[!(is.na(spp.raw.Pilb.df$COLDATEY) | spp.raw.Pilb.df$COLDATEY %in% ""), ]

# Which years are records from?
table(Records_rmYear$COLDATEY)
# There are some clearly incorrect years - remove these
Records_rmYear <- Records_rmYear[!(Records_rmYear$COLDATEY > 2019), ]


# How many records left?
nrow(Records_rmYear)
# 68143

# How many removed?
nrow(spp.raw.Pilb.df) - nrow(Records_rmYear)
# 1269

#############################
# Check coordinate accuracy #
#############################

# Accuracy: 
# Remove or validate records with no accuracy info
table(Records_rmYear$ACCURACY_M)
sum(is.na(Records_rmYear$ACCURACY_M) | Records_rmYear$ACCURACY_M %in% "")

# When were records with missing data sampled?
table(Records_rmYear[is.na(Records_rmYear$ACCURACY_M), "COLDATEY"])
# There aren't so many with missing info - just delete
Records_rmAcc <- Records_rmYear[!(is.na(Records_rmYear$ACCURACY_M) | Records_rmYear$ACCURACY_M %in% ""), ]

# Remove records with accuracy > 1 km
table(Records_rmAcc$ACCURACY_M)
Records_rmAcc <- Records_rmAcc[!(Records_rmAcc$ACCURACY_M > 1000), ]
# How many records left?
nrow(Records_rmAcc)
# 57696
# How many removed?
nrow(Records_rmYear) - nrow(Records_rmAcc)
# 10447


##########################################
#           Check species                #
##########################################

unique(Records_rmAcc$NAME)
unique(Records_rmAcc$FAMILY)

# Remove bats, marine mammals (that have slipped through), humans, ferals

rmSpp <- c("Molossidae", "Balaenopteridae", "Vespertilionidae", "Dugongidae", "Hominidae", "Phocidae", "Kogiidae", "Delphinidae", "Megadermatidae", "Bovidae", "Physeteridae", "Pteropodidae", "Hipposideridae", "Emballonuridae", "Vombatidae")

Records_rmSpp <- Records_rmAcc[!Records_rmAcc$FAMILY %in% rmSpp, ]

unique(Records_rmSpp$FAMILY)
unique(Records_rmSpp$NAME)

# Look up species and find info on common name, known distribution and weight
spp.info <- read.csv("SDM/Data/MammalSpecies_Info.csv", stringsAsFactors = FALSE)

# Add info to main df
Records_rmSpp <- left_join(Records_rmSpp, spp.info[, c(1, 5:7)], by = "NAME")
head(Records_rmSpp)

# Explore
table(Records_rmSpp$Correct_Distribution)
table(Records_rmSpp$CWR.35_4200g.)

# Remove spp outside of range
Records_rmSpp <- Records_rmSpp[!Records_rmSpp$Correct_Distribution %in% c("Islands", "Edge"), ]

# Remove non CWR mammals
Records_rmSpp <- Records_rmSpp[!Records_rmSpp$CWR.35_4200g. %in% c(">", "<"), ]

# Remove some spp unlikely to be caught in same surveys/traps (e.g. wallabies, water rat, echidna - i.e. those not listed in spp.info.csv)
Records_rmSpp <- Records_rmSpp[paste(Records_rmSpp$GENUS, Records_rmSpp$SPECIES) %in% paste(spp.info$GENUS, spp.info$SPECIES), ]

# Because this is for Dhal bias layer, drop this species from df
Records_rmSpp <- Records_rmSpp[!Records_rmSpp$GENUS %in% c("Dasyurus"), ]

# How many records left?
nrow(Records_rmSpp)
# 18669
# How many removed?
nrow(Records_rmAcc) - nrow(Records_rmSpp)
# 39027

###########################
# Remove pre-2000 samples #
###########################

# I'm going to create a 'current' SD using samples from the year 2000 and above
# Therefore, my bias layer needs to match
Records_current <- Records_rmSpp[Records_rmSpp$COLDATEY >= 2000, ]

# How many records left?
nrow(Records_current)
# 18183
# How many removed?
nrow(Records_rmSpp) - nrow(Records_current)
# 486

############################
# Remove duplicate records #
############################

# Thin samples to one record per species per location
Records_rmDup <- Records_current[!duplicated(Records_current[, c("Common_name", "GDA_LONG", "GDA_LAT")]), ]

# How many records left?
nrow(Records_rmDup)
# 7649
# How many removed?
nrow(Records_current) - nrow(Records_rmDup)
# 10534

###############
# Plot on map #
###############

# Convert data to spatial points object:
Records_rmDup.sp <- SpatialPointsDataFrame(coords = Records_rmDup[, c("GDA_LONG", "GDA_LAT")], data = Records_rmDup, proj4string = CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))

# Plot Map
plot(crop(Aus, extent(c(113.6231, 122.7897, -24.90922, -19.53748))), col='gray80')
points(Records_rmDup.sp, bg='mediumseagreen', pch=21)

#################################
# Create a sampling bias raster #
#################################

# multi-species presence locations
CWR.pres <- Records_rmDup.sp

# Convert to spatial points
CWR.pres.sp <- SpatialPoints(CWR.pres[, 22:23])

# Define CRS (GDA94)
crs(CWR.pres.sp) <- crs(crop(Aus, extent(c(113.6231, 122.7897, -24.90922, -19.53748))))

# Convert to UTM to match enviro layers
CWR.pres.sp.UTM <- spTransform(CWR.pres.sp, CRS("+proj=utm +zone=50 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# Rasterize by counting the number of occurences in each pixel
CWR.pres.ras <- rasterize(CWR.pres.sp.UTM, Rasters[[1]], 1, na.rm = F)
# Change NAs to zeros (i.e. absent)
CWR.pres.ras[is.na(CWR.pres.ras)] <- 0
# Mask to coastline and mask out area outside of study buffer
CWR.pres.ras <- mask(CWR.pres.ras, Rasters[[1]])
# Check rasters match
compareRaster(CWR.pres.ras, Rasters[[1]])

# Check by plotting
plot(CWR.pres.ras)
# Make sure resolution is 1000 m x 1000 m
res(CWR.pres.ras)


# Use the kde2d function from the MASS package, which gives a two-dimensional kernel density estimate, based on the coordinates of the occurrence points. (Note: the output of this function is sensitive to the bandwidth selection; if in doubt, use the default.)

# List cells containing presence
presences <- which(values(CWR.pres.ras) == 1)
# List coordinates of cells containing presence
pres.locs <- coordinates(CWR.pres.ras)[presences, ]

# Create coordinate vectors for each corner of the study area
# This is a work around to make sure the density plot has the same extent as the environmental layers
br <- c(xmax(Rasters[[1]]),ymin(Rasters[[1]])) 
bl <- c(xmin(Rasters[[1]]),ymin(Rasters[[1]]))
tr <- c(xmax(Rasters[[1]]),ymax(Rasters[[1]]))
tl <- c(xmin(Rasters[[1]]),ymax(Rasters[[1]]))

# Add corners to sampling locations
pres.locs <- rbind(pres.locs,br,bl,tr,tl) 

# Run density function
dens <- kde2d(pres.locs[,1],
              pres.locs[,2],
              n = c(nrow(Rasters[[1]]),ncol(Rasters[[1]])))

# Convert density to a raster and resample and mask to raster template
dens.ras <- raster(dens)
dens.ras <- resample(dens.ras, Rasters[[1]], method="ngb")
dens.ras <- crop(dens.ras, extent(Rasters[[1]]))
dens.ras <- mask(dens.ras, Rasters[[1]])

# Check rasters align
compareRaster(dens.ras, Rasters[[1]])

# Check by plotting
plot(dens.ras)

# Save raster
writeRaster(dens.ras, 
            "SDM/Data/CWR_Bias_Layer.asc", 
            format = "ascii", 
            datatype = dataType(dens.ras))

