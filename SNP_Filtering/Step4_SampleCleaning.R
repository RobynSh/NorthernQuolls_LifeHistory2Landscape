#############################################
#                 STEP 4                    #  
# Sample Cleaning: Duplicates and Relatives #
#       Data set: Dasyurus hallucatus       #
#            Author: Robyn Shaw             #
#             Date: 14/02/2022              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

############
# Packages #
############

library(stringr)
library(sp)
library(raster)
library(related)
library(dartR)
library(dplyr)
library(sf)
library(rgdal)
library(RColorBrewer)



###########################
# READ IN DATA - ADD UTMS #
###########################

# Read in meta data
Dh.Metadata <- read.csv("SNP_Filtering/Data/Dasyurus_hallucatus_ind.metadata.csv", stringsAsFactors = FALSE)

# Read in final genlight from 
# Step 3: SNP Filtering
load("SNP_Filtering/Filtering_outputs/gl.Dhal_FinalFilt.rdata")

# Remove poor quality individuals from meta-data that were removed during genlight filtering
Dh.Metadata <- Dh.Metadata[Dh.Metadata$id %in% gl.FinalFilt@ind.names, ]

# Convert long/lats to UTMs
# Need to convert so that distances are in metres for
# Creating buffer zones later on
Spatial.P <- SpatialPoints(Dh.Metadata[ , c("lon", "lat")])
crs(Spatial.P) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
Spatial.P.UTM <- spTransform(Spatial.P, CRS("+proj=utm +zone=50 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
colnames(Spatial.P.UTM@coords) <- c("UTMX", "UTMY")
Dh.Metadata <- cbind(Dh.Metadata, Spatial.P.UTM@coords)


############################
# REMOVE DUPLICATE SAMPLES #
############################

# Look for duplicated samples across metadata cols
Dh.lon_lat_sex_year <- paste0(Dh.Metadata$lon, "_", Dh.Metadata$lat, "_", Dh.Metadata$sex, "_", Dh.Metadata$year)
length(Dh.lon_lat_sex_year) - length(unique(Dh.lon_lat_sex_year))
dups <- Dh.Metadata[(duplicated(Dh.lon_lat_sex_year) | 
                       duplicated(Dh.lon_lat_sex_year, fromLast = TRUE)), ]
dups[order(paste0(dups$lon, "_", dups$lat, "_", dups$sex, "_", dups$year)), ]

# Unfortunately, some duplicated samples have clearly been included by mistake in the first plate sent off
# List and remove them from both metadata and genlight
# Will also remove Dh_14.015 as it's labelled as Dolphin Island but clustering with Chichester.
# I believe it's been mis-labelled as it is very unlikely an individual could disperse this far

Dh.dupID <- c("Dh_14.015",
              "Dh_2016.236", 
              "Dh_2016.150", 
              "Dh_2014.012", 
              "Dh_2016.251", 
              "Dh_16.883_R", 
              "Dh_2016.355", 
              "Dh_2012.410", 
              "Dh_2016.479", 
              "Dh_2012.400", 
              "Dh_2016.910")

Dh.Metadata <- Dh.Metadata[!Dh.Metadata$id %in% Dh.dupID, ]

# Remove from genlight
Dh.gl <- gl.FinalFilt[!(gl.FinalFilt@ind.names %in% Dh.dupID), ]
Dh.gl <- Dh.gl[match(Dh.Metadata$id, Dh.gl@ind.names), ]

# Check that they match (should = 0)
sum(!Dh.gl@ind.names == Dh.Metadata$id)


################################
# Calculate Wang's relatedness #
################################

# Run Dolphin Island separately (unlikely individuals are getting over to mainland):
# Prepare file for input into related
Dh.gl.DI <- Dh.gl[Dh.gl@ind.names %in% rownames(Dh.gl@other$ind.metrics)[Dh.gl@other$ind.metrics$pop == "Dolphin Island"], ]
Dh.demerel.DI <- gl2demerelate(Dh.gl.DI)
Dh.demerel.DI <- Dh.demerel.DI[, 3:ncol(Dh.demerel.DI)]
colnames(Dh.demerel.DI) <- NULL

# Export as a txt file
write.table(Dh.demerel.DI, na = "0", file = "SNP_Filtering/SampleClean_outputs/Dh.related.DI.txt", sep = "\t", row.names = TRUE)

# Read in text file
Dh.rel.tab.DI <- readgenotypedata("SNP_Filtering/SampleClean_outputs/Dh.related.DI.txt")

# Calculate relatedness (island pop - so allow inbreeding)
Dh.rel.DI <- coancestry(Dh.rel.tab.DI$gdata, lynchli = 0, lynchrd = 0, quellergt = 0, ritland = 0, wang = 1, allow.inbreeding = TRUE)

# Add meta-data Ind.1
Dh.rel.DI.tab <- cbind(Dh.rel.DI$relatedness[, c(2:3, 6)], Dh.Metadata[match(Dh.rel.DI$relatedness$ind1.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.DI.tab)[4:6] <- c("Ind.1.sex", "Ind.1.pop", "Ind.1.year")
Dh.rel.DI.tab <- cbind(Dh.rel.DI.tab, Dh.Metadata[match(Dh.rel.DI$relatedness$ind2.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.DI.tab)[7:9] <- c("Ind.2.sex", "Ind.2.pop", "Ind.2.year")

# Output results
write.csv(Dh.rel.DI.tab, "SNP_Filtering/SampleClean_outputs/Dh.relatedness.DI.csv", row.names = FALSE)

# Mainland:
# Prepare file for input into related
Dh.gl.ML <- Dh.gl[Dh.gl@ind.names %in% rownames(Dh.gl@other$ind.metrics)[!Dh.gl@other$ind.metrics$pop == "Dolphin Island"], ]

Dh.demerel.ML <- gl2demerelate(Dh.gl.ML)
Dh.demerel.ML <- Dh.demerel.ML[, 3:ncol(Dh.demerel.ML)]
colnames(Dh.demerel.ML) <- NULL

# Export as a txt file
write.table(Dh.demerel.ML, na = "0", file = "SNP_Filtering/SampleClean_outputs/Dh.related.ML.txt", sep = "\t", row.names = TRUE)

# Read in text file
Dh.rel.tab.ML <- readgenotypedata("SNP_Filtering/SampleClean_outputs/Dh.related.ML.txt")

# Calculate relatedness
Dh.rel.ML <- coancestry(Dh.rel.tab.ML$gdata, lynchli = 0, lynchrd = 0, quellergt = 0, ritland = 0, wang = 1, allow.inbreeding = FALSE)

# Add meta-data Ind.1
Dh.rel.ML.tab <- cbind(Dh.rel.ML$relatedness[, c(2:3, 6)], Dh.Metadata[match(Dh.rel.ML$relatedness$ind1.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.ML.tab)[4:6] <- c("Ind.1.sex", "Ind.1.pop", "Ind.1.year")
Dh.rel.ML.tab <- cbind(Dh.rel.ML.tab, Dh.Metadata[match(Dh.rel.ML$relatedness$ind2.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.ML.tab)[7:9] <- c("Ind.2.sex", "Ind.2.pop", "Ind.2.year")


# Output results
write.csv(Dh.rel.ML.tab, "SNP_Filtering/SampleClean_outputs/Dh.relatedness.ML.csv", row.names = FALSE)

######################
# Remove individuals #
######################

# Decide on individuals to remove in spreadsheets

# For all pops, prioritise 2011 - 2015
# Also, priorities retaining as many samples as possible (so remove samples that pair with multiple individuals)
# Abydos Turner Yandee, Callawa, Cattle Nimingarra Shay, Dolphin Island, Quarry Site 2a, Red Hill Chichester, Yarrie Station: priority= F
# Poondano: priority= M

Dh.rm <- c("Dh_14.008",
           "Dh_16.774",
           "Dh_11.345",
           "Dh_11.334",
           "Dh_16.500",
           "Dh_12.081",
           "Dh_12.083",
           "Dh_12.408",
           "Dh_16.050",
           "Dh_16.101",
           "Dh_16.147",
           "Dh_16.269",
           "Dh_16.252",
           "Dh_16.271",
           "Dh_16.354",
           "Dh_16.356",
           "Dh_11.341",
           "Dh_11.415",
           "Dh_16.069",
           "Dh_16.070",
           "Dh_16.278",
           "Dh_16.330",
           "Dh_16.348",
           "Dh_16.346",
           "Dh_17.074")

# Recalculate relatedness for Dh
Dh.gl.DI <- Dh.gl.DI[!Dh.gl.DI@ind.names %in% Dh.rm, ]
Dh.gl.ML <- Dh.gl.ML[!Dh.gl.ML@ind.names %in% Dh.rm, ]

# Dolphin Island
# Prepare file for input into related
Dh.demerel.DI <- gl2demerelate(Dh.gl.DI)
Dh.demerel.DI <- Dh.demerel.DI[, 3:ncol(Dh.demerel.DI)]
colnames(Dh.demerel.DI) <- NULL

# Export as a txt file
write.table(Dh.demerel.DI,  na = "0", file = "SNP_Filtering/SampleClean_outputs/Dh.related.DI.rerun.txt", sep = "\t", row.names = TRUE)

# Read in text file
Dh.rel.tab.DI <- readgenotypedata("SNP_Filtering/SampleClean_outputs/Dh.related.DI.rerun.txt")

# Calculate relatedness (island pop - so allow inbreeding)
Dh.rel.DI <- coancestry(Dh.rel.tab.DI$gdata, lynchli = 0, lynchrd = 0, quellergt = 0, ritland = 0, wang = 1, allow.inbreeding = TRUE)

# Add meta-data Ind.1
Dh.rel.DI.tab <- cbind(Dh.rel.DI$relatedness[, c(2:3, 6)], Dh.Metadata[match(Dh.rel.DI$relatedness$ind1.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.DI.tab)[4:6] <- c("Ind.1.sex", "Ind.1.pop", "Ind.1.year")
Dh.rel.DI.tab <- cbind(Dh.rel.DI.tab, Dh.Metadata[match(Dh.rel.DI$relatedness$ind2.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.DI.tab)[7:9] <- c("Ind.2.sex", "Ind.2.pop", "Ind.2.year")

# Output results
write.csv(Dh.rel.DI.tab, "SNP_Filtering/SampleClean_outputs/Dh.relatedness.DI.rerun.csv", row.names = FALSE)

# Mainland:
# Prepare file for input into related
Dh.demerel.ML <- gl2demerelate(Dh.gl.ML)
Dh.demerel.ML <- Dh.demerel.ML[, 3:ncol(Dh.demerel.ML)]
colnames(Dh.demerel.ML) <- NULL

# Export as a txt file
write.table(Dh.demerel.ML, na = "0", file = "SNP_Filtering/SampleClean_outputs/Dh.related.ML.rerun.txt", sep = "\t", row.names = TRUE)

# Read in text file
Dh.rel.tab.ML <- readgenotypedata("SNP_Filtering/SampleClean_outputs/Dh.related.ML.rerun.txt")

# Calculate relatedness
Dh.rel.ML <- coancestry(Dh.rel.tab.ML$gdata, lynchli = 0, lynchrd = 0, quellergt = 0, ritland = 0, wang = 1, allow.inbreeding = FALSE)

# Add meta-data Ind.1
Dh.rel.ML.tab <- cbind(Dh.rel.ML$relatedness[, c(2:3, 6)], Dh.Metadata[match(Dh.rel.ML$relatedness$ind1.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.ML.tab)[4:6] <- c("Ind.1.sex", "Ind.1.pop", "Ind.1.year")
Dh.rel.ML.tab <- cbind(Dh.rel.ML.tab, Dh.Metadata[match(Dh.rel.ML$relatedness$ind2.id, Dh.Metadata$id), c("Sex", "pop", "Year")])
colnames(Dh.rel.ML.tab)[7:9] <- c("Ind.2.sex", "Ind.2.pop", "Ind.2.year")


# Output results
write.csv(Dh.rel.ML.tab, "SNP_Filtering/SampleClean_outputs/Dh.relatedness.ML.rerun.csv", row.names = FALSE)




########################################
# ID duplicate coordinates and buffers #
########################################

# Remove related individuals from meta-data and genlight
Dh.Metadata.rmRel <- Dh.Metadata[!Dh.Metadata$id %in% Dh.rm, ]
gl.rmRel <- Dh.gl[!Dh.gl@ind.names %in% Dh.rm, ]

# Add column called CoordLocation 
# That groups samples at exact same coordinates
Dh.Metadata.rmRel <- left_join(Dh.Metadata.rmRel, 
                               cbind(unique(Dh.Metadata.rmRel[, c("UTMX", "UTMY")]),
                                     CoordLocation = 1:nrow(unique(Dh.Metadata.rmRel[, c("UTMX", "UTMY")]))),
                               by = c("UTMX", "UTMY"))


# List buffer sizes (in metres) for grouping samples
Buff_m <- c(1000, 2000, 5000, 10000, 15000, 20000)
names(Buff_m) <- c("Buffer_1km", "Buffer_2km", "Buffer_5km", "Buffer_10km", "Buffer_15km", "Buffer_20km")

# Create spatial points to generate buffers
Coords.utm <- SpatialPoints(Dh.Metadata.rmRel[, c("UTMX", "UTMY")], 
                            proj4string = crs("+proj=utm +zone=50 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# Run a loop that groups samples withing the specified buffer distance and 
# saves unique buffer ID as a column in meta data

for (j in 1:length(Buff_m)) {
    
    # Print loop progress
    print(paste0("Generating ", names(Buff_m)[j]))
    
    # Create a buffer around samples
    # Convert to sf so that I can split polygon layer so each buffer has a unique ID
    UTM.Buff <- st_cast(as(buffer(x = Coords.utm, Buff_m[j], dissolve = TRUE), 
                           "sf"),"POLYGON")
    
    # Name with unique ID
    UTM.Buff$BuffID <- 1:dim(UTM.Buff)[1]
    
    # Convert back to sp object
    # Get buffer polygon info at each sample point
    Dh.Metadata.rmRel <- cbind(Dh.Metadata.rmRel, over(Coords.utm, as(UTM.Buff, "Spatial")))
    colnames(Dh.Metadata.rmRel)[which(colnames(Dh.Metadata.rmRel) == "BuffID")] <- names(Buff_m)[j]
    
}

################
# Save outputs #
################

if ((length(intersect(Dh.Metadata.rmRel$id, gl.rmRel@ind.names)) == nrow(Dh.Metadata.rmRel)) & 
    (length(intersect(Dh.Metadata.rmRel$id, gl.rmRel@ind.names)) == nInd(gl.rmRel))) {
  # Write out meta-data and gl in new folder
  write.csv(Dh.Metadata.rmRel, "SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.Ind.metadata.Dh.csv", row.names = FALSE)
  # Rename R object for saving
  Dh.gl <- gl.rmRel
  save(Dh.gl, file = "SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.gl.Dh.rdata")
  
} else {
  print("Samples don't match between genlight and meta-data")
}




