#############################################
#                 STEP 5                    #  
#          Isolation-By-Barrier (IBB)       #
#         Data set: Dasyurus hallucatus     #
#            Author: Robyn Shaw             #
#             Date: 22/02/2022              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

#############################
# LOAD IN REQUIRED PACKAGES #
#############################

library(rgdal)
library(raster)
library(ggplot2)
library(dplyr)
library(ggsn)
library(viridis)
library(dartR)
library(colorblindr)
library(ggpubr)
library(tess3r)
library(LEA)
library(scatterpie)
library(stringr)


################################################
#      SET VARIABLEs AND LOAD IN DATA          #
################################################

# Species name
sp <- "Dh"

# Paths
Out.Path <- "IBB/IBB_outputs/"

################
# PREPARE MAPS #
################

# Load in polygon of IBRA subregions
StudyMap.Pilb <- readOGR("Rasters_Shapefiles/PilbaraIBRA.shp")
StudyMap.WA <- readOGR("Rasters_Shapefiles/IBRA_regions_WA.shp")

# Set map limits 
xlims <- c(113.7, 122.5)
ylims <- c(-24, -19.8)

# WA
StudyMap.WA.crop <- crop(StudyMap.WA, c(xlims, ylims))
StudyMap.WA.crop@data$id <- rownames(StudyMap.WA.crop@data)
StudyMap.WA.points = fortify(StudyMap.WA.crop, region = "id")
StudyMap.WA.df <- left_join(StudyMap.WA.points, StudyMap.WA.crop@data, by="id")
StudyMap.WA.df <- StudyMap.WA.df[, c(7, 9, 11, 1:2)]

# Pilbara
StudyMap.Pilb@data$id <- rownames(StudyMap.Pilb@data)
StudyMap.Pilb.points = fortify(StudyMap.Pilb, region = "id")
StudyMap.Pilb.df <- left_join(StudyMap.Pilb.points, StudyMap.Pilb@data, by="id")
StudyMap.Pilb.df <- StudyMap.Pilb.df[, c(7, 14, 12, 1:2)]

# Rename regions/sub regions so that they are ordered for plotting (so I can assign the correct colours)
StudyMap.WA.df$REG_NAME <- as.character(StudyMap.WA.df$REG_NAME)
StudyMap.WA.df <- StudyMap.WA.df[order(StudyMap.WA.df$REG_NAME), ]
Reg.name <- unique(StudyMap.WA.df$REG_NAME)

# All of this is just a hacky way of renaming regions so that they are ordered first and get assigned the grey-scale colours once I add in colours for sample locations
PlotOrder <- vector()
Reg <- 0

for (Reg in 0:(length(Reg.name)-1)) {
  RegNo <- Reg + 1
  PlotOrder <- c(PlotOrder, rep((RegNo), length(StudyMap.WA.df$REG_NAME[StudyMap.WA.df$REG_NAME == Reg.name[RegNo]])))
}
StudyMap.WA.df$PlotReg <- as.factor(paste0(0, PlotOrder, "_", StudyMap.WA.df$REG_NAME))
unique(StudyMap.WA.df$PlotReg)

StudyMap.Pilb.df$SUB_NAME_7 <- as.character(StudyMap.Pilb.df$SUB_NAME_7)
StudyMap.Pilb.df <- StudyMap.Pilb.df[order(StudyMap.Pilb.df$SUB_NAME_7), ]
Sub.name <- unique(StudyMap.Pilb.df$SUB_NAME)

PlotOrder_Sub <- vector()
SubReg <- 0

for (SubReg in 0:(length(Sub.name)-1)) {
  SubRegNo <- SubReg + 1
  PlotOrder_Sub <- c(PlotOrder_Sub, rep((SubRegNo), length(StudyMap.Pilb.df$SUB_NAME[StudyMap.Pilb.df$SUB_NAME == Sub.name[SubRegNo]])))
}
StudyMap.Pilb.df$PlotSubReg <- as.factor(paste((PlotOrder_Sub + length(unique(StudyMap.WA.df$PlotReg))), StudyMap.Pilb.df$SUB_NAME, sep = "_"))
unique(StudyMap.Pilb.df$PlotSubReg)

# Create palette for background map
MapPal <- c(gray.colors(10)[seq(2,10,2)], gray.colors(10)[seq(1,10,2)])

# Plot
Map.plot <- ggplot() +
  geom_polygon(data = StudyMap.WA.df, 
               aes(x = long, y = lat, group = group, fill = PlotReg)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  geom_polygon(data = StudyMap.Pilb.df, 
               aes(x = long, y = lat, group = group, fill = PlotSubReg)) + 
  scale_fill_manual(values = c(MapPal)) +
  xlab("Longitude") +
  ylab("Latitude") +
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
        legend.position = "none") + 
  scalebar(location = "bottomleft", dist_unit = "km", 
           transform = TRUE, model = 'WGS84',
           x.min = xlims[1] + 0.25, x.max = xlims[2],
           y.min = ylims[1] + 0.25, y.max = ylims[2],
           dist = 100, st.dist = 0.025, st.size = 2.5,
           height = 0.02, border.size = 0.15) +
  north(location = "topleft", scale = 0.06, symbol = 12,
        x.min = xlims[1] + 0.1, x.max = xlims[2], 
        y.min = ylims[1], y.max = ylims[2] - 0.1)



########################
# PREPARE GENETIC DATA #
########################

# This script takes a genlight as input - I'm using the filtered and cleaned genlight from steps 3-4
load("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.gl.Dh.rdata")

# Read in metadata
Dh.MD <- read.csv("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.Ind.metadata.Dh.csv")

# Make sure data sets match/add in meta-data to genlight in "other" slot
Dh.gl@other$ind.metrics <- Dh.MD[match(Dh.gl@ind.names, Dh.MD$id), ]

# Add 15km buffer to pop slot in genlight
# Pre-pend "B15" to numbers
BuffVec.order <- cbind(unique(Dh.MD$Buffer_15km[order(Dh.MD$lon, Dh.MD$lat)]), c(paste0("B15.0", 1:9), paste0("B15.", 10:length(unique(Dh.MD$Buffer_15km[order(Dh.MD$lon, Dh.MD$lat)])))))

# Add to pop slot
Dh.gl@pop <- as.factor(BuffVec.order[,2][match(Dh.gl@other$ind.metrics$Buffer_15km, BuffVec.order[,1])])

# Change factor for Dolphin island individuals so that they are grouped together (15 kms groups them with some mainland samples)
levels(Dh.gl@pop) <- c(levels(Dh.gl@pop), "B15.00")
Dh.gl@pop[Dh.gl@other$ind.metrics$pop == "Dolphin Island"] <- "B15.00"

# Subset dataset so that there is only one sample per location
Dh.gl.rmDupCoords <- Dh.gl[!duplicated(Dh.gl@other$latlong), ]

# Create palette and shape vectors
Pal <- c(viridis(n = 11, option = "D", begin = 0.2), viridis(n = 10, option = "A", direction = -1, begin = 0.5))
Shape <- rep(21:25, 5) 

################################
#          RUN A PCOA          #
################################

# Run a PCoA and output plots
pc <- gl.pcoa(Dh.gl.rmDupCoords)

# Calculate the percentage of variation represented by the axes
Eig <- round(pc$eig * 100/sum(pc$eig), 1)
barplot(height = Eig, cex.axis = 0.8) 

# Compare axes 
# In this case PC1, PC2 and PC3
# Define variables for Pcoa plotting below:
xaxis <- 1
yaxis <- 2
zaxis <- 3
df <- data.frame(cbind(pc$scores[, xaxis], pc$scores[, yaxis], pc$scores[, zaxis]))
xlab <- paste0("PC ", xaxis, " (", Eig[xaxis], "%)")
ylab <- paste0("PC ", yaxis, " (", Eig[yaxis], "%)")
zlab <- paste0("PC ", zaxis, " (", Eig[zaxis], "%)")
ind.pc <- indNames(Dh.gl.rmDupCoords)
pop.pc <- paste0("X", pop(Dh.gl.rmDupCoords)) # Add "x" to the start so plotting order is correct for colours in the map
lat.pc <- Dh.gl.rmDupCoords@other$latlong$lat
long.pc <- Dh.gl.rmDupCoords@other$latlong$lon
df <- cbind(df, ind.pc, as.factor(pop.pc), lat.pc, long.pc)
colnames(df) <- c("PCoAx", "PCoAy", "PCoAz", "ind", "pop", "lat", "long")
PopNo <- length(unique(Dh.gl.rmDupCoords@pop))

# Create plots:
Pcoa_1V2 <- ggplot(df, aes(x = PCoAx, y = PCoAy, fill = pop, shape = pop)) + 
  geom_hline(yintercept = 0, lwd = 0.2, linetype = "dashed") + 
  geom_vline(xintercept = 0, lwd = 0.2, linetype = "dashed") +
  geom_jitter(size = 3, colour = "black", stroke = 0.2, width = 0.2, height = 0.6) + 
  scale_fill_manual(values = Pal) +
  scale_shape_manual(values = Shape) +
  labs(x = xlab, y = ylab) +

  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  
        axis.title.x = element_text(colour="black", face = "bold", 
                                    size=8, margin = margin(10,0,0,0)),
        axis.title.y = element_text(colour="black", face = "bold", 
                                    size=8, margin = margin(0,10,0,0)),
        axis.ticks = element_line(size = 0.2),
        axis.text = element_text(colour="black",size=6),
        legend.position = "none")

Pcoa_1V3 <- ggplot(df, aes(x = PCoAx, y = PCoAz, fill = pop, shape = pop)) + 
  geom_hline(yintercept = 0, lwd = 0.2, linetype = "dashed") + 
  geom_vline(xintercept = 0, lwd = 0.2, linetype = "dashed") +
  geom_point(size = 3, colour = "black", stroke = 0.2) + 
  scale_fill_manual(values = Pal) +
  scale_shape_manual(values = Shape) +
  labs(x = xlab, y = zlab) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  
        axis.title.x = element_text(colour="black", face = "bold", 
                                    size=7, margin = margin(5,0,0,0)),
        axis.title.y = element_text(colour="black", face = "bold", 
                                    size=7, margin = margin(0,5,0,0)),
        axis.ticks = element_line(size = 0.2),
        axis.text = element_text(colour="black",size=5),
        legend.position = "none")

Pcoa_2V3 <- ggplot(df, aes(x = PCoAy, y = PCoAz, fill = pop, shape = pop)) + 
  geom_hline(yintercept = 0, lwd = 0.1, linetype = "dashed") + 
  geom_vline(xintercept = 0, lwd = 0.1, linetype = "dashed") +
  geom_point(size = 3, colour = "black", stroke = 0.2) + 
  scale_fill_manual(values = Pal) +
  scale_shape_manual(values = Shape) +
  labs(x = ylab, y = zlab) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),  
        axis.title.x = element_text(colour="black", face = "bold", 
                                    size=7, margin = margin(5,0,0,0)),
        axis.title.y = element_text(colour="black", face = "bold", 
                                    size=7, margin = margin(0,5,0,0)),
        axis.ticks = element_line(size = 0.2),
        axis.text = element_text(colour="black",size=5),
        legend.position = "none")


Samp.Map <- Map.plot +
  scale_fill_manual(values = c(MapPal, Pal)) +
  guides(fill = "none") +
  geom_point(data = df, aes(x = long, 
                            y = lat, 
                            fill = pop, 
                            shape = pop), 
             size = 4, stroke = 0.2, colour = "black") + 
  scale_shape_manual(values = Shape) +
  theme(legend.position = "none")
Samp.Map


EigPlot <- ggplot(data.frame(PC = 1:length(Eig), Eig = Eig), aes(x = PC, y = Eig)) +
  geom_bar(stat = "identity", fill = "darkgrey", width = 1) +
  ylab(NULL) +
  xlab(NULL) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.1),
        panel.spacing = unit(c(0,0,0,0), "null"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Combine plots
PCoA.Eig <- ggarrange(Samp.Map, ggarrange(Pcoa_1V2, Pcoa_1V3, ncol = 2), ncol = 1, heights = c(5/8,3/8), labels = c("a)", "b)"), font.label = list(size = 8, face = "plain"))

# Check that figures are okay for colourblind
cvd_grid(PCoA.Eig)

# Make eigenvector plot window to print on plot
vp <- viewport(width = 0.178, height = 0.12, x = 0.395, y = 0.3)

# Output to pdf
pdf(paste0(Out.Path, sp, "_Map.PCoA.pdf"), width = 7.1, height = 6.5, useDingbats = FALSE)
print(PCoA.Eig)
print(EigPlot, vp = vp)
dev.off()



##########################
#  Run Tess3R for Total #
##########################

# Define maximum value of K to test
# From the PCoA it doesn't look like there could be much more than 4 groups, so I'll test a maximum of 7
Kmax <- 7

# Prepare data for Tess
Geno.Tess3R <- as.matrix(Dh.gl.rmDupCoords)
Coords.Tess3R <- as.matrix(Dh.gl.rmDupCoords@other$latlong[, c(2,1)])
Pops.Tess3R <- as.character(Dh.gl.rmDupCoords@pop)

# Run Tess3R
set.seed(123)
Tess.obj <- tess3(X = Geno.Tess3R, coord = Coords.Tess3R, K = 1:Kmax, rep = 100, method = "projected.ls", ploidy = 2, openMP.core.num = 4, max.iteration = 1000, keep = "best", mask = 0.1) 

# Save Tess object
save(Tess.obj, file = paste0(Out.Path, "Tess.obj.FullDat.best100reps.K1_", Kmax, ".", sp, ".R"))

# NOTE: Unfortunately, I forgot to set the seed when I ran this originally. This means that the Tess run above is slightly different (just due to stochasticity across runs). Therefore, the run provided in the folder IBB_outputs is the run used in the paper, which I load in below. The parameters etc. used in this run are identical to that documented above.
# load("IBB/IBB_outputs/Tess.obj.FullDat.best100reps.K1_7.Dh.R")

###################################
# Cross validation plot for total #
###################################

# Get cross validation score in format for plotting:
Tess.CrossVal <- data.frame(Iteration = 1:100)
TessCrossvalMeans <- data.frame(PlotOrder = 1:Kmax, K = paste0("K", 1:Kmax), MeanRMSE = rep(NA, Kmax), Min = rep(NA, Kmax), Max = rep(NA, Kmax))
Kval <- 0

for(Kval in 0:(length(Tess.obj)-1)){
  K <- Kval + 1
  CrossScore <- Tess.obj[[K]]$rmse
  Tess.CrossVal$Total <- CrossScore
  Mean <- TessCrossvalMeans[K, "MeanRMSE"] <- mean(Tess.CrossVal$Total)
  TessCrossvalMeans[K, "Min"] <- min(Tess.CrossVal$Total)
  TessCrossvalMeans[K, "Max"] <- max(Tess.CrossVal$Total)
}

CrossValPlot.Tess <- ggplot(TessCrossvalMeans, aes(x=PlotOrder, y=Min)) + 
  geom_point(size = 0.8) +
  geom_line(size = 0.2) +
  geom_errorbar(aes(x = PlotOrder, ymin=Min, ymax=Max), width = 0.1, size = 0.2) +
  labs(x = "K", y = "Cross Validation Score") +
  scale_x_continuous(breaks = seq(1,Kmax,1)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5),
        axis.title.y = element_text(size = 8, margin = margin(0,10,0,0), face ="bold"),
        axis.title.x = element_text(size = 8, margin = margin(10,0,0,0), face ="bold"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none")

# Save as pdf
ggsave(filename = paste0("CrossVal.Tess.", sp, ".pdf"),
       plot = CrossValPlot.Tess,
       device = "pdf",
       path = Out.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# cross val plot for multipanel figure
CrossValPlot.Tess.MP <- ggplot(TessCrossvalMeans, aes(x=PlotOrder, y=Min)) + 
  geom_point(size = 0.8) +
  geom_line(size = 0.2) +
  geom_errorbar(aes(x = PlotOrder, ymin=Min, ymax=Max), width = 0.1, size = 0.2) +
  scale_x_continuous(breaks = seq(1,Kmax,1)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")


#####################
# Run LEA for Total #
#####################

# LEA requires the path to a Geno file, so I'll have to write the Geno matrix (prepared during the Tess run) to file
write.geno(R = Geno.Tess3R, output.file = paste0(Out.Path, sp, ".geno"))

# Run LEA Algorithm
snmf.obj <- snmf(input.file = paste0(Out.Path, sp, ".geno"), K=1:Kmax, alpha = 10, tolerance = 0.00001, percentage = 0.05, iterations = 1000, ploidy=2, entropy=T, project="new", repetitions = 100, seed = 123, CPU = 2) 

# Save snmf object
save(snmf.obj, file = paste0(Out.Path, "snmf.obj.FullDat.best100reps.K1_", Kmax, ".", sp, ".R"))
# NOTE: Unfortunately, I forgot to set the seed when I originally ran snmf. This means that the run above is slightly different (just due to stochasticity across runs). Therefore, the run provided in the folder IBB_outputs is the run used in the paper, which I load in below. The parameters etc. used in this run are identical to that documented above.
# load("IBB/IBB_outputs/snmf.obj.FullDat.best100reps.K1_7.Dh.R")

###################################
# Cross validation plot for total #
###################################

# Get cross validation score in format for plotting:

# Get cross validation score in format for plotting:
# Get summary of the sNMF project
# This will output cross entropy values for each K (min, mean, max)
LEA.sum <- summary(snmf.obj)
LEA.CrossEnt <- as.data.frame(t(LEA.sum$crossEntropy))
LEA.CrossEnt$PlotOrder <- 1:Kmax

# Plot cross-validation

CrossValPlot.LEA <- ggplot(LEA.CrossEnt, aes(x=PlotOrder, y=min)) + 
  geom_point(size = 0.8) +
  geom_line(size = 0.2) +
  geom_errorbar(aes(x = PlotOrder, ymin=min, ymax=max), width = 0.1, size = 0.2) +
  labs(x = "K", y = "Cross Validation Score") +
  scale_x_continuous(breaks = seq(1,Kmax,1)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5),
        axis.title.y = element_text(size = 8, margin = margin(0,10,0,0), face ="bold"),
        axis.title.x = element_text(size = 8, margin = margin(10,0,0,0), face ="bold"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none")

# Save as pdf
ggsave(filename = paste0("CrossVal.LEA.", sp, ".pdf"),
       plot = CrossValPlot.LEA,
       device = "pdf",
       path = Out.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)


#################
# Tess BAR PLOT #
#################

# Define value for K 
K <- 2

# Retrieve Tess3 Qmatrix rep with the lowest cross-entropy criterion
Q.matrix.Tess <- qmatrix(Tess.obj, K = K)

# Reformat data for plotting
# Create a new column by combining longs and lats 
# Order by this column so bar chart matches up with map
Q.matrix.Tess.df <- cbind(sampleID = Dh.gl.rmDupCoords@other$ind.metrics$id, 
                          lon = Dh.gl.rmDupCoords@other$latlong$lon, 
                          lat = Dh.gl.rmDupCoords@other$latlong$lat, 
                          as.data.frame(Q.matrix.Tess[,]))

# Assign individuals to clusters based max ancestry proportion
# If ancestry proportion <0.7 mark as admixed
TessAssign <- Q.matrix.Tess.df
# I want to name cluster 1 as the western-most population 
# Not necessary, but I prefer it aesthetically
# Check which column a western sample assigns to
colnames(TessAssign[order(TessAssign$lon)[1],c(4:(4+K-1))])[apply(TessAssign[order(TessAssign$lon)[1],c(4:(4+K-1))],1,which.max)]
# "V2" - therefore "V2" will be labelled "Cluster 1" and "V1" will be "Cluster 2" (i.e. reverse the order below)
colnames(TessAssign)[c(4:(ncol(TessAssign)))] <- paste0("TessClust", rev(1:K))
Tess.Admixed <- apply(TessAssign[,c(4:ncol(TessAssign))],1,max) < 0.7
TessAssign$TessAdmixed <- ifelse(Tess.Admixed == TRUE, "Y", "N")
TessAssign$TessGenCluster <- colnames(TessAssign[,c(4:(4+K-1))])[apply(TessAssign[,c(4:(4+K-1))],1,which.max)]


TessAssign <- left_join(TessAssign,
                        TessAssign %>%
                          group_by(TessGenCluster) %>%
                          summarise(MeanLon = mean(lon)),
                        by = "TessGenCluster")

# Now sort on coordinates
Qmat.K.sorted <-  TessAssign[order(TessAssign$MeanLon, 
                                   TessAssign[,unique(TessAssign$TessGenCluster[TessAssign$MeanLon == max(TessAssign$MeanLon)])]), 
                             -which(colnames(TessAssign) %in% c("MeanLon", "lon", "lat"))]

# Add new dummy IDs that are in the correct order for plotting
# Adjust maximum number/numbr of zeros acordsingly depending on how many samples
Qmat.K.sorted$PlotID <- c(paste0("Dummy_00", 1:9), 
                          paste0("Dummy_0", 10:99),
                          paste0("Dummy_", 100:nrow(Qmat.K.sorted)))

# Reformat for bar chart
Qmat.K <- tidyr::gather(Qmat.K.sorted, 
                        "popGroup", 
                        "prop",
                        -sampleID,
                        -PlotID,
                        -TessAdmixed,
                        -TessGenCluster)


# Save in workspace for combined plot later on
TessBarPlot <- ggplot(Qmat.K, aes(factor(PlotID), prop, fill = factor(popGroup))) +
  scale_fill_manual(values = Pal[c(1, 19)]) +
  geom_col(color = "transparent", size = 0.1, width = 1) +
  labs(y = " ") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "plain", size = 6, angle = 90),
        strip.background = element_rect(linetype = "solid", size = 0.1, 
                                        colour = "black", fill = "lightgrey"),
        axis.text = element_blank(),
        axis.title.y = element_text(size = 9, margin = margin(0,24,0,0), face ="bold"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Add Tess info into ind.metric slot of genlight
Dh.gl.rmDupCoords@other$ind.metrics <- cbind(Dh.gl.rmDupCoords@other$ind.metrics,
                    TessAssign[match(Dh.gl.rmDupCoords@other$ind.metrics$id, TessAssign$sampleID), 
                               -which(colnames(TessAssign) %in% c("sampleID", "lon", "lat", "MeanLon"))])

######################
# Tess PIE CHART MAP #
######################

# Get mean of individual coordinates so we have a single set of coordinates for each pop
coord.pop.mean <- Dh.gl.rmDupCoords@other$ind.metrics %>% 
  group_by(Pop = Buffer_15km) %>% 
  summarise(Site.lon = mean(lon), Site.lat = mean(lat))

# Calculate population average ancestry proportions and create a df with population coordinates
# Loop through for 1 to K number of clusters
coord.pie.Tess <- coord.pop.mean

for (j in 1:K) {
  coord.pie.col <- Dh.gl.rmDupCoords@other$ind.metrics %>%
    group_by(Pop = Buffer_15km) %>%
    summarise(Mean_AP = mean(get(paste0("TessClust", j))))
  
  coord.pie.Tess <- left_join(coord.pie.Tess, coord.pie.col, by = "Pop")
  colnames(coord.pie.Tess)[ncol(coord.pie.Tess)] <- paste0("ZZ_Mean_AP", j) 
  # Note: Adding "ZZ_" for plotting. 
  # Colour palette alphabetic
  # Adding to plot where map regions and pie charts use the same palette/fill term
}

# Tess pie charts of ancestry proportions on map
Tess.Pie.Map <- Map.plot +
  scale_fill_manual(values = c(MapPal, Pal[c(1,19)])) +
  geom_scatterpie(data = coord.pie.Tess, aes(x=Site.lon, 
                                         y=Site.lat, 
                                         group = Pop), 
                  cols = colnames(coord.pie.Tess[ , 4:ncol(coord.pie.Tess)]),  
                  pie_scale = 1.5, size = 0.1)

Tess.Combined <- ggarrange(TessBarPlot, Tess.Pie.Map, ncol = 1, heights = c(1/5, 4/5), labels = c("a)", "b)"), font.label = list(size = 8, face = "bold"))

# Check that figures are okay for colourblind
cvd_grid(Tess.Combined)

# Add simplified CV plot to figure
vp.TessCV <- viewport(width = 0.3, height = 0.225, x = 0.25, y = 0.228)

# Output to pdf 
pdf(paste0(Out.Path, "TessPieBar_K", K,".", sp, ".pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(Tess.Combined)
print(CrossValPlot.Tess.MP, vp = vp.TessCV)
dev.off()



#################
# LEA BAR PLOT #
#################

# Define value for K 
K <- 2

# get the cross-entropy of all runs for the chosen value of K
ce.LEA <- cross.entropy(snmf.obj, K = K)

# select the run with the lowest cross-entropy for the chosen value of K
best.LEA <- which.min(ce.LEA)

# save the Q-matrix for this run (this will be used for graphical presentation)
Q.matrix.LEA <- as.qmatrix(Q(object = snmf.obj, K = K, run = best.LEA))

# Reformat data for plotting
# Create a new column by combining longs and lats 
# Order by this column so bar chart matches up with map
Q.matrix.LEA.df <- cbind(sampleID = Dh.gl.rmDupCoords@other$ind.metrics$id, 
                          lon = Dh.gl.rmDupCoords@other$latlong$lon, 
                          lat = Dh.gl.rmDupCoords@other$latlong$lat, 
                          as.data.frame(Q.matrix.LEA[,]))

# Assign individuals to clusters based max ancestry proportion
# If ancestry proportion <0.7 mark as admixed
LEAAssign <- Q.matrix.LEA.df
# I want to name cluster 1 as the western-most population 
# Not necessary, but I prefer it aesthetically
# Check which column a western sample assigns to
colnames(LEAAssign[order(LEAAssign$lon)[1],c(4:(4+K-1))])[apply(LEAAssign[order(LEAAssign$lon)[1],c(4:(4+K-1))],1,which.max)]
# "V2" - therefore "V2" will be labelled "Cluster 1" and "V1" will be "Cluster 2" (i.e. reverse the order below)
colnames(LEAAssign)[c(4:(ncol(LEAAssign)))] <- paste0("LEAClust", rev(1:K))
LEA.Admixed <- apply(LEAAssign[,c(4:ncol(LEAAssign))],1,max) < 0.7
LEAAssign$LEAAdmixed <- ifelse(LEA.Admixed == TRUE, "Y", "N")
LEAAssign$LEAGenCluster <- colnames(LEAAssign[,c(4:(4+K-1))])[apply(LEAAssign[,c(4:(4+K-1))],1,which.max)]


LEAAssign <- left_join(LEAAssign,
                        LEAAssign %>%
                          group_by(LEAGenCluster) %>%
                          summarise(MeanLon = mean(lon)),
                        by = "LEAGenCluster")

# Now sort on coordinates
Qmat.K.sorted <-  LEAAssign[order(LEAAssign$MeanLon, 
                                   LEAAssign[,unique(LEAAssign$LEAGenCluster[LEAAssign$MeanLon == max(LEAAssign$MeanLon)])]), 
                             -which(colnames(LEAAssign) %in% c("MeanLon", "lon", "lat"))]

# Add new dummy IDs that are in the correct order for plotting
# Adjust maximum number/numbr of zeros acordsingly depending on how many samples
Qmat.K.sorted$PlotID <- c(paste0("Dummy_00", 1:9), 
                          paste0("Dummy_0", 10:99),
                          paste0("Dummy_", 100:nrow(Qmat.K.sorted)))

# Reformat for bar chart
Qmat.K <- tidyr::gather(Qmat.K.sorted, 
                        "popGroup", 
                        "prop",
                        -sampleID,
                        -PlotID,
                        -LEAAdmixed,
                        -LEAGenCluster)


# Save in workspace for combined plot later on
LEABarPlot <- ggplot(Qmat.K, aes(factor(PlotID), prop, fill = factor(popGroup))) +
  scale_fill_manual(values = Pal[c(1, 19)]) +
  geom_col(color = "transparent", size = 0.1, width = 1) +
  labs(y = " ") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "plain", size = 6, angle = 90),
        strip.background = element_rect(linetype = "solid", size = 0.1, 
                                        colour = "black", fill = "lightgrey"),
        axis.text = element_blank(),
        axis.title.y = element_text(size = 9, margin = margin(0,24,0,0), face ="bold"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Add LEA info into ind.metric slot of genlight
Dh.gl.rmDupCoords@other$ind.metrics <- cbind(Dh.gl.rmDupCoords@other$ind.metrics,
                                             LEAAssign[match(Dh.gl.rmDupCoords@other$ind.metrics$id, LEAAssign$sampleID), 
                                                        -which(colnames(LEAAssign) %in% c("sampleID", "lon", "lat", "MeanLon"))])

######################
# LEA PIE CHART MAP #
######################

# Calculate population average ancestry proportions and create a df with population coordinates
# Loop through for 1 to K number of clusters
coord.pie.LEA <- coord.pop.mean

for (j in 1:K) {
  coord.pie.col <- Dh.gl.rmDupCoords@other$ind.metrics %>%
    group_by(Pop = Buffer_15km) %>%
    summarise(Mean_AP = mean(get(paste0("LEAClust", j))))
  
  coord.pie.LEA <- left_join(coord.pie.LEA, coord.pie.col, by = "Pop")
  colnames(coord.pie.LEA)[ncol(coord.pie.LEA)] <- paste0("ZZ_Mean_AP", j) 
  # Note: Adding "ZZ_" for plotting. 
  # Colour palette alphabetic
  # Adding to plot where map regions and pie charts use the same palette/fill term
}

# LEA pie charts of ancestry proportions on map
LEA.Pie.Map <- Map.plot +
  scale_fill_manual(values = c(MapPal, Pal[c(1,19)])) +
  geom_scatterpie(data = coord.pie.LEA, aes(x=Site.lon, 
                                             y=Site.lat, 
                                             group = Pop), 
                  cols = colnames(coord.pie.LEA[ , 4:ncol(coord.pie.LEA)]),  
                  pie_scale = 1.5, size = 0.1)

LEA.Combined <- ggarrange(LEABarPlot, LEA.Pie.Map, ncol = 1, heights = c(1/5, 4/5), labels = c("a)", "b)"), font.label = list(size = 8, face = "bold"))


# Output to pdf
pdf(paste0(Out.Path, "LEAPieBar_K", K,".", sp, ".pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(LEA.Combined)
dev.off()


#######################
# RUN DAPC (Adegenet) #
#######################

# Transform data using PCA and run K-means clustering to determine the number of groups (K)
# Can keep all PCs for this step, so have just put 1000 (which should be way more than the max)
set.seed(123)
Clusters <- find.clusters(Dh.gl.rmDupCoords, n.iter = 1000, 
                          max.n.clust = Kmax, 
                          n.pca = 1000)

# Save as df for plotting
DAPC.BIC <- data.frame(Clusters = 1:Kmax, BIC = Clusters$Kstat)

K <- as.numeric(str_remove(names(Clusters$stat), "K="))

# Save cross entropy plot in workspace to plot later
DAPC.BIC.plot <- ggplot(DAPC.BIC, aes(x=Clusters, y=BIC)) + 
  geom_point(size = 0.8) +
  geom_line(size = 0.2) +
  labs(x = "K", y = "BIC") +
  scale_x_continuous(breaks = seq(1,10,1)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5),
        axis.title.y = element_text(size = 8, margin = margin(0,10,0,0), face ="bold"),
        axis.title.x = element_text(size = 8, margin = margin(10,0,0,0), face ="bold"),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none")

# Save as Pdf
ggsave(filename = paste0("DAPC.BIC_", sp, ".pdf"),
       plot = DAPC.BIC.plot,
       device = "pdf",
       path = Out.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# Run Discriminate Analysis (of PCA axes)
# Aim to retain a few PCs without sacrificing too much information. I chose to retain 100 PCs which accounted for ~70% of the cumulative genetic variance. 
# Next, to choose the number of discriminant functions to retain, interrogate the barplot of eigenvalues. For a small number of clusters, all eigenvalues can be retained since it won't be too difficult to visualise all discriminant functions. I chose to retain the maximum, which was 1. 
dapc1 <- dapc(x = Dh.gl.rmDupCoords, pop = Clusters$grp)


# Validate DAPC
# Cross validation:
# In cross-validation, the data is divided into two sets: a training set (typically comprising 90% of the data) and a validation set (which contains the remainder - by default, 10% of the data). The function xvalDapc uses stratified random sampling to select the validation set. This ensures that at least one member of each group or population in the original data is represented in both training and validation sets.

# DAPC is carried out on the training set with variable numbers of PCs retained, and the degree to which the analysis is able to accurately predict the group membership of excluded (masked) individuals (i.e. the validation set) is used to identify the optimal number of PCs to retain. At each level of PC retention, the sampling and DAPC procedures are repeated n.rep times. (By default = 30 replicates).

# We first need to remove any missing data (and replace it with the mean genotype)
mat <- tab(Dh.gl.rmDupCoords, NA.method="mean")
grp <- pop(Dh.gl.rmDupCoords)

# Run cross validation:
set.seed(123)
xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)

# Look at the results:
xval[2:6]
crossVal <- xval$`Cross-Validation Results`
crossval.best <- paste0("Highest Mean Success = ", as.numeric(xval$`Number of PCs Achieving Highest Mean Success`), " PCs\nLowest RMSE = ", as.numeric(xval$`Number of PCs Achieving Lowest MSE`), " PCs")

DAPC.crossVal <- ggplot(crossVal, aes(x= n.pca, y= success)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(direction=1) +
  geom_point(size = 0.3) +
  labs(y = "Proportion of successful predictions", x = "Number of PCs") +
  xlim(c(0, max(crossVal$n.pca)+10)) +
  ylim(c(0, max(crossVal$success)+0.1)) +
  geom_text(label = crossval.best, size = 1.5, x = 8, y = 0.3, hjust = 0) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 5, face ="plain"),
        axis.title.x = element_text(size = 6, face ="bold", margin = margin(10,0,0,0)),
        axis.title.y = element_text(size = 6, face ="bold", margin = margin(0,10,0,0)),
        axis.ticks = element_line(size = 0.2),
        legend.position = "none")
DAPC.crossVal

# Save as Pdf
ggsave(filename = paste0("CrossVal.DAPC.", sp, ".pdf"),
       plot = DAPC.crossVal,
       device = "pdf",
       path = Out.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)



# Decide on optimal number of PCs to keep based on cross validation
OptimalPC <- 40
  
# Re-Run DAPC, retaining the appropriate number of PCs
dapc2 <- dapc(x = Dh.gl.rmDupCoords, pop = Clusters$grp, n.pca = OptimalPC)
# retained 1 DF

# Retrieve dapc proportions 
dapc.post <- dapc2$posterior

# Reformat data for plotting
# Create a new column by combining longs and lats 
# Order by this column so bar chart matches up with map
dapc.post.df <- cbind(sampleID = Dh.gl.rmDupCoords$other$ind.metrics$id, 
                        lon = Dh.gl.rmDupCoords$other$latlong$lon, 
                        lat = Dh.gl.rmDupCoords$other$latlong$lat, 
                        as.data.frame(dapc.post[,]))

# Assign individuals to clusters based max ancestry proportion
# If ancestry proportion <0.7 mark as admixed
dapcAssign <- dapc.post.df
# I want to name cluster 1 as the western-most population 
# Not necessary, but I prefer it aesthetically
# Check which column a western sample assigns to
colnames(dapcAssign[order(dapcAssign$lon)[1],c(4:(4+K-1))])[apply(dapcAssign[order(dapcAssign$lon)[1],c(4:(4+K-1))],1,which.max)]
# "2" - therefore "2" will be labelled "Cluster 1" and "1" will be "Cluster 2" (i.e. reverse the order below)
colnames(dapcAssign)[c(4:(ncol(dapcAssign)))] <- paste0("dapcClust", rev(1:K))
dapc.Admixed <- apply(dapcAssign[,c(4:ncol(dapcAssign))],1,max) < 0.7
dapcAssign$dapcAdmixed <- ifelse(dapc.Admixed == TRUE, "Y", "N")
dapcAssign$dapcGenCluster <- colnames(dapcAssign[,c(4:(4+K-1))])[apply(dapcAssign[,c(4:(4+K-1))],1,which.max)]

dapcAssign <- left_join(dapcAssign,
                        dapcAssign %>%
                          group_by(dapcGenCluster) %>%
                          summarise(MeanLon = mean(lon)),
                         by = "dapcGenCluster")
  
# Now sort on coordinates
post.K.sorted <-  dapcAssign[order(dapcAssign$MeanLon, 
                                    dapcAssign[,unique(dapcAssign$dapcGenCluster[dapcAssign$MeanLon == max(dapcAssign$MeanLon)])]),
                              -which(colnames(dapcAssign) %in% c("MeanLon", "lon", "lat"))]
  
# Add new dummy IDs that are in the correct order for plotting
post.K.sorted$PlotID <- c(paste0("Dummy_00", 1:9), 
                              paste0("Dummy_0", 10:99),
                              paste0("Dummy_", 100:nrow(post.K.sorted)))

# Reformat for bar chart
dapc.K <- tidyr::gather(post.K.sorted, 
                        "popGroup", 
                        "prop",
                        -sampleID,
                        -PlotID,
                        -dapcAdmixed,
                        -dapcGenCluster)
  
# Save in workspace for combined plot later on
dapcBarPlot <- ggplot(dapc.K, aes(factor(PlotID), prop, fill = factor(popGroup))) +
    scale_fill_manual(values = Pal[c(1,19)]) +
    geom_col(color = "transparent", size = 0.1, width = 1) +
    labs(y = " ") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", 
                                      fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          strip.text = element_text(face = "plain", size = 6, angle = 90),
          strip.background = element_rect(linetype = "solid", size = 0.1, 
                                          colour = "black", fill = "lightgrey"),
          axis.text = element_blank(),
          axis.title.y = element_text(size = 9, margin = margin(0,24,0,0), face ="bold"),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
  
# Add dapc info into ind.metric slot of genlight
Dh.gl.rmDupCoords@other$ind.metrics <- cbind(Dh.gl.rmDupCoords@other$ind.metrics,
                                             dapcAssign[match(Dh.gl.rmDupCoords@other$ind.metrics$id, dapcAssign$sampleID), 
                                                       -which(colnames(dapcAssign) %in% c("sampleID", "lon", "lat", "MeanLon"))])


######################
# DAPC PIE CHART MAP #
######################

# Calculate population average ancestry proportions and create a df with population coordinates
# Loop through for 1 to K number of clusters
coord.pie.dapc <- coord.pop.mean

for (j in 1:K) {
  coord.pie.col <- Dh.gl.rmDupCoords@other$ind.metrics %>%
    group_by(Pop = Buffer_15km) %>%
    summarise(Mean_AP = mean(get(paste0("dapcClust", j))))
  
  coord.pie.dapc <- left_join(coord.pie.dapc, coord.pie.col, by = "Pop")
  colnames(coord.pie.dapc)[ncol(coord.pie.dapc)] <- paste0("ZZ_Mean_AP", j) 
  # Note: Adding "ZZ_" for plotting. 
  # Colour palette alphabetic
  # Adding to plot where map regions and pie charts use the same palette/fill term
}

# dapc pie charts of ancestry proportions on map
dapc.Pie.Map <- Map.plot +
  scale_fill_manual(values = c(MapPal, Pal[c(1,19)])) +
  geom_scatterpie(data = coord.pie.dapc, aes(x=Site.lon, 
                                            y=Site.lat, 
                                            group = Pop), 
                  cols = colnames(coord.pie.dapc[ , 4:ncol(coord.pie.dapc)]),  
                  pie_scale = 1.5, size = 0.1)

dapc.Combined <- ggarrange(dapcBarPlot, dapc.Pie.Map, ncol = 1, heights = c(1/5, 4/5), labels = c("a)", "b)"), font.label = list(size = 8, face = "bold"))


# Output to pdf
pdf(paste0(Out.Path, "dapcPieBar_K", K,".", sp, ".pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(dapc.Combined)
dev.off()


###################################################
# Output Genlight with updated individual metrics #
#             (i.e. clustering info)              #
###################################################

save(Dh.gl.rmDupCoords, file = "IBB/IBB_outputs/Cleaned.Unrelated.gl.Dh_ClustInfo.rdata")
