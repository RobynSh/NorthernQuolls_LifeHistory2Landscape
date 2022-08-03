library(rgdal)
library(raster)
library(dplyr)
library(ggplot2)
library(ggsn)
library(dartR)
library(ggpubr)
library(colorblindr)
library(tess3r)
library(LEA)
library(scatterpie)
library(adegenet)
library(viridis)

################################################
#      SET VARIABLEs AND LOAD IN DATA          #
################################################

# Species name
sp <- "Dh"

# Paths
Work.Path <- paste0("/Users/robynshaw/Google Drive/Work/Pilbara_Small_Mammals/Manuscripts/Quoll_Paper/Full_Analysis/IBB/")
Fig.Path <- paste0("/Users/robynshaw/Google Drive/Work/Pilbara_Small_Mammals/Manuscripts/Quoll_Paper/Full_Analysis/IBB/Figures/")
gl.Path <- paste0("/Users/robynshaw/Google Drive/Work/Pilbara_Small_Mammals/Manuscripts/Quoll_Paper/Full_Analysis/Final_Datasets/")
Map.Path <- "/Users/robynshaw/Google Drive/Work/Pilbara_Small_Mammals/GIS/Pilbara GIS layers/"

# Set working directory
setwd(Work.Path)


################
# PREPARE MAPS #
################

# Load in polygon of IBRA subregions
StudyMap.Pilb <- readOGR(paste0(Map.Path, "PilbaraIBRA.shp"))
StudyMap.WA <- readOGR(paste0(Map.Path, "IBRA_regions_WA.shp"))

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

# Save study map as an r object so I can load it in when needed
save(Map.plot, file= "Quoll.ggplot.StudyMap.rdata")

# This script takes a genlight as the input format
# I output my genlight as an R object after filtering loci and
# removing duplicated and related individuals
# The corresponding meta-data is saved in the "other" slot

load(paste0(gl.Path, "Cleaned.Unrelated.gl.Dh.rdata"))

###############
# SUBSET DATA #
###############

Dh.gl.rmRel <- Dh.gl

# Read in metadata
Dh.MD <- read.csv(paste0(gl.Path, "Cleaned.Unrelated.Ind.metadata.Dh_New.csv"))
# Add in pop descriptions (but first check datasets match)
Dh.gl.rmRel@other$ind.metrics <- Dh.MD[match(Dh.gl.rmRel@ind.names, Dh.MD$id), ]

BuffVec.order <- cbind(unique(Dh.MD$Buffer.15km[order(Dh.MD$lon, Dh.MD$lat)]), c(paste0("B15.0", 1:9), paste0("B15.", 10:length(unique(Dh.MD$Buffer.15km[order(Dh.MD$lon, Dh.MD$lat)])))))

Dh.gl.rmRel@pop <- as.factor(BuffVec.order[,2][match(Dh.gl.rmRel@other$ind.metrics$Buffer.15km, BuffVec.order[,1])])

levels(Dh.gl.rmRel@pop) <- c(levels(Dh.gl.rmRel@pop), "B15.00")
Dh.gl.rmRel@pop[Dh.gl.rmRel@other$ind.metrics$pop == "Dolphin Island"] <- "B15.00"

# Subset dataset so that there is only one sample per location
Dh.gl.rmDupCoords <- Dh.gl.rmRel[!duplicated(Dh.gl.rmRel@other$latlong), ]

# Create palette and shape vectors
Pal <- c(viridis(n = 11, option = "D", begin = 0.2), viridis(n = 10, option = "A", direction = -1, begin = 0.5))  #c("#b22e21", "#c76c19", "#ffbd5c", "#00858e", "#6ba9e3")
Shape <- rep(21:25, 5) #c(21:25)

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


# Check that figures are okay for colourblind
cvd_grid(Samp.Map)
cvd_grid(Pcoa_1V2)
cvd_grid(Pcoa_1V3)
cvd_grid(Pcoa_2V3)

# Combine plots

PCoA.Eig <- ggarrange(Samp.Map, ggarrange(Pcoa_1V2, Pcoa_1V3, ncol = 2), ncol = 1, heights = c(5/8,3/8), labels = c("a)", "b)"), font.label = list(size = 8, face = "plain"))


#PCoA.Eig <- ggarrange(Samp.Map, ggarrange(Pcoa_1V2, ggarrange(Pcoa_1V3, Pcoa_2V3, ncol = 1), widths = c(3/5, 2/5)), ncol = 1, labels = c("a)", "b)"), font.label = list(size = 8, face = "plain"))

vp <- viewport(width = 0.178, height = 0.12, x = 0.395, y = 0.3)

# Output to pdf and jpeg

pdf(paste0(Fig.Path, sp, "_Map.PCoA.pdf"), width = 7.1, height = 6.5, useDingbats = FALSE)
print(PCoA.Eig)
print(EigPlot, vp = vp)
dev.off()

jpeg(paste0(Fig.Path,sp, "_Map.PCoA.jpg"), width = 7.1, height = 8, units = "in", res = 300)
print(PCoA.Eig)
print(EigPlot, vp = vp)
dev.off()



##########################
#  Run Tess3R for Total #
##########################

# Define maximum value of K to test
# From the PCoA it doesn't look like there could be much more than 4 groups, so I'll test a maximum of 7
Kmax <- 7

Geno.Tess3R <- as.matrix(Dh.gl.rmDupCoords)
Coords.Tess3R <- as.matrix(Dh.gl.rmDupCoords@other$latlong[, c(2,1)])
Pops.Tess3R <- as.character(Dh.gl.rmDupCoords@pop)


## RUN TESS3R ##
Tess.obj <- tess3(X = Geno.Tess3R, coord = Coords.Tess3R, K = 1:Kmax, rep = 100, method = "projected.ls", ploidy = 2, openMP.core.num = 4, max.iteration = 1000, keep = "best", mask = 0.1) 

# Save Tess object
save(Tess.obj, file = paste0(Work.Path, "Tess.obj.FullDat.best100reps.K1_", Kmax, ".", sp, ".R"))

# Load in tess object if replotting
# load(file = paste0(Work.Path, "Tess.obj.FullDat.best100reps.K1_", Kmax, ".", sp, ".R"))


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
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# Save as jpeg
ggsave(filename = paste0("CrossVal.Tess.", sp, ".jpg"),
       plot = CrossValPlot.Tess,
       device = "jpeg",
       path = Fig.Path,
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
write.geno(R = Geno.Tess3R, output.file = paste0(Work.Path, sp, ".geno"))

# Run LEA Algorithm
snmf.obj <- snmf(input.file = paste0(Work.Path, sp, ".geno"), K=1:Kmax, alpha = 10, tolerance = 0.00001, percentage = 0.05, iterations = 1000, ploidy=2, entropy=T, project="new", repetitions = 100) 
  
# Save snmf object
save(snmf.obj, file = paste0(Work.Path, "snmf.obj.FullDat.best100reps.K1_", Kmax, ".", sp, ".R"))

# Load snmf object if just redoing plots
# load(file = paste0(Work.Path, "snmf.obj.FullDat.best100reps.K1_", Kmax, ".", sp, ".R"))

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
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# Save as jpeg
ggsave(filename = paste0("CrossVal.LEA.", sp, ".jpg"),
       plot = CrossValPlot.LEA,
       device = "jpeg",
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)


CrossValPlot.LEA.VP <- ggplot(LEA.CrossEnt, aes(x=PlotOrder, y=min)) + 
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
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

#################
# Tess BAR PLOT #
#################

# Define value for K 
K <- 2

# Retrieve tess3 Q matrix rep with the lowest cross-entropy criterion
Q.matrix.Tess <- qmatrix(Tess.obj, K = K)
  
# Create palette
TessPal <- Pal[c(1,19)]

# Reformat data for plotting
Q.matrix.Tess.df <- cbind(Dh.gl.rmDupCoords@ind.names, Dh.gl.rmDupCoords@other$ind.metrics$Pop_Descriptive, as.data.frame(Q.matrix.Tess[,]))
# Rename columns
# NOTE: Switch around clusters so they match the Tess analysis (sometimes the order changes just because of stochastic reasons)
colnames(Q.matrix.Tess.df) <- c("sampleID", "pop", "Cluster2", "Cluster1")

# Sort by cluster

# Assign individuals to clusters
Q.matrix.Tess.df$GenCluster <- colnames(Q.matrix.Tess.df[, 3:ncol(Q.matrix.Tess.df)])[apply(Q.matrix.Tess.df[, 3:ncol(Q.matrix.Tess.df)],1,which.max)]

# Create a column that specifies if the sample was assigned to the cluster or not (i.e. the q score was less than 0.7)
Q.matrix.Tess.df$Proportion <- ifelse(Q.matrix.Tess.df[, 3] >=0.7 | Q.matrix.Tess.df[, (ncol(Q.matrix.Tess.df)-1)] >=0.7, "Assigned", "Admixed")


# Sort by admixture proportion and cluster
Qmat.sorted.Tess <-  Q.matrix.Tess.df[order(Q.matrix.Tess.df$Cluster1, decreasing = TRUE), -c(which(colnames(Q.matrix.Tess.df) == "GenCluster"), which(colnames(Q.matrix.Tess.df) == "Proportion"))]

# Add new dummy IDs that are in the correct order for plotting
Qmat.sorted.Tess$sampleID.ordered <- c(paste0("Dummy_00", 1:9), paste0("Dummy_0", 10:99), paste0("Dummy_", 100:length(Dh.gl.rmDupCoords@ind.names)))

# Convert to long format for plotting
Qmat.gather.Tess <- tidyr::gather(Qmat.sorted.Tess, "popGroup", "prob", -sampleID, -sampleID.ordered, -pop)


# Plot
TessBarPlot <- ggplot(Qmat.gather.Tess, aes(factor(sampleID.ordered), prob, fill = factor(popGroup))) +
  scale_fill_manual(values = TessPal) +
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

######################
# Tess PIE CHART MAP #
######################

# Get mean of individual coordinates so we have a single set of coordinates for each pop
Tess.order <- match(Qmat.sorted.Tess$sampleID, Dh.gl.rmDupCoords@other$ind.metrics$id)
Q.matrix.Tess.df.pie <- cbind(Qmat.sorted.Tess, Dh.gl.rmDupCoords@other$latlong[Tess.order, ], Dh.gl.rmDupCoords@pop[Tess.order])
colnames(Q.matrix.Tess.df.pie)[8] <- c("Buffer.15km")
coord.pop.mean.Tess <- Q.matrix.Tess.df.pie %>% group_by(Buffer.15km) %>% summarise(Site.lon = mean(lon), Site.lat = mean(lat))

# Calculate population average ancestry proportions and create an array with population coordinates:
Q.matrix.Tess.df.pop <- Q.matrix.Tess.df.pie %>%
  group_by(Buffer.15km) %>%
  summarise(Mean_AP1 = mean(Cluster1), Mean_AP2 = mean(Cluster2))    


# Reformat coordinates and q matrix for plotting
coord.pie.Tess <- left_join(coord.pop.mean.Tess, Q.matrix.Tess.df.pop, by = "Buffer.15km")
coord.pie.Tess$Clust <- ifelse(coord.pie.Tess$Mean_AP1 >= 1/2, "Cluster1", "Cluster2")
coord.pie.Tess <- coord.pie.Tess[order(coord.pie.Tess$Clust), ]
coord.pie.Tess$Group <- c(paste0("Group.0", 1:9), paste0("Group.", 10:nrow(coord.pie.Tess)))


# Tess pie charts of ancestry proportions on map
Tess.Pie.Map <- Map.plot +
  scale_fill_manual(values = c(MapPal, TessPal)) +
  geom_scatterpie(data = coord.pie.Tess, aes(x=Site.lon, y=Site.lat, group = Group), 
                  cols = colnames(coord.pie.Tess[ , 4:(ncol(coord.pie.Tess) - 2)]),  
                  pie_scale = 1.5, size = 0.1)

Tess.Combined <- ggarrange(TessBarPlot, Tess.Pie.Map, ncol = 1, heights = c(1/5, 4/5), labels = c("a)", "b)"), font.label = list(size = 8, face = "bold"))

# Check that figures are okay for colourblind
cvd_grid(Tess.Combined)

# Add simplified CV plot to figure
vp.TessCV <- viewport(width = 0.3, height = 0.225, x = 0.25, y = 0.228)

# Output to pdf and jpeg

pdf(paste0(Fig.Path, "TessPieBar_K", K,".", sp, ".pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(Tess.Combined)
print(CrossValPlot.Tess.MP, vp = vp.TessCV)
dev.off()

jpeg(paste0(Fig.Path, "TessPieBar_K", K,".", sp, ".jpg"), width = 6.8, height = 4.74, units = "in", res = 300)
print(Tess.Combined)
print(CrossValPlot.Tess.MP, vp = vp.TessCV)
dev.off()


pdf(paste0(Fig.Path, "TessPieBar_K", K,".", sp, "_noVP.pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(Tess.Combined)
dev.off()

jpeg(paste0(Fig.Path, "TessPieBar_K", K,".", sp, "_noVP.jpg"), width = 6.8, height = 4.74, units = "in", res = 300)
print(Tess.Combined)
dev.off()


#################
# LEA BAR PLOT #
#################

# List value for K
# Choose by looking at cross validation plots
K <- 2

# get the cross-entropy of all runs for the chosen value of K
ce.LEA <- cross.entropy(snmf.obj, K = K)
  
# select the run with the lowest cross-entropy for the chosen value of K
best.LEA <- which.min(ce.LEA)
  
# save the Q-matrix for this run (this will be used for graphical presentation)
Q.matrix.LEA <- as.qmatrix(Q(object = snmf.obj, K = K, run = best.LEA))
  
# Create palette (rearrange for different spp so clusters are coloured consistently by location)
LEAPal <- TessPal

# Reformat data for plotting
Q.matrix.LEA.df <- cbind(Dh.gl.rmDupCoords@ind.names, Dh.gl.rmDupCoords@other$ind.metrics$Pop_Descriptive, as.data.frame(Q.matrix.LEA[,]))
colnames(Q.matrix.LEA.df) <- c("sampleID", "pop", "Cluster2", "Cluster1")

# Sort by cluster

# Assign individuals to clusters
Q.matrix.LEA.df$GenCluster <- colnames(Q.matrix.LEA.df[, 3:ncol(Q.matrix.LEA.df)])[apply(Q.matrix.LEA.df[, 3:ncol(Q.matrix.LEA.df)],1,which.max)]

# Create a column that specifies if the sample was assigned to the cluster or not (i.e. the q score was less than 0.7)
Q.matrix.LEA.df$Proportion <- ifelse(Q.matrix.LEA.df[, 3] >=0.7 | Q.matrix.LEA.df[, (ncol(Q.matrix.LEA.df)-1)] >=0.7, "Assigned", "Admixed")


# Sort by admixture proportion and cluster
Qmat.sorted.LEA <-  Q.matrix.LEA.df[order(Q.matrix.LEA.df$Cluster1, decreasing = TRUE), -c(which(colnames(Q.matrix.LEA.df) == "GenCluster"), which(colnames(Q.matrix.LEA.df) == "Proportion"))]

# Add new dummy IDs that are in the correct order for plotting
Qmat.sorted.LEA$sampleID.ordered <- c(paste0("Dummy_00", 1:9), paste0("Dummy_0", 10:99), paste0("Dummy_", 100:length(Dh.gl.rmDupCoords@ind.names)))

# Convert to long format for plotting
Qmat.gather.LEA <- tidyr::gather(Qmat.sorted.LEA, "popGroup", "prob", -sampleID, -sampleID.ordered, -pop)


# Plot
LEABarPlot <- ggplot(Qmat.gather.LEA, aes(factor(sampleID.ordered), prob, fill = factor(popGroup))) +
  scale_fill_manual(values = LEAPal) +
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
  
######################
# LEA PIE CHART MAP #
######################

# Get mean of individual coordinates so we have a single set of coordinates for each pop
LEA.order <- match(Qmat.sorted.LEA$sampleID, Dh.gl.rmDupCoords@other$ind.metrics$id)
Q.matrix.LEA.df.pie <- cbind(Qmat.sorted.LEA, Dh.gl.rmDupCoords@other$latlong[LEA.order, ], Dh.gl.rmDupCoords@pop[LEA.order])
colnames(Q.matrix.LEA.df.pie)[8] <- c("Buffer.15km")
coord.pop.mean.LEA <- Q.matrix.LEA.df.pie %>% group_by(Buffer.15km) %>% summarise(Site.lon = mean(lon), Site.lat = mean(lat))

# Calculate population average ancestry proportions and create an array with population coordinates:
Q.matrix.LEA.df.pop <- Q.matrix.LEA.df.pie %>%
  group_by(Buffer.15km) %>%
  summarise(Mean_AP1 = mean(Cluster1), Mean_AP2 = mean(Cluster2))    


# Reformat coordinates and q matrix for plotting
coord.pie.LEA <- left_join(coord.pop.mean.LEA, Q.matrix.LEA.df.pop, by = "Buffer.15km")
coord.pie.LEA$Clust <- ifelse(coord.pie.LEA$Mean_AP1 >= 1/2, "Cluster1", "Cluster2")
coord.pie.LEA <- coord.pie.LEA[order(coord.pie.LEA$Clust), ]
coord.pie.LEA$Group <- c(paste0("Group.0", 1:9), paste0("Group.", 10:nrow(coord.pie.LEA)))


# LEA pie charts of ancestry proportions on map
LEA.Pie.Map <- Map.plot +
  scale_fill_manual(values = c(MapPal, LEAPal)) +
  geom_scatterpie(data = coord.pie.LEA, aes(x=Site.lon, y=Site.lat, group = Group), 
                  cols = colnames(coord.pie.LEA[ , 4:(ncol(coord.pie.LEA) - 2)]),  
                  pie_scale = 1.5, size = 0.1)

LEA.Combined <- ggarrange(LEABarPlot, LEA.Pie.Map, ncol = 1, heights = c(1/5, 4/5), labels = c("a)", "b)"), font.label = list(size = 8, face = "bold"))

# Check that figures are okay for colourblind
cvd_grid(LEA.Combined)

# Add simplified CV plot to figure
vp.LEACV <- viewport(width = 0.3, height = 0.225, x = 0.25, y = 0.228)

# Output to pdf and jpeg

pdf(paste0(Fig.Path, "LEAPieBar_K", K,".", sp, ".pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(LEA.Combined)
print(CrossValPlot.LEA.MP, vp = vp.LEACV)
dev.off()

jpeg(paste0(Fig.Path, "LEAPieBar_K", K,".", sp, ".jpg"), width = 6.8, height = 4.74, units = "in", res = 300)
print(LEA.Combined)
print(CrossValPlot.LEA.MP, vp = vp.LEACV)
dev.off()

pdf(paste0(Fig.Path, "LEAPieBar_K", K,".", sp, "_noVP.pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(LEA.Combined)
dev.off()

jpeg(paste0(Fig.Path, "LEAPieBar_K", K,".", sp, "_noVP.jpg"), width = 6.8, height = 4.74, units = "in", res = 300)
print(LEA.Combined)
dev.off()



#######################
# RUN DAPC (Adegenet) #
#######################

# Step 1: Transform data using PCA and run K-means clustering to determine the number of groups (K)
# Keep ~all PCs for this step (200)
Clusters <- find.clusters(Dh.gl.rmDupCoords, max.n.clust= Kmax, n.iter = 1000, n.pca = 200)

# Choose number of clusters
K <- 2

# Plot BIC so I can save it
DAPC.BIC.df <- data.frame(K = 1:Kmax, BIC = Clusters$Kstat)

DAPC.BIC.plot <- ggplot(DAPC.BIC.df, aes(x=K, y=BIC)) + 
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
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# Save as jpeg
ggsave(filename = paste0("DAPC.BIC_", sp, ".jpg"),
       plot = DAPC.BIC.plot,
       device = "jpeg",
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

DAPC.BIC.plot.VP <- ggplot(DAPC.BIC.df, aes(x=K, y=BIC)) + 
  geom_point(size = 0.8) +
  geom_line(size = 0.2) +
  labs(x = "K", y = "BIC") +
  scale_x_continuous(breaks = seq(1,10,1)) +
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

# Check how well descriptive pops align with the chosen value of K:
# Table:
table(pop(Dh.gl.rmDupCoords), Clusters$grp)


# Step 2: Run Discriminate Analysis (of PCA axes)
# Aim to retain a few PCs without sacrificing too much information. I chose to retain 100 PCs which accounted for ~70% of the cumulative genetic variance. 
# Next, to choose the number of discriminant functions to retain, interrogate the barplot of eigenvalues. For a small number of clusters, all eigenvalues can be retained since it won't be too difficult to visualise all discriminant functions. I chose to retain the maximum, which was 1. 
dapc1 <- dapc(x = Dh.gl.rmDupCoords, pop = Clusters$grp)


# Step 3: Validate DAPC

# Using the a-score
# The a-score measures the trade-off between power of discrimination and over-fitting using the difference between the proportion of successful reassignments of the analysis (observed discrimination) and values obtained using random groups (random discrimination). i.e. The proportion of successful reassignments corrected for the number of retained PCs. 

# The number of retained PCs can be chosen so as to optimise the a-score, using the function optim.a.score. To reduce computational time, this method evaluates the first lot of retained PCs (i.e. because the first few are likely to be the most important), and uses spline interpolation to approximate the optimal number of PCs to retain. Then, you can evaluate all solutions within a restrained range using the argument n.pca.
optim.a.score(dapc1, n.sim = 100)
optim.a.score(dapc1, n.pca = 1:30, n.sim = 100)

# This analysis suggests that keeping 1 PC is the best.


# Cross validation:
# In cross-validation, the data is divided into two sets: a training set (typically comprising 90% of the data) and a validation set (which contains the remainder - by default, 10% of the data). The function xvalDapc uses stratified random sampling to select the validation set. This ensures that at least one member of each group or population in the original data is represented in both training and validation sets.

# DAPC is carried out on the training set with variable numbers of PCs retained, and the degree to which the analysis is able to accurately predict the group membership of excluded (masked) individuals (i.e. the validation set) is used to identify the optimal number of PCs to retain. At each level of PC retention, the sampling and DAPC procedures are repeated n.rep times. (By default = 30 replicates).

# We first need to remove any missing data (and replace it with the mean genotype)
mat <- tab(Dh.gl.rmDupCoords, NA.method="mean")
grp <- pop(Dh.gl.rmDupCoords)

# Run cross validation:
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
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# Save as jpeg
ggsave(filename = paste0("CrossVal.DAPC.", sp, ".jpg"),
       plot = DAPC.crossVal,
       device = "jpeg",
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)


# Cross validation suggests that retaining 60 PCs results in the lowest RMSE and the highest predictive success


# Compare the recommendations for the number of PCs to retain based on the a-score and on cross validation. They give different numbers, so I'm going to stick with cross validation
OptimalPC <- 20

# Step 4: Re-Run Discriminate Analysis, retaining the appropriate number of PCs
dapc2 <- dapc(x = Dh.gl.rmDupCoords, pop = Clusters$grp, n.pca = OptimalPC)
# retained 1 discriminant function

# The principal components of a DAPC can be plotted to provide a direct visual assessment of between-group structures.


# look at Discriminant Function 1 of the DAPC:
DAPC.density <- as.data.frame(cbind(dapc2$ind.coord, dapc2$assign))
DAPC.density$Axes <- "DF1"
colnames(DAPC.density)[1:2] <- c("LD.1", "Cluster")
DAPC.density$Cluster[DAPC.density$Cluster == 1] <- 9
DAPC.density$Cluster[DAPC.density$Cluster == 2] <- 1
DAPC.density$Cluster[DAPC.density$Cluster == 9] <- 2

DAPC.density$Cluster <- as.factor(DAPC.density$Cluster)

# Create palette (rearrange for different spp so clusters are coloured consistently by location)
DAPCPal <- TessPal

DAPC.dens <- ggplot(DAPC.density, aes(x = LD.1*-1, fill = Cluster)) +
  geom_density(size = 0.2, alpha = 0.8) +
  scale_fill_manual(values = DAPCPal) +
  xlim(c(-9,5.5)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        strip.text = element_text(face = "bold", size = 7),
        strip.background = element_rect(linetype = "solid", size = 0.2, 
                                        colour = "black", fill = "transparent"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

DAPC.dens

# Save as Pdf
ggsave(filename = paste0("DF1.DAPC.", sp, ".pdf"),
       plot = DAPC.dens,
       device = "pdf",
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# Save as jpeg
ggsave(filename = paste0("DF1.DAPC.", sp, ".jpg"),
       plot = DAPC.dens,
       device = "jpeg",
       path = Fig.Path,
       scale = 1,
       width = 3.2,
       height = 2,
       units = "in",
       dpi = 300)

# Prepare dapc as dataframe comparable to Q.matrix.df for Tess/LEA:
DAPC.df <- data.frame(sampleID = rownames(dapc2$ind.coord), pop = Q.matrix.LEA.df$pop[match(rownames(dapc2$ind.coord), Q.matrix.LEA.df$sampleID)], Cluster1 = as.matrix(dapc2$posterior)[,2], Cluster2 = as.matrix(dapc2$posterior)[,1], GenCluster = paste0("Cluster", dapc2$assign))

# Create a column that specifies if the sample was assigned to the cluster or not (i.e. the score was less than 0.7)
DAPC.df$Proportion <- ifelse(DAPC.df[, 3] >=0.7 | DAPC.df[, (ncol(DAPC.df)-1)] >=0.7, "Assigned", "Admixed")

# Sort by admixture proportion and cluster
DAPC.df.sorted <-  DAPC.df[order(DAPC.df$Cluster1, decreasing = TRUE), -c(which(colnames(DAPC.df) == "GenCluster"), which(colnames(DAPC.df) == "Proportion"))]

# Add new dummy IDs that are in the correct order for plotting
DAPC.df.sorted$sampleID.ordered <- c(paste0("Dummy_00", 1:9), paste0("Dummy_0", 10:99), paste0("Dummy_", 100:length(Dh.gl.rmDupCoords@ind.names)))

# Convert to long format for plotting
DAPC.df.gather <- tidyr::gather(DAPC.df.sorted, "popGroup", "prob", -sampleID, -sampleID.ordered, -pop)


# Plot
DAPCBarPlot <- ggplot(DAPC.df.gather, aes(factor(sampleID.ordered), prob, fill = factor(popGroup))) +
  scale_fill_manual(values = DAPCPal) +
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


######################
# DAPC PIE CHART MAP #
######################
Qmat.sorted.DAPC <- DAPC.df.sorted

# Get mean of individual coordinates so we have a single set of coordinates for each pop
DAPC.order <- match(Qmat.sorted.DAPC$sampleID, Dh.gl.rmDupCoords@other$ind.metrics$id)
Q.matrix.DAPC.df.pie <- cbind(Qmat.sorted.DAPC, Dh.gl.rmDupCoords@other$latlong[DAPC.order, ], Dh.gl.rmDupCoords@pop[DAPC.order])
colnames(Q.matrix.DAPC.df.pie)[8] <- c("Buffer.15km")
coord.pop.mean.DAPC <- Q.matrix.DAPC.df.pie %>% group_by(Buffer.15km) %>% summarise(Site.lon = mean(lon), Site.lat = mean(lat))

# Calculate population average ancestry proportions and create an array with population coordinates:
Q.matrix.DAPC.df.pop <- Q.matrix.DAPC.df.pie %>%
  group_by(Buffer.15km) %>%
  summarise(Mean_AP1 = mean(Cluster1), Mean_AP2 = mean(Cluster2))    


# Reformat coordinates and q matrix for plotting
coord.pie.DAPC <- left_join(coord.pop.mean.DAPC, Q.matrix.DAPC.df.pop, by = "Buffer.15km")
coord.pie.DAPC$Clust <- ifelse(coord.pie.DAPC$Mean_AP1 >= 1/2, "Cluster1", "Cluster2")
coord.pie.DAPC <- coord.pie.DAPC[order(coord.pie.DAPC$Clust), ]
coord.pie.DAPC$Group <- c(paste0("Group.0", 1:9), paste0("Group.", 10:nrow(coord.pie.DAPC)))


# DAPC pie charts of ancestry proportions on map
DAPC.Pie.Map <- Map.plot +
  scale_fill_manual(values = c(MapPal, DAPCPal)) +
  geom_scatterpie(data = coord.pie.DAPC, aes(x=Site.lon, y=Site.lat, group = Group), 
                  cols = colnames(coord.pie.DAPC[ , 4:(ncol(coord.pie.DAPC) - 2)]),  
                  pie_scale = 1.5, size = 0.1)

DAPC.Combined <- ggarrange(DAPCBarPlot, DAPC.Pie.Map, ncol = 1, heights = c(1/5, 4/5), labels = c("a)", "b)"), font.label = list(size = 8, face = "bold"))

# Check that figures are okay for colourblind
cvd_grid(DAPC.Combined)

# Add simplified CV plot to figure
vp.DAPCBIC <- viewport(width = 0.3, height = 0.225, x = 0.25, y = 0.228)

# Output to pdf and jpeg

pdf(paste0(Fig.Path, "DAPCPieBar_K", K,".", sp, ".pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(DAPC.Combined)
print(DAPC.BIC.plot.VP, vp = vp.DAPCBIC)
dev.off()

jpeg(paste0(Fig.Path, "DAPCPieBar_K", K,".", sp, ".jpg"), width = 6.8, height = 4.74, units = "in", res = 300)
print(DAPC.Combined)
print(DAPC.BIC.plot.VP, vp = vp.DAPCBIC)
dev.off()

pdf(paste0(Fig.Path, "DAPCPieBar_K", K,".", sp, "_NoVP.pdf"), width = 6.8, height = 4.74, useDingbats = FALSE)
print(DAPC.Combined)
dev.off()

jpeg(paste0(Fig.Path, "DAPCPieBar_K", K,".", sp, "_NoVP.jpg"), width = 6.8, height = 4.74, units = "in", res = 300)
print(DAPC.Combined)
dev.off()


###################################
# Merge dataframes across methods #
###################################

colnames(Q.matrix.Tess.df)[3:6] <- c("Clust1.Tess", "Clust2.Tess", "GenClust.Tess", "Prop.Tess")
colnames(Q.matrix.LEA.df)[3:6] <- c("Clust1.LEA", "Clust2.LEA", "GenClust.LEA", "Prop.LEA")
colnames(DAPC.df)[3:6] <- c("Clust1.DAPC", "Clust2.DAPC", "GenClust.DAPC", "Prop.DAPC")

Dh.Structure <- left_join(Q.matrix.Tess.df, Q.matrix.LEA.df[, -2], by = "sampleID")
Dh.Structure <- left_join(Dh.Structure, DAPC.df[, -2], by = "sampleID")

# Save Data frame with info across all analyses
write.csv(Dh.Structure, "Dh.Structure.Tess.LEA.DAPC.csv", row.names = FALSE)
