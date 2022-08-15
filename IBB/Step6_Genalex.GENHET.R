#############################################
#                 STEP 6                    #  
#              Filter on HWE                #
#      Output Genalex formatted file        #
#             Plot Mantel results           #
#    Calculate individual heterozygosity    #
#         Data set: Dasyurus hallucatus     #
#            Author: Robyn Shaw             #
#             Date: 22/02/2022              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

#############################
# LOAD IN REQUIRED PACKAGES #
#############################

library(dartR)
library(splitstackshape)
library(readxl)
library(ggplot2)
library(stringr)
library(gtools)


############################
# Filter SNPs based on HWE #
############################

# Load in genlight (post-filtering, related individuals removed) with clustering info/metadata
load("IBB/IBB_outputs/Cleaned.Unrelated.RmDup.gl.Dh_ClustInfo.rdata")

# Remove admixed individuals (based on Tess)
Dh.gl.rmAd <- Dh.gl.rmDupCoords[!(Dh.gl.rmDupCoords@other$ind.metrics$TessAdmixed == "Y"), ]
# Change pop to Genetic Clusters
Dh.gl.rmAd@pop <- as.factor(Dh.gl.rmAd@other$ind.metrics$TessGenCluster)

# Calculate HWE (with Bonferroni Corrected significance)
HWE_Clusts <- gl.report.hwe(x = Dh.gl.rmAd,
                            plot.out = FALSE,
                            multi_comp = TRUE, 
                            multi_comp_method = "bonferroni")
HWE_Clusts$Locus <- as.character(HWE_Clusts$Locus)

# Replace "sig" (significant) with a one, and everything else with a zero Note that I am using the bonferroni significance, rather than just the 0.05 significance which is why there are ns results.
HWE_Clusts$BonSig_1.0 <- ifelse(HWE_Clusts$Sig.adj == "sig", 1, 0)

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
sum(Locus_HWE$PropPopsSig.HWE >= 0.5) # Loci just out in one pop 
sum(Locus_HWE$PropPopsSig.HWE > 0.5) # Loci out in two pops

# In this case I'm going to remove everything because there are only two (real) pops
gl.HWE <- Dh.gl.rmAd[, !(Dh.gl.rmAd@loc.names %in% Locus_HWE$Locus[Locus_HWE$PropPopsSig.HWE >= 0.5])]


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
Genalex.Final <- cbind(Genalexout_IDs, NA_Col, gl.HWE@other$ind.metrics[, c("UTMX", "UTMY")], NA_Col, gl.HWE@other$ind.metrics[, c("id", "pop", "lat", "lon", "Sex", "Year", "IBRA_SubRegion", "Original_id", "Buffer_15km")])

# Write genalex formatted file
write.csv(x = Genalex.Final, file = "IBB/IBB_outputs/Dh.HWE.Filt.genalex.csv", row.names = FALSE, na = "")
# NOTE: Population summary statistics and Mantel statistics calculated in Genalex using this file (saved as excel workbook)


###########################################################
# Make a nice plot for Mantel test carried out in Genalex #
###########################################################

# Read in Genalex output
mantel <- read_xlsx("IBB/IBB_outputs/Dh.HWE.Filt.genalex.xlsx", sheet = "Mantel_Data")

# Convert metres to kms
mantel$Geographic_Distance <- mantel$Geographic_Distance/1000

# Format labels for plot
poplab <- c("Cluster 1", "Cluster 2")
names(poplab) <- c("Cluster1", "Cluster2")

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


pdf("IBB/IBB_outputs/Mantel.Pop.Plots.pdf", width = 8, height = 5, useDingbats=FALSE)
mantel.plot
dev.off()



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

# This function calculates several metrics.
# I will use the proportion of heterozygous loci (PHt), standardized heterozygosity relative to the mean expected heterozygosity (Hs_exp, Coltman 1999)

# Load in total dataset (i.e. before duplicated coordinates were removed)
# Then remove loci out of HWE
load("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.gl.Dh.rdata")
Dh.gl@other$ind.metrics <- read.csv("SNP_Filtering/SampleClean_outputs/Cleaned.Unrelated.Ind.metadata.Dh.csv")
sum(!Dh.gl@other$ind.metrics$id == Dh.gl@ind.names) # 0 if all match
Dh.gl <- Dh.gl[, Dh.gl@loc.names %in% gl.HWE@loc.names]

# Reformat data frame for GENHET
GenDat.genhet <- as.matrix(Dh.gl[,])

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

# Convert to data frame
IndHet.df <- data.frame(IndHet, stringsAsFactors = FALSE)
# Bring over Buffer info
IndHet.df <- cbind(IndHet.df, Dh.gl@other$ind.metrics[match(IndHet.df$sampleid, Dh.gl@other$ind.metrics$id), c("pop", "Buffer_15km")])
# Pre-pend with "B15" to denote 15km buffer 
IndHet.df$Buffer_15km <- ifelse(IndHet.df$Buffer_15km < 10, paste0("B15.0", IndHet.df$Buffer_15km), paste0("B15.", IndHet.df$Buffer_15km))
# Rename buffer for dolphin island samples, as 15kms groups them with mainland individuals, and I want to visualise them separately
IndHet.df$Buffer_15km[IndHet.df$pop == "Dolphin Island"] <- "B15.00"
# Classify as numeric for calculating the mean (currently "character")
IndHet.df$PHt <- as.numeric(IndHet.df$PHt)

# Group by pair and then take the mean
IndHet.df.mean <- IndHet.df[, c("Buffer_15km", "PHt")] %>%
  group_by(Buffer_15km) %>%
  summarise_all(mean)

# Get SD for means
IndHet.df.mean$PHt.HotCold <- (IndHet.df.mean$PHt - mean(IndHet.df.mean$PHt))/sd(IndHet.df.mean$PHt)
# Any hotspots?
sum(IndHet.df.mean$PHt.HotCold > 1.5) # 0 - no
# Any coldspots?
sum(IndHet.df.mean$PHt.HotCold < -1.5) # 2 locations
# Which Buffer locations
IndHet.df.mean$Buffer_15km[IndHet.df.mean$PHt.HotCold < -1.5]

# How many individuals sampled in the coldspot locations
sum(Dh.gl@other$ind.metrics$Buffer_15km %in% c("1")) # Mean over 2 individuals
sum(Dh.gl@other$ind.metrics$pop %in% c("Dolphin Island")) # Mean over 14 individuals

# Get SD for Inds
IndHet.df$PHt.HotCold <- (IndHet.df$PHt - mean(IndHet.df$PHt))/sd(IndHet.df$PHt)
# Any individuals with higher than average genetic diversity?
sum(IndHet.df$PHt.HotCold > 1.5) # 1 individual
# Any individuals with lower than average genetic diversity?
sum(IndHet.df$PHt.HotCold < -1.5) # 19 individuals

# Output results
write.csv(IndHet.df, "IBB/IBB_outputs/IndHet.df15km.csv", row.names = FALSE)
write.csv(IndHet.df.mean, "IBB/IBB_outputs/IndHet.df.15km_mean.csv", row.names = FALSE)
