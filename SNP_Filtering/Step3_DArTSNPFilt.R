#############################################
#                 STEP 3                    #  
# DArTSeq SNPs: Visualisation and Filtering #
#         Data set: Dasyurus hallucatus     #
#            Author: Robyn Shaw             #
#             Date: 03/06/2022              #
#############################################

setwd("/Users/robynshaw/Google Drive/Work/Rcode_repos/NorthernQuolls_LifeHistory2Landscape/")

############
# Packages #
############

library(dartR)
library(ggplot2)
library(flextable)
library(officer)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(viridis)
library(gridExtra)
library(ggpubr)
library(grid)
library(corrplot)
library(SNPRelate)
library(rgdal)
library(ggsn)
library(sessioninfo)
library(details)


#####################################################
# Input variables and check population sample sizes #
#####################################################

# Define path to data folder, file names and variables 

# Path to folder containing data
datpath <- "SNP_Filtering/Data/"

# Define prefix for naming output files
Outfile_Desc <- "Dhal"

# Create a new directory/sub folder for saving outputs (if it doesn't already exist)
outpath <- "SNP_Filtering/Filtering_outputs/"
ifelse(!dir.exists(file.path(outpath)), dir.create(file.path(file.path(outpath))), FALSE)


# Name of DArT SNP data (1-row)
DartGenos_1Row <- "Report_DDasy19-4717_1Row_NameEdit.csv"

# Name of last metadata column in SNP files before samples start
LastGeno.MD.Col <- "RepAvg"

# Name of DArT read count data and 2-row genotype file
# Make sure that sample names are in the same row as the rest of the headers before loading in
DartGenos_2Row <- "Report_DDasy19-4717_2Row_NameEdit.csv"
DartReadCounts <- "ReadCounts_DDasy19-4717_NameEdit.csv" 

# Name of last metadata column in read count file before samples start
LastRC.MD.Col <- "TotalPicRepSnpTest"

# Name of individual meta data file
# Must include column for: sample IDs (that match IDs in DArT SNP file) and populations
Ind_Metadat <- "Dasyurus_hallucatus_ind.metadata.csv"
pop_col <- "pop"

# Have DArT SNPs been blasted to a reference? (Y/N)
BLAST <- "Y"

# Add contig info from NCBI if available
# This will allow us to see which chromosome the SNPs are on (rather than just the contig name)
# Also add DArT column names for the chromosome/contig ID, position, alignment count and e value
if(BLAST == "Y") {
  # If sex chromosomes are known, list here (so they can be filtered out)
  SexChrom <- "x"
  ChromInfo.file <- "TasDevil7.1_Scaffold_Info.csv"
  ContigCol_dart <- "Chrom_Tasmanian_devil_v7"
  ContigPos_dart <- "ChromPos_Tasmanian_devil_v7"
  ContigCol_ncbi <- "Devil7.0.SCAFFOLD"
  ChromPos_ncbi <- "CHROMOSOME.POSITION"
  EVal.col <- "AlnEvalue_Tasmanian_devil_v7"
  Cnt.col <- "AlnCnt_Tasmanian_devil_v7"
}


# Path to shapefile for sample map
Map.shp <- "PilbaraIBRA.shp"

# Define mapping variables

# Map extent 
# This will be used to place the north arrow and the scale bar
# So make it slightly smaller than the actual map extent
xmin <- 115
xmax <- 122
ymin <- -24
ymax <- -19.5

# Define fill variable for map (e.g. here it is IBRA sub-regions)
MapFill <- "SUB_NAME_7"

# Define variables for map scale bar (distance unit, crs, distance increments)
distU <- "km" 
CoordRef <- 'WGS84'
dist <- 100


# Determine population sample size cutoff
# Read in individual metadata
ind.metrics <- read.csv(paste0(datpath, Ind_Metadat), row.names = 1)
PopNo <- data.frame(table(ind.metrics[, pop_col]))

# Create plot
popNo.plot <- ggplot(data = PopNo, aes(x = Var1, y = Freq, fill = Freq)) +
  geom_bar(stat = "identity", colour = "black", size = 0.2) +
  scale_fill_viridis(discrete = F, option = "D", direction = -1) +
  labs(x = "Population", y = "No. Individuals", fill = "Pop Size") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 6))

# View
popNo.plot


# Create Table
colnames(PopNo) <- c("Population", "N")
PopTab <- flextable(data = PopNo)
PopTab <- autofit(PopTab)
PopTab <- bg(PopTab, bg = "lightgrey", part = "header") 
PopTab <- bold(PopTab, part = "header")
PopTab <- align(PopTab, align = "left")
PopTab <- border_outer(PopTab, part = "header", border = fp_border(color = "black",style = "solid", width = 0.5))
PopTab <- border_outer(PopTab, part = "body", border = fp_border(color = "black",style = "solid", width = 0.5))
PopTab <- fontsize(PopTab, size = 10)

# View
PopTab


# Use figure/table above to decide which populations to include in calculating population level metrics. If populations have a small sample size, the population genetic metrics are less robust and will throw off mean values. Choose a threshold that makes sense for the data set, at least ~10 if possible (ideally at least three populations with sample sizes at or above 30). It's worth noting that the population summary statistics calculated in this script are only used to visualise the data (in other words, filtering is not based on these population groupings and their associated summary statistics). However, it's still a good idea to decide on reasonable populations, to see if any strange population genetic patterns have real biological meaning, or if they are instead correlated with other measures such as quality.  

# In this case, there are quite a few populations with very good sample sizes. Therefore, the sampling site will be used as the 'population'.    

# If the sample sizes in defined pops are too low, 
# choose a different way to group samples 
# (i.e. merge some pops, use regional info etc.)

# Do you want to choose a different population grouping? (Y/N)
Change_pop <- "N"

if(Change_pop == "Y"){
  # Choose a new column to group samples by (in individual metadata file)
  pop_col <- "IBRA_SubRegion"
  PopNo <- data.frame(table(ind.metrics[, pop_col]))

  popNo.plot <- ggplot(data = PopNo, aes(x = Var1, y = Freq, fill = Freq)) +
  geom_bar(stat = "identity", colour = "black", size = 0.2) +
  scale_fill_viridis(discrete = F, option = "D", direction = -1) +
  labs(x = "Population", y = "No. Individuals", fill = "Pop Size") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 6))
  
  popNo.plot
  
  
  colnames(PopNo) <- c("Population", "N")
  PopTab <- flextable(data = PopNo)
  PopTab <- autofit(PopTab)
  PopTab <- bg(PopTab, bg = "lightgrey", part = "header") 
  PopTab <- bold(PopTab, part = "header")
  PopTab <- align(PopTab, align = "left")
  PopTab <- border_outer(PopTab, part = "header", 
                         border = fp_border(color = "black",style = "solid", width = 0.5))
  PopTab <- border_outer(PopTab, part = "body", 
                         border = fp_border(color = "black",style = "solid", width = 0.5))
  PopTab <- fontsize(PopTab, size = 10)

}


# Define population number threshold
PopCutOff <- 10

# List pop names above this threshold
PopKeep <- names(table(ind.metrics[, pop_col]) >= PopCutOff)[table(ind.metrics[, pop_col]) >= PopCutOff]

# This leaves:
cat(length(PopKeep), "populations")
# to use in the calculation of population genetic summary statistics.


###########################################################
# Data preparation and filtering on individual read count #
###########################################################

# Is it necessary to filter on individual read count, or is it okay to filter on the mean locus values routinely provided by DArT in their meta data instead?

# Prepare data
# Note that data were obtained from DArT in 2019-2020, and formats may have changed
# DArT provides read counts for each individual at each locus. This data is in two-row format, with read counts provided for the reference and SNP allele. The data is similar to the two-row SNP format file, however the cloneIDs do not match between files. Therefore, a bit of data wrangling is needed to compare individual read counts to called genotypes.

# Add a unique column (sequence and SNP position) to match loci across data sets
# Read in data as dataframe (instead of genlight)
# DArT read count data csv (skip top 'no data' rows)
ReadCounts <- read.csv(file = paste0(datpath, DartReadCounts), skip = sum(read.csv(paste0(datpath, DartReadCounts), header = FALSE)[,1] == "*", na.rm = TRUE), header = TRUE, stringsAsFactors = FALSE)

# DArT 1-row Genotype csv (skips top 'no data' rows)- easier to manipulate as data frame versus genlight in this case
Genos_1row <- read.csv(file = paste0(datpath, DartGenos_1Row), na.strings = "-", skip = sum(read.csv(paste0(datpath, DartGenos_1Row), header = FALSE)[,1] == "*", na.rm = TRUE), header = TRUE, stringsAsFactors = FALSE)

# DArT 2-row Genotype csv (skips top 'no data' rows)- easier to manipulate as data frame versus genlight in this case
Genos_2row <- read.csv(file = paste0(datpath, DartGenos_2Row), na.strings = "-", skip = sum(read.csv(paste0(datpath, DartGenos_2Row), header = FALSE)[,1] == "*", na.rm = TRUE), header = TRUE, stringsAsFactors = FALSE)

# The clone/allele ids do not match between the geno and read count data sets.
# Add a unique sequence column that combines the sequence with the SNP position.
# This will ensure that the correct ref sequence is paired with the 
# corresponding SNP sequence when there are multiple SNPs per sequence.
ReadCounts <- add_column(ReadCounts, UniqueSeq = paste(ReadCounts$AlleleSequence, ReadCounts$SnpPosition, sep = "_"), .before = 1)
Genos_1row <- add_column(Genos_1row, UniqueSeq = paste(Genos_1row$AlleleSequence, Genos_1row$SnpPosition, sep = "_"), .before = 1)
Genos_2row <- add_column(Genos_2row, UniqueSeq = paste(Genos_2row$AlleleSequence, Genos_2row$SnpPosition, sep = "_"), .before = 1)

# Check that the number of loci matches across data sets
# Check if these unique sequences match between csv files
length(setdiff(ReadCounts$UniqueSeq, Genos_2row$UniqueSeq))
length(setdiff(Genos_2row$UniqueSeq, ReadCounts$UniqueSeq))


# There are some unique allele sequences that are in one data set, but not the other. I contacted DArT to ask about this mismatch and got this reply from Jason Carling (Lead, Analytical Services):  

# "There will sometimes be relatively small differences in the marker sets present in the files which have 1 data column per assay vs those which have 1 data column per sample, due to the application of automated marker quality filtering. In these two file types the underlying genotype frequencies measured across the markers can be very slightly different due to the difference in sample number (reps present vs consensus only). In these cases the markers which are not in common across the files will be those which are close to the filtering threshold value in one or more metadata parameters used in filtering."

# Therefore, just filter these out.  

# Filter out mismatched loci, create a new cloneID that matches across all data sets (so don't have to use the entire sequence) and create a column specifying if reference/SNP allele (for matching with one-row format)

# Read count data frame: Only keep Unique Sequences that are in both files
ReadCounts.SeqFilt <- ReadCounts[ReadCounts$UniqueSeq %in% intersect(ReadCounts$UniqueSeq, Genos_2row$UniqueSeq), ]

# 1-row Genotype data frame: Only keep Unique Sequences that are in both files
Genos_1row.SeqFilt <- Genos_1row[Genos_1row$UniqueSeq %in% intersect(Genos_1row$UniqueSeq, ReadCounts$UniqueSeq), ]

# 2-row Genotype data frame: Only keep Unique Sequences that are in both files
Genos_2row.SeqFilt <- Genos_2row[Genos_2row$UniqueSeq %in% intersect(Genos_2row$UniqueSeq, ReadCounts$UniqueSeq), ]

# Create unique cloneID per SNP by appending SNP position
# Identify Snps versus reference alleles
# Copy this identifier over to read count DF, matching by sequence (for ease of comparing DFs)

Genos_2row.SeqFilt <- add_column(Genos_2row.SeqFilt, Locus = paste(Genos_2row.SeqFilt$CloneID, Genos_2row.SeqFilt$SnpPosition, sep = "_"), .before = 1)
Genos_2row.SeqFilt <- add_column(Genos_2row.SeqFilt, SNP.Ref = Genos_2row.SeqFilt$SNP, .after = "Locus")
Genos_2row.SeqFilt$SNP.Ref[Genos_2row.SeqFilt$SNP.Ref == ""] <- "Ref" 

# The SNP column provided by Dart is empty if it's a ref allele, but contains a position/base if it's a SNP allele. 
# Therefore, rename everything that's empty as "Ref", everything else as "SNP"

Genos_2row.SeqFilt$SNP.Ref[!Genos_2row.SeqFilt$SNP.Ref == "Ref"] <- "SNP" 

table(Genos_2row.SeqFilt$SNP.Ref) # Make sure even numbers of each

# Add unique CloneIDs, and Ref/SNP info (from the 2-row geno DF) to 1-row geno DF
Geno2Geno_Match <- match(Genos_1row.SeqFilt$UniqueSeq, Genos_2row.SeqFilt$UniqueSeq)
Genos_1row.SeqFilt <- add_column(Genos_1row.SeqFilt, Locus = Genos_2row.SeqFilt$Locus[Geno2Geno_Match], .before = 1)
sum(is.na(Genos_1row.SeqFilt$Locus))

# Add unique CloneIDs, and Ref/SNP info (from the 2-row geno DF) to the read count DF
Geno2Count_Match <- match(ReadCounts.SeqFilt$UniqueSeq, Genos_2row.SeqFilt$UniqueSeq)
ReadCounts.SeqFilt <- add_column(ReadCounts.SeqFilt, Locus = Genos_2row.SeqFilt$Locus[Geno2Count_Match], .before = 1)
ReadCounts.SeqFilt <- add_column(ReadCounts.SeqFilt, SNP.Ref = Genos_2row.SeqFilt$SNP.Ref[Geno2Count_Match], .after = "Locus")

table(ReadCounts.SeqFilt$SNP.Ref) # Make sure even numbers of each (ref vs snp)


# Remove replicate samples by aggregating duplicates and keeping the max read count  

# DArT runs 30 technical replicates, which means that there are 30 duplicate samples (per plate) in the read count data set. However, there is only one genotype called for each duplicated sample in the SNP data sets. Unsure how DArT comes up with a consensus across repeated samples, so take the maximum value (i.e. for each row, take the maximum read count between the original and the repeat).  


# READ COUNT DF

# Remove extra metadata info from read count data frame, order rows by CloneID and columns by sample name
LastRC.MD.No <- which(colnames(ReadCounts.SeqFilt) == LastRC.MD.Col)
RC_Rm_MD <- ReadCounts.SeqFilt[ , c("Locus", "SNP.Ref", colnames(ReadCounts.SeqFilt)[(LastRC.MD.No + 1):ncol(ReadCounts.SeqFilt)])]
RC_order <- RC_Rm_MD[order(RC_Rm_MD$Locus), c(1:2, order(colnames(RC_Rm_MD)[3:ncol(RC_Rm_MD)])+2)]

# Find column number where sample names start in the 2-row data set (for comparing sample names)
LastGeno.MD.No.2 <- which(colnames(Genos_2row.SeqFilt) == LastGeno.MD.Col)

# Find technical replicates (i.e. duplicated samples)
# Duplicate sample IDs are automatically appended with ". + dup number" in R (so that column names aren't repeated). 
# Remove these characters to find duplicated sample IDs
col.rename <- colnames(RC_order)

# Find samples names that don't match between data sets (for example, because they've been appended with ".1")
reps <- which(col.rename %in% setdiff(col.rename[3:length(col.rename)], colnames(Genos_2row.SeqFilt[, LastGeno.MD.No.2:ncol(Genos_2row.SeqFilt)])))

# Remove the last two characters and replace these sample names
col.rename[reps] <- substr(col.rename[reps], 1, nchar(col.rename[reps])-2)

# Create loop to find the column numbers for identical samples (technical replicates)
# Take the maximum read count value (across columns) for each row and add this aggregate column to the end of the data frame

ColNo <- 0
col_ID <- colnames(RC_order)
DupNo <- as.data.frame(table(col.rename), stringsAsFactors = FALSE)
rmCol_Vect <- vector()

for(ColNo in (1:nrow(DupNo)-1)) {
  Col <- ColNo + 1
  if (DupNo[Col, 2] == 1) {
  } else {
    rmCol_Vect <- c(rmCol_Vect, which(col.rename %in% DupNo[Col, 1]))
    RC_order <- cbind(RC_order, apply(RC_order[, which(col.rename %in% DupNo[Col, 1])], 1, max, na.rm = TRUE))
    colnames(RC_order)[ncol(RC_order)] <- DupNo[Col, 1]
  }
}

# Keep only unique column numbers
rmCol_Vect <- unique(rmCol_Vect)

# Remove replicate and original sample from data frame (just keep new max read count column), then order sample columns
RC_MaxRep <- RC_order[, -(rmCol_Vect)]
RC_MaxRep <- RC_MaxRep[, c(1:2, order(colnames(RC_MaxRep)[3:ncol(RC_MaxRep)])+2)]
colnames(RC_MaxRep) # Check

# 2-ROW GENO DF

# Order rows by CloneID and columns by sample name
G2_Final <- Genos_2row.SeqFilt[order(Genos_2row.SeqFilt$Locus), c(1:2, order(colnames(Genos_2row.SeqFilt)[(LastGeno.MD.No.2 + 1):ncol(Genos_2row.SeqFilt)]) + LastGeno.MD.No.2)]

# Check that order of both data frames match exactly 
sum(colnames(G2_Final) == colnames(RC_MaxRep)) == dim(G2_Final)[2]
sum(G2_Final$Locus == RC_MaxRep$Locus) == dim(G2_Final)[1]

# 1-ROW GENO DF
LastGeno1.MD.No <- which(colnames(Genos_1row.SeqFilt) == LastGeno.MD.Col)
G1_Final <- Genos_1row.SeqFilt[order(Genos_1row.SeqFilt$Locus), c(1, order(colnames(Genos_1row.SeqFilt)[(LastGeno1.MD.No + 1):ncol(Genos_1row.SeqFilt)])+LastGeno1.MD.No)]

# Check that SNP numbers add up/sample IDs match
2*dim(G1_Final)[1] == dim(G2_Final)[1] # Should be 2x as many SNPs (ref+snp) in 2-row Geno DF
sum(!colnames(G1_Final) == colnames(G2_Final)[-2])

# Screen out read counts for NA genotypes

# Create a new data frame where read count data for genotypes that haven't been called are screened out (i.e. replace count value with NA if same sample/loci is missing in SNP data set). Don't include these in the read count summary statistics, instead account for these through missing data filters.

RC_Final.Raw <- map2_df(G2_Final, RC_MaxRep, ~ifelse(is.na(.x), .x, .y))
RC_Final.Raw <- as.data.frame(RC_Final.Raw)

# Create sample/loci/metric data frames

# Reshape data, calculate summary statistics and visualise individual read count.

# Make new data frame with total read counts
Tot.Counts <- RC_Final.Raw[,-2]
Tot.Counts <- group_by(Tot.Counts, Locus)
Tot.Counts <- summarise_each(Tot.Counts, sum)
Tot.Counts <- Tot.Counts[order(Tot.Counts$Locus), ]

# Make new data frame with ref read counts
Ref.Counts <- RC_Final.Raw[RC_Final.Raw$SNP.Ref == "Ref", -2]
Ref.Counts <- Ref.Counts[order(Ref.Counts$Locus), ]

# Make new data frame with SNP read counts
SNP.Counts <- RC_Final.Raw[RC_Final.Raw$SNP.Ref == "SNP", -2]
SNP.Counts <- SNP.Counts[order(SNP.Counts$Locus), ]

# Make a new data frame with locus metrics to append to genlight
loc.metrics <- Genos_1row.SeqFilt[order(Genos_1row.SeqFilt$Locus), 1:which(colnames(Genos_1row.SeqFilt)==LastGeno.MD.Col)]
sum(!loc.metrics$Locus == Tot.Counts$Locus) + sum(!loc.metrics$Locus == Ref.Counts$Locus) + sum(!loc.metrics$Locus == SNP.Counts$Locus)

# Calculate summary metrics across individuals for each locus (total read counts)
loc.metrics$RC_MeanTot <- apply(Tot.Counts[, 2:ncol(Tot.Counts)], 1, mean, na.rm = TRUE)
loc.metrics$RC_MeanRef <- apply(Ref.Counts[, 2:ncol(Ref.Counts)], 1, mean, na.rm = TRUE)
loc.metrics$RC_MeanSNP <- apply(SNP.Counts[, 2:ncol(SNP.Counts)], 1, mean, na.rm = TRUE)

# INDIVIDUAL METRICS

# Reshape each data frame so that 1 row per locus per sample (i.e. from 'wide' to 'long' format)
Locus.Samp.TotRC <- gather(data = Tot.Counts, Sample_ID, Total_ReadCount, -Locus)
Locus.Samp.RefRC <- gather(data = Ref.Counts, Sample_ID, Ref_ReadCount, -Locus)
Locus.Samp.SNPRC <- gather(data = SNP.Counts, Sample_ID, SNP_ReadCount, -Locus)
Locus.Samp.Geno <- gather(data = G1_Final, Sample_ID, Genotype, -Locus)

# Create individual level read count data frame
IndRC_DF <- left_join(Locus.Samp.Geno, Locus.Samp.TotRC, by = c("Locus", "Sample_ID"))
IndRC_DF <- left_join(IndRC_DF, Locus.Samp.RefRC, by = c("Locus", "Sample_ID"))
IndRC_DF <- left_join(IndRC_DF, Locus.Samp.SNPRC, by = c("Locus", "Sample_ID"))

# Reshape data for comparing read count by proportion of each genotype
Reads.by.geno <- table(IndRC_DF$Total_ReadCount, IndRC_DF$Genotype, useNA = "no")
colnames(Reads.by.geno) <- c("Ref.Hom", "SNP.Hom", "Het")
Reads.by.geno <- cbind(Reads.by.geno, Total.Called.Genos = rowSums(Reads.by.geno))
Reads.by.geno <- as.data.frame(Reads.by.geno)
Reads.by.geno <- add_column(Reads.by.geno, Total_ReadCount = as.integer(row.names(Reads.by.geno)), .before = 1)
Reads.by.geno <- cbind(Reads.by.geno, Prop.Ref.Hom = Reads.by.geno$Ref.Hom/Reads.by.geno$Total.Called.Genos, Prop.SNP.Hom = Reads.by.geno$SNP.Hom/Reads.by.geno$Total.Called.Genos, Prop.Het = Reads.by.geno$Het/Reads.by.geno$Total.Called.Genos)

# Calculate allele balance (compare read count by looking at proportion of ref reads compared to total across hets)

# INDIVIDUAL METRICS

HetRC <- filter(IndRC_DF, Genotype == 2)
HetRC.T <- HetRC %>% group_by(Locus) %>% summarise(AverageT = mean(Total_ReadCount))
HetRC.R <- HetRC %>% group_by(Locus) %>% summarise(AverageR = mean(Ref_ReadCount))
HetRC.Final <- left_join(HetRC.T, HetRC.R, by = "Locus")
HetRC.Final$AB <- HetRC.Final$AverageR/HetRC.Final$AverageT
loc.metrics <- left_join(loc.metrics, HetRC.Final[, c(1,4)], by = "Locus")

# Reformat data for facet_wrap in ggplot
Reads.by.geno.reshape <- pivot_longer(data = Reads.by.geno[, c("Total_ReadCount", "Prop.Ref.Hom", "Prop.SNP.Hom", "Prop.Het")], cols = c("Prop.Ref.Hom", "Prop.SNP.Hom", "Prop.Het"), names_to = "Prop_Geno")

# Rename variables in plotting order
Reads.by.geno.reshape$Prop_Geno[Reads.by.geno.reshape$Prop_Geno == "Prop.Ref.Hom"] <- "a_Prop.Ref.Hom"
Reads.by.geno.reshape$Prop_Geno[Reads.by.geno.reshape$Prop_Geno == "Prop.SNP.Hom"] <- "b_Prop.Ref.Hom"
Reads.by.geno.reshape$Prop_Geno[Reads.by.geno.reshape$Prop_Geno == "Prop.Het"] <- "c_Prop.Ref.Het"

PropGeno.labs <- c("Ref Homozygotes", "SNP Homozygotes", "Heterozygotes")
names(PropGeno.labs) <- unique(Reads.by.geno.reshape$Prop_Geno)

RC.geno.Raw <- ggplot(Reads.by.geno.reshape, aes(x = Total_ReadCount, y = value, colour = as.factor(Prop_Geno))) +
  geom_point(size = 0.5, na.rm = T) +
  facet_wrap(~ Prop_Geno, nrow = 3, dir = "h", strip.position = "right", labeller = labeller(Prop_Geno = PropGeno.labs)) +
  labs(x = "\nTotal read count", y = "Proportion of called genotypes (raw)\n") +
  scale_x_continuous(breaks = seq(0, max(Reads.by.geno.reshape$Total_ReadCount), 50)) +
  scale_colour_brewer(palette="Dark2") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        legend.position = "none")

loc.metrics.raw <- loc.metrics

# View
RC.geno.Raw + labs(title = "Raw individual read counts at each locus") + theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5))


# Decide on individual read count threshold
# Inspect above figure to decide at which min/max values the genotypes are no longer strongly impacted by read count, then filter individual genotypes on these values. It might be expected that the proportion of each genotype should be stable, regardless of read count. However, at the extremes, low read counts increase the likelihood that one of the alleles hasn't been sequenced (likely the minor allele, which is usually the SNP), and high read counts might indicate that there are multiple copies of this sequence spread throughout the genome (i.e. paralogues).  

# In this data set, there is a strong relationship between individual read count and the proportion of each called genotype at the lower and upper extremes. At the lower end of the range, heterozygotes are under called, while SNP homozygotes are over called. After about 350 reads, the genotypes become fixed for the reference or the SNP allele (indicating that these may represent paralogous regions).  

# DEFINE MIN MAX VALUES
MinRC <- 4
MaxRC <- 300

# Create a new data frame where count data for genotypes that haven't been called are screened out (i.e. replace count value with NA if same sample/loci is missing in geno file)

Tot.Counts <- as.data.frame(RC_Final.Raw[,-2] %>% group_by(Locus) %>% summarise_each(sum))
rownames(Tot.Counts) <- Tot.Counts[, 1]
Tot.Counts <- Tot.Counts[, -1]

G1_Final_screen <- as.matrix(G1_Final[ ,-1])
rownames(G1_Final_screen) <- G1_Final$Locus

# Check that data frames match exactly
sum(!rownames(Tot.Counts) == rownames(G1_Final_screen))
sum(!colnames(Tot.Counts) == colnames(G1_Final_screen))

# Screen out genotypes above/below threshold
Geno_RC.Filt <- ifelse(Tot.Counts < MinRC | Tot.Counts > MaxRC, NA, G1_Final_screen)

# Screen out read counts using genotype file above (screen out counts for genotypes = NA)
# Cant just do this by removing counts below/above thresholds because that will replace read counts of 0 with NA
# Read counts of 0 are important, because there might be 0 reads of the SNP, but 10 reads of the reference
# Therefore, have to use a slightly convoluted workaround

# Save screened/filtered genotype file as data frame
Geno_RC.Filt <- as.data.frame(Geno_RC.Filt)

# Separate read count data frame into ref and SNPs
RC_MaxRep.SNP <- RC_MaxRep[RC_MaxRep$SNP.Ref == "SNP", -c(1,2)]
RC_MaxRep.Ref <- RC_MaxRep[RC_MaxRep$SNP.Ref == "Ref", -c(1,2)]

row.names(RC_MaxRep.SNP) <- RC_MaxRep[RC_MaxRep$SNP.Ref == "SNP", 1]
row.names(RC_MaxRep.Ref) <- RC_MaxRep[RC_MaxRep$SNP.Ref == "Ref", 1]

# Make sure loci match/data frames are identical
sum(!row.names(RC_MaxRep.SNP) == row.names(RC_MaxRep.Ref))
sum(!row.names(RC_MaxRep.SNP) == row.names(Geno_RC.Filt))
sum(!row.names(RC_MaxRep.Ref) == row.names(Geno_RC.Filt))

sum(!colnames(RC_MaxRep.SNP) == colnames(RC_MaxRep.Ref))
sum(!colnames(RC_MaxRep.SNP) == colnames(Geno_RC.Filt))
sum(!colnames(RC_MaxRep.Ref) == colnames(Geno_RC.Filt))

# Create a new data frame where count data for genotypes that were screened out is replaced with NA
RC_Final_RefFilt <- map2_df(Geno_RC.Filt, RC_MaxRep.Ref, ~ifelse(is.na(.x), .x, .y))
RC_Final_SNPFilt <- map2_df(Geno_RC.Filt, RC_MaxRep.SNP, ~ifelse(is.na(.x), .x, .y))

RC_Final_RefFilt <- cbind(RC_MaxRep[RC_MaxRep$SNP.Ref == "Ref", 1:2], RC_Final_RefFilt)
RC_Final_SNPFilt <- cbind(RC_MaxRep[RC_MaxRep$SNP.Ref == "SNP", 1:2], RC_Final_SNPFilt)

# Merge data frames
RC_Final_RCFilt <- rbind(RC_Final_RefFilt, RC_Final_SNPFilt)
RC_Final_RCFilt <- RC_Final_RCFilt[order(RC_Final_RCFilt$Locus), ]

# Get Geno DF in correct format for function
Geno_RC.Filt.loc <- as.data.frame(Geno_RC.Filt)
Geno_RC.Filt.loc$Locus <- rownames(Geno_RC.Filt.loc)


# Reshape data, calculate summary statistics and visualise individual read count.

# Make new data frame with total read counts
Tot.Counts <- RC_Final_RCFilt[,-2]
Tot.Counts <- group_by(Tot.Counts, Locus)
Tot.Counts <- summarise_each(Tot.Counts, sum)
Tot.Counts <- Tot.Counts[order(Tot.Counts$Locus), ]

# Make new data frame with ref read counts
Ref.Counts <- RC_Final_RCFilt[RC_Final_RCFilt$SNP.Ref == "Ref", -2]
Ref.Counts <- Ref.Counts[order(Ref.Counts$Locus), ]

# Make new data frame with SNP read counts
SNP.Counts <- RC_Final_RCFilt[RC_Final_RCFilt$SNP.Ref == "SNP", -2]
SNP.Counts <- SNP.Counts[order(SNP.Counts$Locus), ]

# Make a new data frame with locus metrics to append to genlight
loc.metrics <- Genos_1row.SeqFilt[order(Genos_1row.SeqFilt$Locus), 1:which(colnames(Genos_1row.SeqFilt) == LastGeno.MD.Col)]

# Check that order of loci matches among data frames
sum(!loc.metrics$Locus == Tot.Counts$Locus) + sum(!loc.metrics$Locus == Ref.Counts$Locus) + sum(!loc.metrics$Locus == SNP.Counts$Locus)

# Calculate summary metrics across individuals for each locus (total read counts)
loc.metrics$RC_MeanTot <- apply(Tot.Counts[, 2:ncol(Tot.Counts)], 1, mean, na.rm = TRUE)
loc.metrics$RC_MeanRef <- apply(Ref.Counts[, 2:ncol(Ref.Counts)], 1, mean, na.rm = TRUE)
loc.metrics$RC_MeanSNP <- apply(SNP.Counts[, 2:ncol(SNP.Counts)], 1, mean, na.rm = TRUE)

# INDIVIDUAL METRICS

# Reshape each data frame so that 1 row per locus per sample (i.e. from 'wide' to 'long' format)
Locus.Samp.TotRC <- gather(data = Tot.Counts, Sample_ID, Total_ReadCount, -Locus)
Locus.Samp.RefRC <- gather(data = Ref.Counts, Sample_ID, Ref_ReadCount, -Locus)
Locus.Samp.SNPRC <- gather(data = SNP.Counts, Sample_ID, SNP_ReadCount, -Locus)
Locus.Samp.Geno <- gather(data = Geno_RC.Filt.loc, Sample_ID, Genotype, -Locus)

# Create individual level read count data frame
IndRC_DF <- left_join(Locus.Samp.Geno, Locus.Samp.TotRC, by = c("Locus", "Sample_ID"))
IndRC_DF <- left_join(IndRC_DF, Locus.Samp.RefRC, by = c("Locus", "Sample_ID"))
IndRC_DF <- left_join(IndRC_DF, Locus.Samp.SNPRC, by = c("Locus", "Sample_ID"))

# Reshape data for comparing read count by proportion of each genotype
Reads.by.geno <- table(IndRC_DF$Total_ReadCount, IndRC_DF$Genotype, useNA = "no")
colnames(Reads.by.geno) <- c("Ref.Hom", "SNP.Hom", "Het")
Reads.by.geno <- cbind(Reads.by.geno, Total.Called.Genos = rowSums(Reads.by.geno))
Reads.by.geno <- as.data.frame(Reads.by.geno)
Reads.by.geno <- add_column(Reads.by.geno, Total_ReadCount = as.integer(row.names(Reads.by.geno)), .before = 1)
Reads.by.geno <- cbind(Reads.by.geno, Prop.Ref.Hom = Reads.by.geno$Ref.Hom/Reads.by.geno$Total.Called.Genos, Prop.SNP.Hom = Reads.by.geno$SNP.Hom/Reads.by.geno$Total.Called.Genos, Prop.Het = Reads.by.geno$Het/Reads.by.geno$Total.Called.Genos)

# Calculate allele balance (compare read count by looking at proportion of ref reads compared to total across hets)

# INDIVIDUAL METRICS

HetRC <- filter(IndRC_DF, Genotype == 2)
HetRC.T <- HetRC %>% group_by(Locus) %>% summarise(AverageT = mean(Total_ReadCount))
HetRC.R <- HetRC %>% group_by(Locus) %>% summarise(AverageR = mean(Ref_ReadCount))
HetRC.Final <- left_join(HetRC.T, HetRC.R, by = "Locus")
HetRC.Final$AB <- HetRC.Final$AverageR/HetRC.Final$AverageT
loc.metrics <- left_join(loc.metrics, HetRC.Final[, c(1,4)], by = "Locus")

# Reformat data for facet_wrap in ggplot
Reads.by.geno.reshape <- pivot_longer(data = Reads.by.geno[, c("Total_ReadCount", "Prop.Ref.Hom", "Prop.SNP.Hom", "Prop.Het")], cols = c("Prop.Ref.Hom", "Prop.SNP.Hom", "Prop.Het"), names_to = "Prop_Geno")

# Rename variables in plotting order
Reads.by.geno.reshape$Prop_Geno[Reads.by.geno.reshape$Prop_Geno == "Prop.Ref.Hom"] <- "a_Prop.Ref.Hom"
Reads.by.geno.reshape$Prop_Geno[Reads.by.geno.reshape$Prop_Geno == "Prop.SNP.Hom"] <- "b_Prop.Ref.Hom"
Reads.by.geno.reshape$Prop_Geno[Reads.by.geno.reshape$Prop_Geno == "Prop.Het"] <- "c_Prop.Ref.Het"

PropGeno.labs <- c("Ref Homozygotes", "SNP Homozygotes", "Heterozygotes")
names(PropGeno.labs) <- unique(Reads.by.geno.reshape$Prop_Geno)

RC.geno.Filt <- ggplot(Reads.by.geno.reshape, aes(x = Total_ReadCount, y = value, colour = as.factor(Prop_Geno))) +
  geom_point(size = 0.5, na.rm = T) +
  facet_wrap(~ Prop_Geno, nrow = 3, dir = "h", strip.position = "right", labeller = labeller(Prop_Geno = PropGeno.labs)) +
  labs(x = "\nTotal read count", y = "Proportion of called genotypes (filtered)\n") +
  scale_x_continuous(breaks = seq(0, max(Reads.by.geno.reshape$Total_ReadCount), 50)) +
  scale_colour_brewer(palette="Dark2") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        legend.position = "none")

loc.metrics.filt <- loc.metrics

# View
RC.geno.Filt + labs(title = "Filtered individual read counts at each locus") + theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5))

# The above figure shows that individual calls that were most strongly impacted by read count have been removed. The stability of the called genotypes does start to dissolve at around 200-300 reads. However, there is a trade off between removing the effect of read count and introducing too much missing data; i.e. be pragmatic.  


# Create genlight objects for raw and filtered data sets

# Raw data set
# Recode Genotypes in standard format (i.e. not Dart format)
G1_Final_screen[G1_Final_screen== 1] <- 9
G1_Final_screen[G1_Final_screen== 2] <- 1
G1_Final_screen[G1_Final_screen== 9] <- 2

# Transpose so in correct format to load as genlight
G1_Final.t <- t(G1_Final_screen)

# Create new genlight with genotypes
gl.RC.raw <- new("genlight", G1_Final.t)

# Check that locus order matches between loc.metrics and new genlight so can append
sum(!gl.RC.raw@loc.names == loc.metrics.filt$Locus)
gl.RC.raw@other$loc.metrics <- loc.metrics.filt

# Match individual metrics to genlight
ind.metrics <- ind.metrics[match(gl.RC.raw@ind.names, rownames(ind.metrics)),]
sum(!rownames(ind.metrics) == gl.RC.raw@ind.names)

# Check that individual order matches between meta data and new genlight so can append
sum(!rownames(ind.metrics) == gl.RC.raw@ind.names)

# Append to ind.metrics
gl.RC.raw@other$ind.metrics <- ind.metrics
gl.RC.raw@other$latlong <- ind.metrics[ , 2:3]

# Fill in variables to match dartR's genlight
pop(gl.RC.raw) <- gl.RC.raw@other$ind.metrics[, pop_col]
gl.RC.raw@loc.all <- gsub(".*:","", gsub(">","/", gl.RC.raw@other$loc.metrics$SNP))
if(BLAST == "Y") {
  gl.RC.raw@chromosome <- as.factor(gl.RC.raw@other$loc.metrics[ , ContigCol_dart])
  gl.RC.raw@position <- gl.RC.raw@other$loc.metrics[ , ContigPos_dart]
}
gl.RC.raw@ploidy <- as.integer(rep(2, nInd(gl.RC.raw)))

# Recalculate other metrics (like call rate)
recalc.flags <- c("AvgPIC", "OneRatioRef", "OneRatioSnp", 
                  "PICRef", "PICSnp", "CallRate", "maf", "FreqHets", "FreqHomRef", 
                  "FreqHomSnp", "monomorphs", "OneRatio", "PIC")
gl.RC.raw@other$loc.metrics.flags <- data.frame(matrix(TRUE, 
                                                   nrow = 1, ncol = length(recalc.flags)))
names(gl.RC.raw@other$loc.metrics.flags) <- recalc.flags
gl.RC.raw <- gl.recalc.metrics(gl.RC.raw, mono.rm = FALSE)


# Filtered data set

# Recode Genotypes in standard format (i.e. not Dart format)
Geno_RC.Filt[Geno_RC.Filt== 1] <- 9
Geno_RC.Filt[Geno_RC.Filt== 2] <- 1
Geno_RC.Filt[Geno_RC.Filt== 9] <- 2

# Transpose so in correct format to load as genlight
Geno_RC.Filt.t <- t(Geno_RC.Filt)

# Create new genlight with filtered genotypes
gl.RC.filt <- new("genlight", Geno_RC.Filt.t)

# Check that locus order matches between loc.metrics and new genlight so can append
sum(!gl.RC.filt@loc.names == loc.metrics.filt$Locus)
gl.RC.filt@other$loc.metrics <- loc.metrics.filt

# Read in individual metadata
ind.metrics <- read.csv(paste0(datpath, Ind_Metadat), row.names = 1)
ind.metrics <- ind.metrics[match(gl.RC.filt@ind.names, rownames(ind.metrics)),]
sum(!rownames(ind.metrics) == gl.RC.filt@ind.names)

# Check that individual order matches between meta data and new genlight so can append
sum(!rownames(ind.metrics) == gl.RC.filt@ind.names)

# Append to ind.metrics
gl.RC.filt@other$ind.metrics <- ind.metrics
gl.RC.filt@other$latlong <- ind.metrics[ , 2:3]

# Fill in variables to match dartR's genlight
pop(gl.RC.filt) <- gl.RC.filt@other$ind.metrics[, pop_col]
gl.RC.filt@loc.all <- gsub(".*:","", gsub(">","/", gl.RC.filt@other$loc.metrics$SNP))
if(BLAST == "Y") {
  gl.RC.raw@chromosome <- as.factor(gl.RC.raw@other$loc.metrics[ , ContigCol_dart])
  gl.RC.raw@position <- gl.RC.raw@other$loc.metrics[ , ContigPos_dart]
}
gl.RC.filt@ploidy <- as.integer(rep(2, nInd(gl.RC.filt)))

# Recalculate other metrics (like call rate)
gl.RC.filt@other$loc.metrics.flags <- data.frame(matrix(TRUE, nrow = 1, ncol = length(recalc.flags)))
names(gl.RC.filt@other$loc.metrics.flags) <- recalc.flags
gl.RC.filt <- gl.recalc.metrics(gl.RC.filt)


##########################################
#     Create function to calculate       #
# population genetics summary statistics #
##########################################

# For each locus, the following function calculates:  
  
## Allele frequencies (total and within populations)
## Mean observed heterozygosity (within populations)
## Mean expected heterozygosity (HS - over populations)
## HT (He over total using mean allele frequencies over pops)
## FIS
## FIT
## FST
## The proportion of populations significantly out of HWE (using a Bonferroni correction) 

### Note that calculations were based on GenAlEx formulae (Peakall et al. 2006; Peakall et al. 2012)

# A genlight object is returned, with these metrics appended to the loc.metrics slot. The population-level metrics are only calculated for the populations specified above.  

PopGenMetrics <- function(gl, pops){
  
  # Calculate SNP/ref allele frequency over total
  gl@other$loc.metrics$AlleleFreq1 <- gl.alf(gl)[, 1]
  gl@other$loc.metrics$AlleleFreq2 <- gl.alf(gl)[, 2]
  
  # Calculate sample number (excluding missing data) for each locus, for pops above or equal to pop cutoff threshold
  SampleNo <- aggregate(gl[gl@pop %in% pops], list(gl[gl@pop %in% pops]@pop), function(x) {sum(!is.na(x))})
  
  # count number of heterozygotes (excluding missing data) for each locus, for pops above or equal to pop cutoff threshold
  HetCount <- aggregate(gl[gl@pop %in% pops], list(gl[gl@pop %in% pops]@pop), function(x) {sum(x[x=="1"], na.rm=TRUE)})
  
  # Calculate allele frequencies (excluding missing data) for each locus, for pops above or equal to pop cutoff threshold and append to genlight loc.metrics
  PopAlleleFreq1 <- aggregate(gl[gl@pop %in% pops], list(gl[gl@pop %in% pops]@pop), function(x) {gl.alf(x)[, 1]})
  
  PopAlleleFreq2 <- aggregate(gl[gl@pop %in% pops], list(gl[gl@pop %in% pops]@pop), function(x) {gl.alf(x)[, 2]})
  
  # Calculate observed and expected heterozygosity (excluding missing data) for each locus, for pops above or equal to pop cutoff threshold
  PopHo <- HetCount[,-1]/SampleNo[,-1]
  PopHe <- 1-((PopAlleleFreq1[,-1]^2) + (PopAlleleFreq2[,-1]^2))
  
  # Calculate HS (Mean He over pops) and append to locus metrics in genlight object
  gl$other$loc.metrics$HS <- colMeans(PopHe, na.rm=TRUE)
  
  # Calculate mean Ho over pops and append to locus metrics in genlight object
  gl$other$loc.metrics$meanHo <- colMeans(PopHo, na.rm=TRUE)
  
  # Calculate HT (He over total using mean allele freq over pops) and append to locus metrics in genlight object
  # First append mean allele frequencies over pops to genlight loc.metrics
  gl$other$loc.metrics$PopAlleleFreq1 <- colMeans(PopAlleleFreq1[,-1])^2
  gl$other$loc.metrics$PopAlleleFreq2 <- colMeans(PopAlleleFreq2[,-1])^2
  
  # Then replace NaN with 0 in pop allele frequency df, so matches the way that Genalex calculates HT
  is.nan.df <- function(x) {
    do.call(cbind, lapply(x, is.nan))
  }
  
  PopAlleleFreq1[is.nan.df(PopAlleleFreq1)] <- 0
  PopAlleleFreq2[is.nan.df(PopAlleleFreq2)] <- 0
  
  gl$other$loc.metrics$HT <- 1 - ((colMeans(PopAlleleFreq1[,-1])^2) + (colMeans(PopAlleleFreq2[,-1])^2))
  
  # Calculate FIS and append to locus metrics in genlight object
  gl$other$loc.metrics$FIS <- 1 - (gl$other$loc.metrics$meanHo / gl$other$loc.metrics$HS)
  
  # Calculate FIT and append to locus metrics in genlight object
  gl$other$loc.metrics$FIT <-(gl$other$loc.metrics$HT - gl$other$loc.metrics$meanHo) / gl$other$loc.metrics$HT
  
  # Calculate FST and append to locus metrics in genlight object
  gl$other$loc.metrics$FST <- (gl$other$loc.metrics$HT - gl$other$loc.metrics$HS) / gl$other$loc.metrics$HT
  
  # Calculate HWE (with Bonferroni Corrected significance) and append proportion of pops sig out of HWE to locus metrics in genlight object
  # Calculate HWE with Bonferroni correction for multiple testing usng dartR function only over pops with sample size over threshold
  HWE_Raw <- gl.report.hwe(x = gl[gl@pop %in% pops, ], 
                           method_sig = "ChiSquare", 
                           plot.out = FALSE,
                           multi_comp = TRUE, 
                           multi_comp_method = "bonferroni")
  HWE_Raw$Locus <- as.character(HWE_Raw$Locus)
  
  # Replace "ns" (non-significant) with a zero, and everything else (i.e. significant loci "*")
  # with a one. Note that I am using the bonferroni significance, rather than just the 0.05 significance
  # which is why there are ns results.
  HWE_Raw$BonSig_1.0 <- ifelse(HWE_Raw$Sig.adj == "ns", 0, 1)
  
  # Calculate the proportion of populations where that locus is out of HWE
  HWE.sig.pop.count <- HWE_Raw %>% 
    group_by(Locus) %>% 
    summarise(PropPopsSig.HWE = sum(BonSig_1.0)/length(pops))

  # Create a new data frame with all loci (rather than just the subset out of HWE)
  Locus_HWE <- data.frame(Locus = locNames(gl))
  Locus_HWE$Locus <- as.character(Locus_HWE$Locus)
  Locus_HWE <- left_join(Locus_HWE, HWE.sig.pop.count, by = "Locus")
  
  # Replace NA (i.e. those that were not out of HWE) with a 0
  Locus_HWE$PropPopsSig.HWE[is.na(Locus_HWE$PropPopsSig.HWE)] <- 0
  
  # Check locus order matches so can append to genlight
  sum(!Locus_HWE$Locus == locNames(gl))
  
  # Append to genlight
  gl@other$loc.metrics$PropPopsSig.HWE <- Locus_HWE$PropPopsSig.HWE
  
  return(gl)
  
}


# Calculate summary stats for raw versus filtered SNP data sets 
gl.RC.raw.PG <- PopGenMetrics(gl = gl.RC.raw, pops = PopKeep)
gl.RC.filt.PG <- PopGenMetrics(gl = gl.RC.filt, pops = PopKeep)


# How does individual read count impact population genetics statistics (across loci)?

# Note that mean read count and call rate thresholds were based on default values purely for exploring the impact of different filtering strategies. Real thresholds will be assigned later during the filtering steps.

# Use the following cutoffs for individual read count: 
## min: 
MinRC
## max: 
MaxRC

# Visualise the different data sets:  
# 1) Raw SNPs vs. filtered on individual read count 
# 2) Filtered on mean total read count (min= 20, max= 200) vs. filtered on individual read count
# 3) Filtered on mean total read count and call rate (min= 20, max= 200, CR = 90%) vs. filtered on individual read count and call rate (CR = 90%)
# Visualise how the distribution of population genetics metrics change with different read count filters.  

# Also visualise allele balance, the proportion of reference reads compared to total reads across heterozygotes (expected = 0.5), as large deviations may indicate false heterozygotes due to coverage effects, multilocus contigs or other artifacts O'Leary et al. 2018.  

# Compare raw versus filtered RC summary stats:
RawRCSummStats <- gl.RC.raw.PG@other$loc.metrics
RawRCSummStats <- as.data.frame(add_column(RawRCSummStats, ReadCount = "a", .before = 1))
RawRCSummStats.AvFilt <- RawRCSummStats[!is.na(RawRCSummStats$RC_MeanTot) & (RawRCSummStats$RC_MeanTot >= 20) & (RawRCSummStats$RC_MeanTot <= 200), ]
RawRCSummStats.Av.CR.Filt <- RawRCSummStats.AvFilt[RawRCSummStats.AvFilt$CallRate >= 0.9, ]

FiltRCSummStats <- gl.RC.filt.PG@other$loc.metrics
FiltRCSummStats <- add_column(FiltRCSummStats, ReadCount = "b", .before = 1)
FiltRCSummStats.CR.Filt <- FiltRCSummStats[FiltRCSummStats$CallRate >= 0.9, ]

# Merge data frames, add column with filtering variable specified
RawRCSummStats <- add_column(RawRCSummStats, DF = "a", .before = 1)
RawRCSummStats.AvFilt <- add_column(RawRCSummStats.AvFilt, DF = "b", .before = 1)
RawRCSummStats.Av.CR.Filt <- add_column(RawRCSummStats.Av.CR.Filt, DF = "c", .before = 1)

FiltRCSummStats <- add_column(FiltRCSummStats, DF = "a", .before = 1)
FiltRCSummStats.2 <- FiltRCSummStats
FiltRCSummStats.2$DF[FiltRCSummStats.2$DF == "a"] <- "b"
FiltRCSummStats.CR.Filt <- add_column(FiltRCSummStats.CR.Filt, DF = "c", .before = 1)

RCSummStats <- rbind(RawRCSummStats, RawRCSummStats.AvFilt, RawRCSummStats.Av.CR.Filt, FiltRCSummStats, FiltRCSummStats.2, FiltRCSummStats.CR.Filt)

# Reshape by metric
RCSummStats.reshape <- pivot_longer(data = RCSummStats[, c("DF", "ReadCount", "meanHo", "HS", "HT", "FIS", "FST")], cols = c("meanHo", "HS", "HT", "FIS", "FST"), names_to = "Metric")

RCSummStats.reshape$Metric[RCSummStats.reshape$Metric == "meanHo"] <- "a"
RCSummStats.reshape$Metric[RCSummStats.reshape$Metric == "HS"] <- "b"
RCSummStats.reshape$Metric[RCSummStats.reshape$Metric == "HT"] <- "c"
RCSummStats.reshape$Metric[RCSummStats.reshape$Metric == "FIS"] <- "d"
RCSummStats.reshape$Metric[RCSummStats.reshape$Metric == "FST"] <- "e"

RC.labs <- c("Unfiltered vs. Ind RC Filter", "Locus Mean RC Filter vs. Ind RC Filter", "Locus Mean RC Filter vs. Ind RC Filter (Call Rate >= 0.9)")
names(RC.labs) <- unique(RCSummStats.reshape$DF)


RC.Plot <- ggplot(RCSummStats.reshape, aes(x = Metric, y = value, fill = ReadCount)) +
  facet_wrap(~ DF, ncol = 1, strip.position = "top", labeller = labeller(DF = RC.labs)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.2) +
  geom_boxplot(colour = "black", size = 0.2, outlier.colour="black", outlier.size = 0.2) +
  scale_x_discrete(breaks = unique(RCSummStats.reshape$Metric), labels = c("Mean Ho", "HS", "HT", "FIS", "FST")) +
  scale_fill_manual(name = "Individual Read Count Filter:", labels = c("Raw", paste0("Min = ", MinRC, ", Max = ", MaxRC)), values = c("#00A087FF", "#91D1C2FF")) +
  labs(x = "Summary Statistic", y = "Value") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(face = "bold", size = 6),
        legend.position = "bottom")


# Prepare for Facet wrap in AB Plot
RCSummStats <- add_column(RCSummStats, RC_DF = paste(RCSummStats$ReadCount, RCSummStats$DF, sep = "_"), .before = 1)

AB.labs <- c("Unfiltered", "Locus Mean RC Filter", "Locus Mean RC Filter (CR >= 0.9)", "Ind RC Filter", "Ind RC Filter", "Ind RC Filter (CR >= 0.9)")
names(AB.labs) <- sort(unique(RCSummStats$RC_DF))

AB.plot <- ggplot(data = RCSummStats, aes(x = AB)) +
  geom_histogram(bins = 100, colour = "black", fill = "#7570B3", size = 0.2) +
  facet_wrap(~ RC_DF, ncol = 2, dir = "v", strip.position = "top", labeller = labeller(RC_DF = AB.labs)) +
  geom_vline(xintercept = c(0.2, 0.8), color = "darkblue", linetype = "dotted", size = 0.4) +
  geom_vline(xintercept = 0.5, color= "red", linetype = "dashed", size = 0.4) +
  labs(x = "\nAllele Balance", y = "No. Loci\n") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5))

# View
grid.arrange(RC.Plot, AB.plot, ncol = 2)

# Output figure
jpeg(filename = paste0(outpath, Outfile_Desc, "Ind.ReadCount.jpg"), units =  "cm", width = 20, height = 28, res = 300)
grid.arrange(ggarrange(RC.geno.Raw, RC.geno.Filt, ncol = 2, labels = c("a)", "b)"), font.label = list(size = 8, face = "plain")), ggarrange(RC.Plot, AB.plot, ncol = 2, labels = c("c)", "d)"), font.label = list(size = 8, face = "plain")), heights = c(1/3, 2/3), top = textGrob("Read count summary statistics\n", gp = gpar(fontsize = 10, fontface = "bold")))
dev.off()

# The above figure demonstrates that filtering on mean read count is roughly equivalent to filtering on individual read count, when combined with a call rate filter (in terms of population genetic metrics).


#############################
# Remove failed individuals #
#############################

# After exploring read count, go through the visualisation/filtering process. The first step involves removing individuals that have failed, i.e. individuals that have such high levels of missing data that they will lower the quality of the entire data set if kept in. Retain as many individuals as possible, as they may be usable after downstream locus filtering. However, it's important to remove any individuals that will have too much missing data at the start, because metrics will need to be recalculated every time samples are removed.  

# Visualise missing data for individuals

# Map showing location of samples, coloured by individual call rate

# Read in map
StudyMap <- readOGR(paste0("../Rasters_Shapefiles/", Map.shp))

# Process map for plotting
StudyMap@data$id <- rownames(StudyMap@data)
StudyMap.points = fortify(StudyMap, region = "id")
StudyMap.df <- left_join(StudyMap.points, StudyMap@data, by="id")

# Match individual call rates with coordinates
IndCR <- data.frame(ID = rownames(as.matrix(gl.RC.raw)), IndCR = 1 - rowSums(is.na(as.matrix(gl.RC.raw)))/nLoc(gl.RC.raw))
IndCRxy <- left_join(x = cbind(ID = rownames(gl.RC.raw@other$latlong), gl.RC.raw@other$latlong), IndCR, by = "ID")

# Plot
Samp.Map <- ggplot() +
  geom_polygon(data = StudyMap.df, aes(x = long, y = lat, group = group, fill = get(MapFill))) +
  scale_fill_manual(values = gray.colors(4)) +
  geom_path(data = StudyMap.df, aes(x = long, y = lat, group = group), colour = "transparent", size = 0.2) +
  scalebar(location = "bottomleft", dist_unit = distU, transform = TRUE, model = CoordRef,
           x.min = xmin, x.max = xmax, y.min = ymin, y.max = ymax, 
           dist = dist, st.dist = 0.04, st.size = 1.5, height = 0.03, border.size = 0.15) +
  north(location = "topleft", scale = 0.1, symbol = 4,
        x.min = xmin, x.max = xmax, y.min = ymin, y.max = ymax) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 10, margin = margin(0,16,0,0), face ="bold"),
        axis.title.x = element_text(size = 10, margin = margin(16,0,0,0), face ="bold"),
        axis.ticks = element_line(size = 0.2)) +
guides(fill = "none") +
geom_jitter(data = IndCRxy, aes(x = lon, y = lat, colour = IndCR), size = 1, width = 0.05, height = 0.05) +
  labs(colour = "Individual\ncall rate") +
  scale_color_viridis(option = "C", direction = -1) +
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))

# Bar plot (call rate threshold vs number of individuals filtered out)
ind.call.rate <- data.frame(CR = 1 - rowSums(is.na(as.matrix(gl.RC.raw)))/nLoc(gl.RC.raw))
ind.call.rate <- data.frame(Threshold = seq(0,1,0.05), 
                            Filtered = cumsum(table(cut(ind.call.rate$CR, 
                                                        seq(-0.05, 1, 0.05)))))

IndCR.plot <- ggplot(data = ind.call.rate, aes(x = Threshold, y = Filtered)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Filtered), vjust = -0.3, color = "black", size = 2) +
  labs(x = "Individual call rate threshold", y = "No. Individuals removed") +
  scale_x_continuous(breaks = seq(0,1,0.05)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text = element_text(size = 4.5))


# View smear plot
smear.plot <- gl.smearplot(x = gl.RC.raw, group_pop = TRUE, plot_colors = c("#482677FF", "#20A387FF", "#FDE725FF", "white"))
smear.plot

# View samp plot
ggarrange(IndCR.plot, Samp.Map, ncol = 1, heights = c(1/3, 2/3), labels = c("a)", "b)"), font.label = list(size= 10, face = "bold"))

# Output plots
jpeg(filename = paste0(outpath, Outfile_Desc, ".Ind.Pop.jpg"), units =  "cm", width = 28, height = 20, res = 300)
ggarrange(smear.plot, labels = "a)", font.label = list(size = 8, face = "plain"), ggarrange(IndCR.plot, popNo.plot, ncol = 2, labels = c("b)", "c)"), font.label = list(size = 8, face = "plain"), widths = c(2/5, 3/5)), ncol = 1)
dev.off()

jpeg(filename = paste0(outpath, Outfile_Desc, "SampleMap.IndCR.jpg"), units =  "cm", width = 20, height = 13, res = 300)
Samp.Map
dev.off()

# Inspect above Figure and smear plot to decide if any individuals should be removed initially, by deciding on an individual missing data threshold. The aim here is just to remove any samples that are obvious mis-ids (for example, a different species where only a small number of loci have amplified that are obviously fixed for a different allele) or samples that have mainly failed (where there will never be enough loci for them to be usable). Use the map to decide if any individuals are particularly important (i.e. they are the only sample in a priority location) so that the call rate threshold can be tailored to retain these samples if possible.  

# Filter failed individuals from data set

IndCallRate_Cutoff <- 0.64

# Remove failed samples from SNP data frame and recalculate metadata
gl.RC.DropInds <- gl.filter.callrate(x = gl.RC.raw, method = "ind", recalc = TRUE, mono.rm = FALSE, threshold = IndCallRate_Cutoff)

# Use a default locus filtering value to check if any samples drop below an acceptable level of missing data
# Carry out the locus filtering properly later on. However, be sure to get the individual filtering right here so that the population genetic summary statistics only need to be calculated once (because they take a long time and have to be recalculated every time individuals are removed)
gl.checkCR <- gl.filter.callrate(gl.RC.DropInds, method = "loc", threshold = 0.85)
IndCR.2 <- data.frame(ID = rownames(as.matrix(gl.checkCR)), IndCR = 1 - rowSums(is.na(as.matrix(gl.checkCR)))/nLoc(gl.checkCR))


# Define individual missing data threshold. This depends on the data set! Data sets with lots of SNPs can afford a lower threshold (e.g. in a data set of 100k snps, 10% still leaves 10k SNPs to play with). 

# Unfortunately, many individuals did not sequence well in this data set (likely due to DNA quality/quantity). This species also yielded low numbers of loci (with only ~12k SNPs called in the unfiltered data set). This means that the individual call rate threshold will have to be quite stringent so that the data set has enough power to detect the (presumed) low genetic structure present across the Pilbara. 

# I've decided on a threshold of:
IndCallRate_Cutoff
# which has removed a massive:
cat(nInd(gl.RC.raw) - nInd(gl.RC.DropInds), "individuals")
# If data were (hypothetically) filtered on a locus call rate of at least 85%, 
cat(sum(IndCR.2$IndCR < 0.85), "individuals")
# would drop below an 85% individual call rate threshold. However, 85% of:
nLoc(gl.RC.DropInds)
# still leaves:
cat(round(nLoc(gl.RC.DropInds)*0.85, 0), "loci")
# to play with, so I'm hopeful that I'll end up with enough overlapping loci after further downstream filtering.  


# Check population cutoff
# Check if the number of individuals in each population have dropped below the sample size cutoff now that failed individuals have been filtered out of the data set. If they have changed, decide whether or not to include these populations when calculating the population genetic metrics for each locus.

# Look at the raw genlight, as both have the same metadata
PopNo <- data.frame(table(pop(gl.RC.DropInds)))

popNo.plot <- ggplot(data = PopNo, aes(x = Var1, y = Freq, fill = Freq)) +
  geom_bar(stat = "identity", colour = "black", size = 0.2) +
  scale_fill_viridis(discrete = F, option = "D", direction = -1) +
  labs(x = "Population", y = "No. Individuals", fill = "Pop Size") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 6))

# View
popNo.plot

# Check if any of the pops have dropped below threshold sample size number
PopKeep.IndFilt <- names(table(gl.RC.DropInds@pop) >= PopCutOff)[table(gl.RC.DropInds@pop) >= PopCutOff]

if(!length(PopKeep.IndFilt) == length(PopKeep)){
  PopKeep <- PopKeep.IndFilt
}

# Here, one population has dropped below the sample size cutoff of 10 (so will be removed from population genetic summary statistic calculations).  


####################################
# Add extra chromosome information #
####################################

# Some DArTseq data sets are blasted against a reference genome (if one has been provided). The information that comes back is often in the form of contig IDs, rather than the actual chromosome. Searching the genome on NCBI often provides extra information relating these contig IDs to the chromosome (depending on how good the reference is). If available, add this information to explore how many SNPs there are per Chromosome.

if(BLAST == "Y"){

ChromInfo.ncbi <- read.csv(file = paste0(datpath, ChromInfo.file), stringsAsFactors = FALSE)
  
## Rename column so I can join the two data frames by this variable
colnames(ChromInfo.ncbi)[which(colnames(ChromInfo.ncbi) == ContigCol_ncbi)] <- ContigCol_dart

## Create a new data frame with all loci
gl.RC.DropInds@other$loc.metrics[, c("Locus", ContigCol_dart, ContigPos_dart)] %>% mutate_if(is.factor, as.character) -> Locus_Chrom

## Join data frames for the chromosome number that corresponds to the contig ID
Locus_Chrom <- left_join(Locus_Chrom, ChromInfo.ncbi, by = ContigCol_dart)

## Add the contig position to the chromosome position to get the actual position along the chromosome
Locus_Chrom$TrueChromPos <- Locus_Chrom[, ContigPos_dart] + Locus_Chrom[, ChromPos_ncbi]

## Rename NA as "UnMapped
Locus_Chrom$Chromosome[is.na(Locus_Chrom$Chromosome)] <- "UnMapped"

## Check locus order matches so can append to genlight
sum(!Locus_Chrom$Locus == locNames(gl.RC.DropInds))

## Append to genlight
gl.RC.DropInds@other$loc.metrics$Chromosome <- Locus_Chrom$Chromosome
gl.RC.DropInds@other$loc.metrics$TrueChromPos <- Locus_Chrom$TrueChromPos
}

########################################
# Calculate population genetic metrics #
########################################

# Run the function to calculate population genetic metrics across the data set filtered on failed individuals.
gl.RC.DI.PG <- PopGenMetrics(gl = gl.RC.DropInds, pops = PopKeep)


#######################################
# Calculate number of private alleles #
#######################################

# How does the minor allele frequency affect the number of private alleles. This is a useful way of visualising the interplay between real differences among populations (which may be important for detecting migrants, for example), and the prevalence of error in the data set (i.e. a low frequency SNP may indicate sequencing error). For now, just calculate private alleles across the raw data set, and dig into this after filtering. Note that this script calculates the mean number of private alleles between pairwise population comparisons.

# Use dartR script to calculate number of private alleles between each population pair
PrivAllele_raw <- gl.report.pa(gl.RC.DI.PG[gl.RC.DI.PG@pop %in% PopKeep, ], plot.out = FALSE)

# Drop pop factors that have been filtered out using drop levels and only keep relevant columns (rename so can combine)
priv1_raw <- PrivAllele_raw[ , c("pop1", "priv1")]
priv2_raw <- PrivAllele_raw[ , c("pop2", "priv2")]
colnames(priv1_raw) <- c("pop", "priv")
colnames(priv2_raw) <- c("pop", "priv")

# Combine into one data set so all pop comparisons are in one column
priv_raw <- rbind(priv1_raw, priv2_raw)
# Will plot this later

# Output all summary statistics for pre-filtered data set
write.csv(gl.RC.DI.PG@other$loc.metrics, paste0(outpath, Outfile_Desc, ".RawSumStats.csv"), row.names = F)


#################################
# Visualise and filter raw data #
#################################

# Thoroughly interrogate the raw data set to tailor filtering thresholds to the study. The first step is to look at the patterns in the unfiltered, raw data set using a PCoA.

# Specify column for colouring/grouping PCoA
GroupCol_PCoA <- "IBRA_SubRegion"

# Run PCoA
gl.RC.DI.PG.pcoa <- gl.RC.DI.PG
pop(gl.RC.DI.PG.pcoa) <- gl.RC.DI.PG.pcoa$other$ind.metrics[, GroupCol_PCoA]
pc <- gl.pcoa(gl.RC.DI.PG.pcoa)
gl.pcoa.plot(pc, gl.RC.DI.PG.pcoa, 
             scale = TRUE, 
             pt.size = 1,
             pop.labels = "pop", 
             save2tmp = TRUE)

# Figure out which temp file is the PCoA plot and save it
files_tempdir_plot <- list.files(tempdir())[which(str_match(list.files(tempdir()), "Plot") == "Plot")]
plot.last <- files_tempdir_plot[order(file.info(paste0(tempdir(), "/", files_tempdir_plot))$atime)][length(files_tempdir_plot[order(file.info(paste0(tempdir(), "/", files_tempdir_plot))$atime)])]
plot.no <- which(list.files(tempdir()) == plot.last)

Pcoa.filt <- gl.print.reports(print_report = plot.no)$plot
Pcoa.raw <- gl.print.reports(print_report = plot.no)$plot

# Plot
Pcoa.raw <- Pcoa.raw +
  labs(title = "Raw PCoA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", 
                                    colour = "black", 
                                    fill = "transparent",
                                    size = 0.5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(colour="black",
                                  size=12, 
                                  face = "bold", 
                                  hjust = 0.5),
        axis.title = element_text(colour="black",
                                  size=10, 
                                  face = "bold"),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(colour="black", 
                                 size=8, 
                                 face= "plain"),
        axis.text.y = element_text(colour="black", 
                                   size=8, 
                                   face= "plain"),
        legend.position = "none")

# Add map coloured by PCoA pop category
IndPop <- cbind(gl.RC.DI.PG@other$ind.metrics[, GroupCol_PCoA], gl.RC.DI.PG@other$latlong)
colnames(IndPop)[1] <- "PCoA_Grp"

Samp.Map.pop <- ggplot() +
  geom_polygon(data = StudyMap.df, aes(x = long, y = lat, group = group, fill = get(MapFill))) +
  scale_fill_manual(values = gray.colors(4)) +
  geom_path(data = StudyMap.df, aes(x = long, y = lat, group = group), colour = "transparent", size = 0.2) +
  scalebar(location = "bottomleft", dist_unit = distU, transform = TRUE, model = CoordRef,
           x.min = xmin, x.max = xmax, y.min = ymin, y.max = ymax, 
           dist = dist, st.dist = 0.04, st.size = 1.5, height = 0.03, border.size = 0.15) +
  north(location = "topleft", scale = 0.1, symbol = 4,
        x.min = xmin, x.max = xmax, y.min = ymin, y.max = ymax) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", 
                                    fill = "transparent", size = 0.5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 10, margin = margin(0,16,0,0), face ="bold"),
        axis.title.x = element_text(size = 10, margin = margin(16,0,0,0), face ="bold"),
        axis.ticks = element_line(size = 0.2)) +
guides(fill = "none") +
geom_jitter(data = IndPop, aes(x = lon, y = lat, colour = PCoA_Grp), size = 1, width = 0.05, height = 0.05) +
  labs(colour = "IBRA\nSub-Region") +
  theme(legend.position = "right",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8, hjust = 0.5),
        legend.key = element_rect(fill = "white"))

# Put PCoA plot and Map together
PCoA.raw.plot <- ggarrange(Pcoa.raw, Samp.Map.pop, nrow = 2, labels = c("a)", "b)"), font.label = list(size= 10, face = "bold"))

# View 
PCoA.raw.plot

# Output result
jpeg(filename = paste0(outpath, Outfile_Desc, "Raw.PCoA_Map.jpg"), units =  "cm", width = 20, height = 28, res = 300)
PCoA.raw.plot
dev.off()


# There are some strong regional groupings in the PCoA, with 5% of the variation described by PC1 (suggesting that there is a high level of genetic differentiation between dolphin island and the rest of the Pilbara). Note that the patterns shown in the unfiltered PCoA often do not change substantially after filtering (which is a good indication that the filtering hasn't biased results, but has just cleaned the data set so that inferences are based on high quality, true loci).


######################
# Summary statistics #
######################

# How do quality metrics impact population genetic summary statistics? The goal is to avoid biasing 'real' results by artifacts in the data (e.g. quality). Therefore, choose filters based on a threshold where data quality scores are no longer driving/correlated with summary statistics.  

# Histograms

# Save loc.metrics as data frame
RawRCSummStats <- gl.RC.DI.PG@other$loc.metrics

# Add a column specifying if locus has blasted to the reference (mapped/unmapped)
if(BLAST == "Y") {
  RawRCSummStats$Chrom.Mapped <- ifelse(RawRCSummStats$Chromosome == "UnMapped", "UnMapped", "Mapped")
  Chrom_Mapped <- RawRCSummStats[!RawRCSummStats$Chromosome == "UnMapped", ]

# Reshape by metric for facet_wrap in ggplot
prefiltDF.hist <- RawRCSummStats[, c("Locus", "Chrom.Mapped", "CallRate", "RepAvg", "RC_MeanTot", "maf", "meanHo", "HS", "FIS", "FST")]
prefiltDF.hist <- pivot_longer(data = prefiltDF.hist, cols = colnames(select_if(prefiltDF.hist, is.double)), names_to = "MetaDat_Var")
} else{
# Reshape by metric for facet_wrap in ggplot
prefiltDF.hist <- RawRCSummStats[, c("Locus", "CallRate", "RepAvg", "RC_MeanTot", "maf", "meanHo", "HS", "FIS", "FST")]
prefiltDF.hist <- pivot_longer(data = prefiltDF.hist, cols = colnames(select_if(prefiltDF.hist, is.double)), names_to = "MetaDat_Var")

}

# Rename variables so that they are in the right order in the plot
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "CallRate"] <- "a"
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "RepAvg"] <- "b"
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "RC_MeanTot"] <- "c"
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "maf"] <- "d"
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "meanHo"] <- "e"
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "HS"] <- "f"
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "FIS"] <- "g"
prefiltDF.hist$MetaDat_Var[prefiltDF.hist$MetaDat_Var == "FST"] <- "h"

hist.labs <- c("Call Rate", "Mean Repeatability", "Mean Total Read Count", "MAF", "Mean Ho", "Mean He", "FIS", "FST")
names(hist.labs) <- letters[1:8]

# Plot histograms- total data set
Prefilt.hist.plot <- ggplot(data = prefiltDF.hist, aes(x = value, fill = MetaDat_Var)) +
  geom_histogram(bins = 50, colour = "black", size = 0.2) +
  facet_wrap(~ MetaDat_Var, ncol = 2, scales = "free", labeller = labeller(MetaDat_Var = hist.labs)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "\nValue", y = "No. Loci\n") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "None")


# Mean Ho, HS, HT plots

# Prepare individual plots for grid arrange
mHoHS.CR_Plot <- ggplot(RawRCSummStats, aes(x = HS, y = meanHo, col = CallRate)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Call Rate", y= " \n ", x = " ") +
  scale_colour_viridis(option = "D", direction = -1, begin = 0.2) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

HS.HT.CR_Plot <- ggplot(RawRCSummStats, aes(x = HT, y = HS, col = CallRate)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y= " \n ", x = " ") +
  scale_colour_viridis(option = "D", direction = -1, begin = 0.2) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

scatter.CR <- ggarrange(mHoHS.CR_Plot, HS.HT.CR_Plot, ncol = 2, labels = c("a)", "b)"), font.label = list(size = 8, face = "plain"))

mHoHS.RC_Plot <- ggplot(RawRCSummStats, aes(x = HS, y = meanHo, col = RC_MeanTot)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Mean\nRead Count", y= " \n ", x = " ") +
  scale_colour_viridis(option = "A", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

HS.HT.RC_Plot <- ggplot(RawRCSummStats, aes(x = HT, y = HS, col = RC_MeanTot)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y= " \n ", x = " ") +
  scale_colour_viridis(option = "A", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

scatter.RC <- ggarrange(mHoHS.RC_Plot, HS.HT.RC_Plot, ncol = 2, labels = c("c)", "d)"), font.label = list(size = 8, face = "plain"))

mHoHS.Rep_Plot <- ggplot(RawRCSummStats, aes(x = HS, y = meanHo, col = RepAvg)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Mean\nRepeatability", y = " \nMean Ho", x = " ") +
  scale_colour_viridis(option = "B", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.title.y = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

HS.HT.Rep_Plot <- ggplot(RawRCSummStats, aes(x = HT, y = HS, col = RepAvg)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y = " \nMean He (HS)", x = " ") +
  scale_colour_viridis(option = "B", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.title.y = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

scatter.Rep <- ggarrange(mHoHS.Rep_Plot, HS.HT.Rep_Plot, ncol = 2, labels = c("e)", "f)"), font.label = list(size = 8, face = "plain"))

mHoHS.HWE_Plot <- ggplot(RawRCSummStats, aes(x = HS, y = meanHo, col = PropPopsSig.HWE)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "HWE:\nProp. Pops Sig", y= " \n ", x = " ") +
  scale_colour_viridis(option = "C") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

HS.HT.HWE_Plot <- ggplot(RawRCSummStats, aes(x = HT, y = HS, col = PropPopsSig.HWE)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y = " \n ", x = " ") +
  scale_colour_viridis(option = "C") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

scatter.HWE <- ggarrange(mHoHS.HWE_Plot, HS.HT.HWE_Plot, ncol = 2, labels = c("g)", "h)"), font.label = list(size = 8, face = "plain"))

mHoHS.FIS_Plot <- ggplot(RawRCSummStats, aes(x = HS, y = meanHo, col = FIS)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "FIS", y= " \n ", x = "Mean He (HS)") +
  scale_colour_viridis(option = "B") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_line(size = 0.2),
        axis.text.x = element_text(colour="black",size=5,face="plain"),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.09, 0.6),
        legend.key.size = unit(0.2, "cm"))

HS.HT.FST_Plot <- ggplot(RawRCSummStats, aes(x = HT, y = HS, col = FST)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "FST", y = " \n ", x = "HT") +
  scale_colour_viridis(option = "B") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_line(size = 0.2),
        axis.text.x = element_text(colour="black",size=5,face="plain"),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.09, 0.6),
        legend.key.size = unit(0.2, "cm"))

scatter.Fstats <- ggarrange(mHoHS.FIS_Plot, HS.HT.FST_Plot, ncol = 2, labels = c("i)", "j)"), font.label = list(size = 8, face = "plain"))


# Correlation Plot
if(BLAST == "Y"){
prefiltDF.cor <- RawRCSummStats[, c("CallRate", "RC_MeanTot", "RC_MeanRef", "RC_MeanSNP", "AB", "RepAvg", EVal.col, Cnt.col, "AvgPIC", "maf", "HS", "meanHo", "HT", "FIS", "FST", "PropPopsSig.HWE")]
prefilt.cor <- cor(prefiltDF.cor, use="pairwise.complete", method = "spearman")

colnames(prefilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "E Val.", "Aln. Count", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
rownames(prefilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "E Val.", "Aln. Count", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
} else {
prefiltDF.cor <- RawRCSummStats[, c("CallRate", "RC_MeanTot", "RC_MeanRef", "RC_MeanSNP", "AB", "RepAvg", "AvgPIC", "maf", "HS", "meanHo", "HT", "FIS", "FST", "PropPopsSig.HWE")]
prefilt.cor <- cor(prefiltDF.cor, use="pairwise.complete", method = "spearman")

colnames(prefilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
rownames(prefilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
}

# View
Prefilt.hist.plot
mHoHSHT.plot <- ggarrange(scatter.CR, scatter.RC, scatter.Rep, scatter.HWE, scatter.Fstats, ncol = 1)
mHoHSHT.plot
corrplot.mixed(prefilt.cor, upper = "circle", lower = "number", order = "FPC", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = 0.5, number.font = 1, cl.cex = 0.5, lower.col = "black", upper.col = c(viridis(n = 10, option = "D", begin = 0.1, end = 0.9), viridis(n = 10, option = "A", direction = -1, begin = 0.3)), outline = TRUE)

# Output figures
jpeg(filename = paste0(outpath, Outfile_Desc, "Corr.Raw.jpg"), units =  "cm", width = 20, height = 20, res = 300)
corrplot.mixed(prefilt.cor, upper = "circle", lower = "number", order = "FPC", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = 0.5, number.font = 1, cl.cex = 0.5, lower.col = "black", upper.col = c(viridis(n = 10, option = "D", begin = 0.1, end = 0.9), viridis(n = 10, option = "A", direction = -1, begin = 0.3)), outline = TRUE)
dev.off()

jpeg(filename = paste0(outpath, Outfile_Desc, "SummStat.Raw.jpg"), units =  "cm", width = 28, height = 20, res = 300)
ggarrange(Prefilt.hist.plot, mHoHSHT.plot, ncol = 2)
dev.off()


# The correlation plot shows that call rate is strongly correlated with read count. This is not surprising, as a low read count results in missing data. Call rate is also strongly correlated with FIS, FST, HWE and heterozygosity (which can also be seen in the scatter plots). 

# It is not particularly interesting to look into correlations between the different population genetic metrics, because it makes sense for the different heterozygosity measures and read count metrics to be correlated. Instead, look for correlations between quality metrics (read count, call rate, repeatability) and biological inferences (heterozygosity, FIS, FST, HWE). In this way, the data set can be filtered based on justifiable thresholds that correspond to data quality rather than biological hypotheses/expectations. After filtering on these metrics, see if any strange patterns in the data set disappear by monitoring the shape of the distributions in the histograms and by checking if the quality and population genetic metrics are no longer correlated. As an interesting side note, here, mean repeatability appears to be correlated with heterozygosity and MAF.  

# How are SNPs distributed among/along chromosomes? 

# If a reference genome was provided, explore how SNPs are mapped to the reference by visualising how SNPs are distributed among the different chromosomes.  

if(BLAST == "Y"){
# Plot distribution of SNPs along chromosomes (if you have this info)
ChromDensity <- ggplot(RawRCSummStats, aes(x = Chromosome, fill = Chromosome)) + 
                geom_bar(show.legend = FALSE, col = "black", size = 0.2) +
                scale_fill_viridis(discrete = TRUE, option = "D") +
                xlab("Chromosome") + 
                ylab("Number of SNPs") + 
                theme(panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype = "solid", 
                                                  colour = "black", 
                                                  fill = "transparent", 
                                                  size = 0.2),
                      panel.grid.minor = element_blank(), 
                      panel.grid.major = element_blank(),       
                      axis.title.x = element_text(face = "bold", size = 6),
                      axis.title.y = element_text(face = "bold", size = 6),
                      axis.ticks = element_line(size = 0.2),
                      axis.text.x = element_text(size = 5),
                      axis.text.y = element_text(size = 5),
                      legend.position = "None")

SNPDensity <- ggplot(Chrom_Mapped, aes(fill = Chromosome)) + 
    geom_histogram(aes(x = TrueChromPos), show.legend = FALSE, col = "black", size = 0.2) + 
    facet_wrap( ~ Chromosome, ncol = 2, scales = "free_x") +
    scale_fill_viridis(discrete = TRUE, option = "D") +
    xlab("Position on Chromosome") + 
    ylab("SNP density") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "None")
  
# Output as jpeg
jpeg(filename = paste0(outpath, Outfile_Desc, "Chrom.jpg"), units =  "cm", width = 20, height = 28, res = 300)
grid.arrange(ChromDensity, SNPDensity, nrow = 2,  top = textGrob("Pre-Filtering SNP distributions on Chromosomes", gp = gpar(fontsize = 10, fontface = "bold")), heights = c(1/4, 3/4))
dev.off()
}

# View
if(BLAST == "Y"){
  grid.arrange(ChromDensity, SNPDensity, nrow = 2,  top = textGrob("Pre-Filtering SNP distributions on Chromosomes", gp = gpar(fontsize = 10, fontface = "bold")), heights = c(1/5, 4/5))
}

# In this case:
cat(round((nrow(RawRCSummStats[RawRCSummStats$Chromosome == "UnMapped",])/nrow(RawRCSummStats))*100, 0), "%")
# of SNPs mapped to the reference genome. This suggests that the Tasmanian devil genome is fairly divergent from D. hallucatus, although still seems to make an adequate reference. In addition, some loci have BLASTED to sex chromosomes and unknown contigs.  

# Mapping quality
# Visualise the distributions of the different population genetic summary statistics and quality/frequency metrics to see if there's a difference between SNPs that were mapped to the reference genome, versus those that were not. There might be a difference if the reference genome is closely related (in which case, unmapped reads may represent poorer quality SNPs), or if the reference is very distant (in which case, mapped reads may be less variable). In either situation, it may be useful to filter based on whether a SNP was mapped to the reference or not (although the former may be more justifiable than the latter).

# Plot histograms- mapped data set
if(BLAST == "Y"){
  prefiltDF.hist.mapped <- prefiltDF.hist[prefiltDF.hist$Chrom.Mapped == "Mapped", ]
  
  Prefilt.hist.M.plot <- ggplot(data = prefiltDF.hist.mapped, 
                                aes(x = value, fill = MetaDat_Var)) +
    geom_histogram(bins = 50, colour = "black", size = 0.2) +
    facet_wrap(~ MetaDat_Var, 
               ncol = 2, 
               scales = "free", 
               labeller = labeller(MetaDat_Var = hist.labs)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "\nValue", y = "No. Loci\n", title = "Mapped SNPs") +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", 
                                      colour = "black", 
                                      fill = "transparent", 
                                      size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),       
          strip.text = element_text(face = "bold", size = 6),
          strip.background = element_rect(linetype = "solid", 
                                          colour = "black", 
                                          fill = "white", 
                                          size = 0.2),
          plot.title = element_text(face = "bold", 
                                    size = 8, 
                                    hjust = 0.5, 
                                    colour = "DarkBlue"),
          axis.title.x = element_text(face = "bold", size = 6),
          axis.title.y = element_text(face = "bold", size = 6),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "None")
  
  
  # Plot histograms- unmapped data set
  prefiltDF.hist.unmapped <- prefiltDF.hist[prefiltDF.hist$Chrom.Mapped == "UnMapped", ]
  
  Prefilt.hist.UM.plot <- ggplot(data = prefiltDF.hist.unmapped, 
                                 aes(x = value, fill = MetaDat_Var)) +
    geom_histogram(bins = 50, colour = "black", size = 0.2) +
    facet_wrap(~ MetaDat_Var, ncol = 2, 
               scales = "free", 
               labeller = labeller(MetaDat_Var = hist.labs)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "\nValue", y = "No. Loci\n", title = "Un-Mapped SNPs") +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", 
                                      colour = "black", 
                                      fill = "transparent", 
                                      size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),       
          strip.text = element_text(face = "bold", size = 6),
          strip.background = element_rect(linetype = "solid", 
                                          colour = "black", 
                                          fill = "white", 
                                          size = 0.2),
          plot.title = element_text(face = "bold", 
                                    size = 8, 
                                    hjust = 0.5, 
                                    colour = "DarkBlue"),
          axis.title.x = element_text(face = "bold", size = 6),
          axis.title.y = element_text(face = "bold", size = 6),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.position = "None")
  
  # Output figures
  jpeg(filename = paste0(outpath, Outfile_Desc, ".RawHistograms.Mapped-UnMapped.jpg"), units =  "cm", width = 28, height = 20, res = 300)
  
  grid.arrange(Prefilt.hist.M.plot + theme(plot.background = element_rect(colour = "black", size = 0.5)), Prefilt.hist.UM.plot + theme(plot.background = element_rect(colour = "black", size = 0.5)), ncol = 2)
  
  dev.off()
}

# View
if(BLAST == "Y"){
  grid.arrange(Prefilt.hist.M.plot + theme(plot.background = element_rect(colour = "black", size = 0.5)), Prefilt.hist.UM.plot + theme(plot.background = element_rect(colour = "black", size = 0.5)), ncol = 2)
}


# Looks like there's no real difference between mapped versus unmapped reads in terms of quality and population genetic metrics. Therefore, it isn't worth filtering on mapped/unmapped loci.  


#####################
# Quality filtering #
#####################

# Now that the data has been inspected, it's time to start filtering.

CR.report <- data.frame(CallRate = seq(0,1,0.05), Filtered = cumsum(table(cut(gl.RC.DI.PG@other$loc.metrics$CallRate,seq(-0.05, 1, 0.05)))))

CR.tab <- flextable(data = CR.report)
CR.tab <- autofit(CR.tab)
CR.tab <- bg(CR.tab, bg = "lightgrey", part = "header") 
CR.tab <- bold(CR.tab, part = "header")
CR.tab <- align(CR.tab, align = "left")
CR.tab <- border_outer(CR.tab, part = "header", border = fp_border(color = "black",style = "solid", width = 0.5))
CR.tab <- border_outer(CR.tab, part = "body", border = fp_border(color = "black",style = "solid", width = 0.5))
CR.tab <- fontsize(CR.tab, size = 10)

# View
CR.tab


# Choose a call rate filtering threshold

# After inspecting raw data figures, choose a call rate filter that makes sense for the data. This will take some playing around and is always a balance between retaining the highest quality SNPs and maintaining a large enough data set for adequate statistical power. It will likely take a few runs to get the balance right. 

# In this case, I chose a call rate threshold of:
CallRate <- 0.95
# which, after filtering:
gl.CR.Filt <- gl.filter.callrate(gl.RC.DI.PG, method = "loc", threshold = CallRate, mono.rm = TRUE, recalc = TRUE)
# has reduced the data set down to:
cat(nLoc(gl.CR.Filt), "SNPs.") 


# Choose read count filtering thresholds:

# Mean read count appears to represent individual read counts well (see initial read count exploration). Use raw data figures to help decisions about minimum and maximum read counts (while keeping in mind that genotypes were strongly impacted by read count below/above individual read count thresholds in initial exploration). Low read counts could represent erroneous SNP calls (or loci that will have high levels of allelic dropout and therefore poor estimates of heterozygosity, etc.). Very high read counts often represent paralogues (two different regions of the genome, that have been collapsed as a single locus due to having similar sequences). 

# In this case, I chose a minimum/maximum mean read count of:
RClower <- 20
RCupper <- 200
# and a maximum of:
# which, after filtering:
gl.CR.RC.Filt <- gl.CR.Filt[, !is.na(gl.CR.Filt@other$loc.metrics$RC_MeanTot) & gl.CR.Filt@other$loc.metrics$RC_MeanTot > RClower & gl.CR.Filt@other$loc.metrics$RC_MeanTot < RCupper] # Subset genlight rather than use dart function, because I've recalculated mean call rate
gl.CR.RC.Filt@other$loc.metrics <- gl.CR.Filt@other$loc.metrics[!is.na(gl.CR.Filt@other$loc.metrics$RC_MeanTot) & gl.CR.Filt@other$loc.metrics$RC_MeanTot > RClower & gl.CR.Filt@other$loc.metrics$RC_MeanTot < RCupper, ] # Have to subset loc.metrics too (annoyingly it doesn't do it automatically)

# has reduced the data set down to 
cat(nLoc(gl.CR.RC.Filt), "SNPs.")



# Choose repeatability filter:

# It's important to be stringent with repeatability filters to avoid including poor quality loci. However, being too stringent may result in filtering out SNPs with high levels of heterozygosity (as heterozygote genotype calls may vary between replicates due to slight variations in ref/SNP reads, as suggested by the correlation between mean repeatability and heterozygosity in the raw correlation plot. 

# For this reason, I've chosen a mean repeatability threshold of:
RepFilt <- 0.95
# which, after filtering:
gl.CR.RC.Rp.Filt <- gl.filter.reproducibility(gl.CR.RC.Filt, threshold = RepFilt)
# has reduced our data set down to:
cat(nLoc(gl.CR.RC.Rp.Filt), "SNPs")



# Now that the 'quality' filters are complete, visualise the data again.  

# Histograms

# Save loc.metrics as data frame
QualFiltSummStats <- gl.CR.RC.Rp.Filt@other$loc.metrics

if(BLAST == "Y"){
  # Add a column specifying if locus has blasted to the reference (mapped/unmapped)
  QualFiltSummStats$Chrom.Mapped <- ifelse(QualFiltSummStats$Chromosome == "UnMapped", "UnMapped", "Mapped")
  
  # Reshape by metric for facet_wrap in ggplot
  QualFiltDF.hist <- QualFiltSummStats[, c("Locus", "Chrom.Mapped", "CallRate", "RepAvg", "RC_MeanTot", "maf", "meanHo", "HS", "FIS", "FST")]
}

# Reshape by metric for facet_wrap in ggplot
QualFiltDF.hist <- QualFiltSummStats[, c("Locus", "CallRate", "RepAvg", "RC_MeanTot", "maf", "meanHo", "HS", "FIS", "FST")]

QualFiltDF.hist <- pivot_longer(data = QualFiltDF.hist, cols = colnames(select_if(QualFiltDF.hist, is.double)), names_to = "MetaDat_Var")

# Rename variables so that they are in the right order in the plot
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "CallRate"] <- "a"
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "RepAvg"] <- "b"
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "RC_MeanTot"] <- "c"
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "maf"] <- "d"
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "meanHo"] <- "e"
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "HS"] <- "f"
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "FIS"] <- "g"
QualFiltDF.hist$MetaDat_Var[QualFiltDF.hist$MetaDat_Var == "FST"] <- "h"

QualFilt.hist.labs <- c("Call Rate", "Mean Repeatability", "Mean Total Read Count", "MAF", "Mean Ho", "Mean He", "FIS", "FST")
names(QualFilt.hist.labs) <- letters[1:8]

# Plot histograms- total data set
QualFilt.hist.plot <- ggplot(data = QualFiltDF.hist, aes(x = value, fill = MetaDat_Var)) +
  geom_histogram(bins = 50, colour = "black", size = 0.2) +
  facet_wrap(~ MetaDat_Var, ncol = 2, scales = "free", labeller = labeller(MetaDat_Var = QualFilt.hist.labs)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "\nValue", y = "No. Loci\n") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "None")


## Mean Ho, HS, HT plots


# Prepare individual plots for grid arrange
QualFilt.mHoHS.CR_Plot <- ggplot(QualFiltSummStats, aes(x = HS, y = meanHo, col = CallRate)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Call Rate", y= " \n ", x = " ") +
  scale_colour_viridis(option = "D", direction = -1, begin = 0.2) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

QualFilt.HS.HT.CR_Plot <- ggplot(QualFiltSummStats, aes(x = HT, y = HS, col = CallRate)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y= " \n ", x = " ") +
  scale_colour_viridis(option = "D", direction = -1, begin = 0.2) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

QualFilt.scatter.CR <- ggarrange(QualFilt.mHoHS.CR_Plot, QualFilt.HS.HT.CR_Plot, ncol = 2, labels = c("a)", "b)"), font.label = list(size = 8, face = "plain"))



QualFilt.mHoHS.RC_Plot <- ggplot(QualFiltSummStats, aes(x = HS, y = meanHo, col = RC_MeanTot)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Mean\nRead Count", y= " \n ", x = " ") +
  scale_colour_viridis(option = "A", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

QualFilt.HS.HT.RC_Plot <- ggplot(QualFiltSummStats, aes(x = HT, y = HS, col = RC_MeanTot)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y= " \n ", x = " ") +
  scale_colour_viridis(option = "A", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")


QualFilt.scatter.RC <- ggarrange(QualFilt.mHoHS.RC_Plot, QualFilt.HS.HT.RC_Plot, ncol = 2, labels = c("c)", "d)"), font.label = list(size = 8, face = "plain"))



QualFilt.mHoHS.Rep_Plot <- ggplot(QualFiltSummStats, aes(x = HS, y = meanHo, col = RepAvg)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Mean\nRepeatability", y = " \nMean Ho", x = " ") +
  scale_colour_viridis(option = "B", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.title.y = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

QualFilt.HS.HT.Rep_Plot <- ggplot(QualFiltSummStats, aes(x = HT, y = HS, col = RepAvg)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y = " \nMean He (HS)", x = " ") +
  scale_colour_viridis(option = "B", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.title.y = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

QualFilt.scatter.Rep <- ggarrange(QualFilt.mHoHS.Rep_Plot, QualFilt.HS.HT.Rep_Plot, ncol = 2, labels = c("e)", "f)"), font.label = list(size = 8, face = "plain"))




QualFilt.mHoHS.HWE_Plot <- ggplot(QualFiltSummStats, aes(x = HS, y = meanHo, col = PropPopsSig.HWE)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "HWE:\nProp. Pops Sig", y= " \n ", x = " ") +
  scale_colour_viridis(option = "C") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

QualFilt.HS.HT.HWE_Plot <- ggplot(QualFiltSummStats, aes(x = HT, y = HS, col = PropPopsSig.HWE)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y = " \n ", x = " ") +
  scale_colour_viridis(option = "C") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

QualFilt.scatter.HWE <- ggarrange(QualFilt.mHoHS.HWE_Plot, QualFilt.HS.HT.HWE_Plot, ncol = 2, labels = c("g)", "h)"), font.label = list(size = 8, face = "plain"))


QualFilt.mHoHS.FIS_Plot <- ggplot(QualFiltSummStats, aes(x = HS, y = meanHo, col = FIS)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "FIS", y= " \n ", x = "Mean He (HS)") +
  scale_colour_viridis(option = "B") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_line(size = 0.2),
        axis.text.x = element_text(colour="black",size=5,face="plain"),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.09, 0.6),
        legend.key.size = unit(0.2, "cm"))


QualFilt.HS.HT.FST_Plot <- ggplot(QualFiltSummStats, aes(x = HT, y = HS, col = FST)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "FST", y = " \n ", x = "HT") +
  scale_colour_viridis(option = "B") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_line(size = 0.2),
        axis.text.x = element_text(colour="black",size=5,face="plain"),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.09, 0.6),
        legend.key.size = unit(0.2, "cm"))

QualFilt.scatter.Fstats <- ggarrange(QualFilt.mHoHS.FIS_Plot, QualFilt.HS.HT.FST_Plot, ncol = 2, labels = c("i)", "j)"), font.label = list(size = 8, face = "plain"))


## Correlation Plot
if(BLAST == "Y"){
  QualFiltDF.cor <- QualFiltSummStats[, c("CallRate", "RC_MeanTot", "RC_MeanRef", "RC_MeanSNP", "AB", "RepAvg", EVal.col, Cnt.col, "AvgPIC", "maf", "HS", "meanHo", "HT", "FIS", "FST", "PropPopsSig.HWE")]
  QualFilt.cor <- cor(QualFiltDF.cor, use="pairwise.complete", method = "spearman")
  
  colnames(QualFilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "E Val.", "Aln. Count", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
  rownames(QualFilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "E Val.", "Aln. Count", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
} else {
  QualFiltDF.cor <- QualFiltSummStats[, c("CallRate", "RC_MeanTot", "RC_MeanRef", "RC_MeanSNP", "AB", "RepAvg", "AvgPIC", "maf", "HS", "meanHo", "HT", "FIS", "FST", "PropPopsSig.HWE")]
  QualFilt.cor <- cor(QualFiltDF.cor, use="pairwise.complete", method = "spearman")
  
  colnames(QualFilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
  rownames(QualFilt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
}

QualFilt.mHoHSHT.plot <- ggarrange(QualFilt.scatter.CR, QualFilt.scatter.RC, QualFilt.scatter.Rep, QualFilt.scatter.HWE, QualFilt.scatter.Fstats, ncol = 1)

# View
QualFilt.hist.plot

QualFilt.mHoHSHT.plot

corrplot.mixed(QualFilt.cor, upper = "circle", lower = "number", order = "FPC", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = 0.5, number.font = 1, cl.cex = 0.5, lower.col = "black", upper.col = c(viridis(n = 10, option = "D", begin = 0.1, end = 0.9), viridis(n = 10, option = "A", direction = -1, begin = 0.3)), outline = TRUE)

# Call rate is no longer having such an effect on population summary statistics (although FIS still seems to be somewhat influenced). Many of the loci with high FST and FIS values have been filtered out and the remaining loci follow the diagonal Ho/He line much more closely. 

# Output figures
jpeg(filename = paste0(outpath, Outfile_Desc, "Corr.QualFilt.jpg"), units =  "cm", width = 20, height = 20, res = 300)
corrplot.mixed(QualFilt.cor, upper = "circle", lower = "number", order = "FPC", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = 0.5, number.font = 1, cl.cex = 0.5, lower.col = "black", upper.col = c(viridis(n = 10, option = "D", begin = 0.1, end = 0.9), viridis(n = 10, option = "A", direction = -1, begin = 0.3)), outline = TRUE)
dev.off()

jpeg(filename = paste0(outpath, Outfile_Desc, "SummStat.QualFilt.jpg"), units =  "cm", width = 28, height = 20, res = 300)
ggarrange(QualFilt.hist.plot, QualFilt.mHoHSHT.plot, ncol = 2)
dev.off()


#######################
# Frequency filtering #
#######################

# Minor Allele Frequency (MAF)
# It is important not to explicitly filter on allele frequencies, in order not to bias results/shape the data set based on the expected biological patterns. One exception to this is MAF. While choosing your MAF threshold, carefully think through what this filter actually means. For example, a MAF of 0.05 will mean very different things in a data set of 10 individuals compared to a data set of 1000 (i.e. 1 copy of this allele, versus 100). In the former case, this allele could be a sequencing error (which would be unlikely in the latter given how many copies of the allele are present). However, low frequency alleles may be missed at low sample sizes. All of this needs to be weighed up when coming up with the MAF threshold (what is the acceptable error? How does this impact private alleles? How will this impact heterozygosity? What are the specific requirements/assumptions of flow on analyses?). 

# In this data set, the lowest MAF possible (1 copy of the allele) is:
1/nInd(gl.CR.RC.Rp.Filt)*2
# Therefore, if I want to filter so that an allele appears at least twice in the data set, my threshold must be greater than:
2/nInd(gl.CR.RC.Rp.Filt)*2


# Choose MAF filtering threshold:

# I chose a MAF threshold of 
MAFFilt <- 0.02
# which, after filtering:
gl.CR.RC.Rp.MAF.Filt <- gl.filter.maf(gl.CR.RC.Rp.Filt, threshold = MAFFilt)
# has reduced the data set down to:
cat(nLoc(gl.CR.RC.Rp.MAF.Filt), "SNPs") 


# Visualise how this influences private alleles
# Use dartR script to calculate number of private alleles between each population pair for filtered data set
PrivAllele_filt <- gl.report.pa(gl.CR.RC.Rp.MAF.Filt[gl.CR.RC.Rp.MAF.Filt@pop %in% PopKeep, ], plot.out = FALSE)

# Drop pop factors that have been filtered out using drop levels and only keep relevant columns (rename so can combine)
priv1_filt <- PrivAllele_filt[ , c("pop1", "priv1")]
priv2_filt <- PrivAllele_filt[ , c("pop2", "priv2")]
colnames(priv1_filt) <- c("pop", "priv")
colnames(priv2_filt) <- c("pop", "priv")

# Combine into one data set so all pop comparisons are in one column
priv_filt <- rbind(priv1_filt, priv2_filt)

if(nrow(priv_raw) & nrow(priv_filt) > 2){
  # Plot
  PrivAl.Plot.raw <- ggplot(priv_raw, aes(x = pop, y = priv, fill = pop)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.2) +
    geom_boxplot(colour = "black", size = 0.2, outlier.colour="black", outlier.size = 0.2) +
    scale_fill_viridis(discrete = T) +
    labs(x = "Population", y = "No. private alleles (pairwise comparisons among pops", title = "Raw") +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),   
          plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 6),
          axis.title.y = element_text(face = "bold", size = 6),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(face = "bold", size = 6),
          legend.position = "none")
  
  PrivAl.Plot.filt <- ggplot(priv_filt, aes(x = pop, y = priv, fill = pop)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.2) +
    geom_boxplot(colour = "black", size = 0.2, outlier.colour="black", outlier.size = 0.2) +
    scale_fill_viridis(discrete = T) +
    labs(x = "Population", y = "No. private alleles (pairwise comparisons among pops", title = paste0("MAF = ", MAFFilt)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 6),
          axis.title.y = element_text(face = "bold", size = 6),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(face = "bold", size = 6),
          legend.position = "none")
} else{
  # Plot
  PrivAl.Plot.raw <- ggplot(priv_raw, aes(x = pop, y = priv, fill = pop)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.2) +
    geom_bar(stat = "identity", colour = "black", size = 0.2) +
    scale_fill_viridis(discrete = T) +
    labs(x = "Population", y = "No. private alleles (pairwise comparisons among pops", title = "Raw") +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),   
          plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 6),
          axis.title.y = element_text(face = "bold", size = 6),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(face = "bold", size = 6),
          legend.position = "none")
  
  PrivAl.Plot.filt <- ggplot(priv_filt, aes(x = pop, y = priv, fill = pop)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.2) +
    geom_bar(stat = "identity", colour = "black", size = 0.2) +
    scale_fill_viridis(discrete = T) +
    labs(x = "Population", y = "No. private alleles (pairwise comparisons among pops", title = paste0("MAF = ", MAFFilt)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 6),
          axis.title.y = element_text(face = "bold", size = 6),
          axis.ticks = element_line(size = 0.2),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size = 5),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(face = "bold", size = 6),
          legend.position = "none")
}

# View:
ggarrange(PrivAl.Plot.raw, PrivAl.Plot.filt)

# This figure shows that the pattern of private alleles among populations doesn't really change, but the scale is drastically reduced after filtering (so filtering didn't change patterns of uniqueness/differences across populations - if it had dramatically changed,  explore further/find out why)

# Output figure
jpeg(filename = paste0(outpath, Outfile_Desc, "Priv.Al.jpg"), units =  "cm", width = 28, height = 20, res = 300)
ggarrange(PrivAl.Plot.raw, PrivAl.Plot.filt)
dev.off()


# Hardy Weinberg Equilibrium (HWE)
# Does it make sense to filter on HWE? For example, are there discrete populations? Will real results be removed that represent a biological question being tested? Do analyses assume HWE? If it's clear that the locus is out of HWE across populations, it might make sense to remove these loci for certain analyses, but in many cases, the SNPs out of HWE will be removed by quality filters anyway (i.e. if it's driven by missing data/low read count). Correct for multiple tests (Bonferroni), otherwise 5% of SNPs will be significant due to random chance, and in a large data set this means that many perfectly good SNPs will be removed for no reason.


#####################
# Linkage filtering #
#####################

# Linkage means that SNPs are close enough to each other that they are not acting independently. It's important to remove non-independent SNPs from the data set, so that they don't bias the results by strengthening false patterns. However, it's a good idea to run this filter last, as many of these non-independent SNPs will be removed by earlier filters. This means that performing the LD filter last maximises the number of quality loci that are retained.

# The first way to filter on LD is to remove multiple SNPs per sequence/locus (these are the SNPs that are known to be close enough to each other that they are not acting independently). After this, test for LD (possible even without reference genome) by determining if they are acting in a similar way (correlated). It is also important to remove sex-linked loci from the data set, as they have different patterns of inheritance to autosomal loci.

# Filter so 1 SNP per fragment
# Have to rename locus names as cloneID for this function to work
locNames(gl.CR.RC.Rp.MAF.Filt) <- gl.CR.RC.Rp.MAF.Filt@other$loc.metrics$CloneID
gl.CR.RC.Rp.MAF.lD.Filt <- gl.filter.secondaries(gl.CR.RC.Rp.MAF.Filt, method = "random")

# Filter out sex chromosomes
if(BLAST == "Y"){
  gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics$Chromosome[gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics$Chromosome == "UnMapped"] <- NA
  gl.CR.RC.Rp.MAF.lD.Filt <- gl.CR.RC.Rp.MAF.lD.Filt[ , !gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics$Chromosome %in% SexChrom]
  
  # Have to subset loc.metrics too (annoyingly it doesn't do it automatically)
  gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics <- gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics[!gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics$Chromosome %in% SexChrom, ]
}


# After removing multiple SNPs per locus and sex chromosomes, there are:
cat(nLoc(gl.CR.RC.Rp.MAF.lD.Filt), "SNPs")


# Now use SNP relate to remove linked loci
if(BLAST == "Y"){
  Chrom <- as.integer(gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics$Chromosome)
  genmat.m <- t(as.matrix(gl.CR.RC.Rp.MAF.lD.Filt)[,!is.na(Chrom)])
  SampID.m <- colnames(genmat.m)
  SNP_ID.m <- rownames(genmat.m)
  ChromPos <- gl.CR.RC.Rp.MAF.lD.Filt@other$loc.metrics$TrueChromPos[!is.na(Chrom)]
  
  snpgdsCreateGeno(paste0(result.dir, "Mapped_gl.CR.RC.Rp.MAF.lD.gds"), genmat = genmat.m, sample.id = SampID.m, snp.id = SNP_ID.m, snp.chromosome = Chrom[!is.na(Chrom)], snp.position = ChromPos, snpfirstdim = TRUE)
  
  genmat.um <- t(as.matrix(gl.CR.RC.Rp.MAF.lD.Filt)[,is.na(Chrom)])
  SampID.um <- colnames(genmat.um)
  SNP_ID.um <- rownames(genmat.um)
  
  snpgdsCreateGeno(paste0(result.dir, "UnMapped_gl.CR.RC.Rp.MAF.lD.gds"), genmat = genmat.um, sample.id = SampID.um, snp.id = SNP_ID.um, snpfirstdim = TRUE)

  gds.m <- snpgdsOpen(paste0(result.dir, "Mapped_gl.CR.RC.Rp.MAF.lD.gds"))
  gds.um <- snpgdsOpen(paste0(result.dir, "UnMapped_gl.CR.RC.Rp.MAF.lD.gds"))
  
  # LD pruning - Try different LD thresholds for sensitivity analysis
  gds_LD.m <- snpgdsLDpruning(gdsobj = gds.m, sample.id = SampID.m, snp.id = SNP_ID.m, ld.threshold = 0.5, method = "corr", autosome.only = FALSE)

  gds_LD.um <- snpgdsLDpruning(gdsobj = gds.um, sample.id = SampID.um, snp.id = SNP_ID.um, ld.threshold = 0.5, method = "corr", autosome.only = FALSE)
  
  #Subset last Dart file with these snps
  SNP.List <- c(unlist(gds_LD.m, use.names = FALSE), unlist(gds_LD.um, use.names = FALSE))

  #Filter
  gl.CR.RC.Rp.MAF.lD.SR.Filt <- gl.CR.RC.Rp.MAF.lD.Filt[ , gl.CR.RC.Rp.MAF.lD.Filt$loc.names %in% SNP.List]

  gl.CR.RC.Rp.MAF.lD.SR.Filt@other$loc.metrics <- gl.CR.RC.Rp.MAF.lD.SR.Filt@other$loc.metrics[gl.CR.RC.Rp.MAF.lD.SR.Filt@other$loc.metrics$CloneID %in% SNP.List, ]

  snpgdsClose(gds.um)
  snpgdsClose(gds.m)
  rm(gds.m, gds.um)


} else {
  genmat <- t(as.matrix(gl.CR.RC.Rp.MAF.lD.Filt))
  SampID <- colnames(genmat)
  SNP_ID <- rownames(genmat)
  
  snpgdsCreateGeno(paste0(result.dir, "gl.CR.RC.Rp.MAF.lD.gds"), genmat = genmat, sample.id = SampID, snp.id = SNP_ID, snpfirstdim = TRUE)
  gds <- snpgdsOpen(paste0(result.dir, "gl.CR.RC.Rp.MAF.lD.gds"))
  
  # LD pruning - Try different LD thresholds for sensitivity analysis
  gds_LD <- snpgdsLDpruning(gdsobj = gds, sample.id = SampID, snp.id = SNP_ID, ld.threshold = 0.5, method = "corr", autosome.only = FALSE)

  #Subset last Dart file with these snps
  SNP.List <- unlist(gds_LD, use.names = FALSE)
  
  #Filter
  gl.CR.RC.Rp.MAF.lD.SR.Filt <- gl.CR.RC.Rp.MAF.lD.Filt[ , gl.CR.RC.Rp.MAF.lD.Filt$loc.names %in% SNP.List]
  
  gl.CR.RC.Rp.MAF.lD.SR.Filt@other$loc.metrics <- gl.CR.RC.Rp.MAF.lD.SR.Filt@other$loc.metrics[ gl.CR.RC.Rp.MAF.lD.SR.Filt@other$loc.metrics$CloneID %in% SNP.List, ]
  
  snpgdsClose(gds)
  rm(gds)
}

# Filtering on linkage disequilibrium has reduced the data set down to:
cat(nLoc(gl.CR.RC.Rp.MAF.lD.SR.Filt), "SNPs")  


#######################
# Inspect final plots #
#######################

# Histograms

# Save loc.metrics as data frame
FiltSummStats <- gl.CR.RC.Rp.MAF.lD.SR.Filt@other$loc.metrics

if(BLAST == "Y"){
  # Add a column specifying if locus has blasted to the reference (mapped/unmapped)
  FiltSummStats$Chrom.Mapped <- ifelse(FiltSummStats$Chromosome == "UnMapped", "UnMapped", "Mapped")
  
  # Reshape by metric for facet_wrap in ggplot
  FiltDF.hist <- FiltSummStats[, c("Locus", "Chrom.Mapped", "CallRate", "RepAvg", "RC_MeanTot", "maf", "meanHo", "HS", "FIS", "FST")]
}

# Reshape by metric for facet_wrap in ggplot
FiltDF.hist <- FiltSummStats[, c("Locus", "CallRate", "RepAvg", "RC_MeanTot", "maf", "meanHo", "HS", "FIS", "FST")]

FiltDF.hist <- pivot_longer(data = FiltDF.hist, cols = colnames(select_if(FiltDF.hist, is.double)), names_to = "MetaDat_Var")

# Rename variables so that they are in the right order in the plot
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "CallRate"] <- "a"
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "RepAvg"] <- "b"
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "RC_MeanTot"] <- "c"
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "maf"] <- "d"
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "meanHo"] <- "e"
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "HS"] <- "f"
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "FIS"] <- "g"
FiltDF.hist$MetaDat_Var[FiltDF.hist$MetaDat_Var == "FST"] <- "h"

Filt.hist.labs <- c("Call Rate", "Mean Repeatability", "Mean Total Read Count", "MAF", "Mean Ho", "Mean He", "FIS", "FST")
names(Filt.hist.labs) <- letters[1:8]

# Plot histograms- total data set
Filt.hist.plot <- ggplot(data = FiltDF.hist, aes(x = value, fill = MetaDat_Var)) +
  geom_histogram(bins = 50, colour = "black", size = 0.2) +
  facet_wrap(~ MetaDat_Var, ncol = 2, scales = "free", labeller = labeller(MetaDat_Var = Filt.hist.labs)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "\nValue", y = "No. Loci\n") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        strip.text = element_text(face = "bold", size = 6),
        strip.background = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        axis.title.x = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 6),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "None")


## Mean Ho, HS, HT plots


# Prepare individual plots for grid arrange
Filt.mHoHS.CR_Plot <- ggplot(FiltSummStats, aes(x = HS, y = meanHo, col = CallRate)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Call Rate", y= " \n ", x = " ") +
  scale_colour_viridis(option = "D", direction = -1, begin = 0.2) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

Filt.HS.HT.CR_Plot <- ggplot(FiltSummStats, aes(x = HT, y = HS, col = CallRate)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y= " \n ", x = " ") +
  scale_colour_viridis(option = "D", direction = -1, begin = 0.2) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

Filt.scatter.CR <- ggarrange(Filt.mHoHS.CR_Plot, Filt.HS.HT.CR_Plot, ncol = 2, labels = c("a)", "b)"), font.label = list(size = 8, face = "plain"))



Filt.mHoHS.RC_Plot <- ggplot(FiltSummStats, aes(x = HS, y = meanHo, col = RC_MeanTot)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Mean\nRead Count", y= " \n ", x = " ") +
  scale_colour_viridis(option = "A", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

Filt.HS.HT.RC_Plot <- ggplot(FiltSummStats, aes(x = HT, y = HS, col = RC_MeanTot)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y= " \n ", x = " ") +
  scale_colour_viridis(option = "A", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")


Filt.scatter.RC <- ggarrange(Filt.mHoHS.RC_Plot, Filt.HS.HT.RC_Plot, ncol = 2, labels = c("c)", "d)"), font.label = list(size = 8, face = "plain"))



Filt.mHoHS.Rep_Plot <- ggplot(FiltSummStats, aes(x = HS, y = meanHo, col = RepAvg)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "Mean\nRepeatability", y = " \nMean Ho", x = " ") +
  scale_colour_viridis(option = "B", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.title.y = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

Filt.HS.HT.Rep_Plot <- ggplot(FiltSummStats, aes(x = HT, y = HS, col = RepAvg)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y = " \nMean He (HS)", x = " ") +
  scale_colour_viridis(option = "B", direction = -1) +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.title.y = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

Filt.scatter.Rep <- ggarrange(Filt.mHoHS.Rep_Plot, Filt.HS.HT.Rep_Plot, ncol = 2, labels = c("e)", "f)"), font.label = list(size = 8, face = "plain"))




Filt.mHoHS.HWE_Plot <- ggplot(FiltSummStats, aes(x = HS, y = meanHo, col = PropPopsSig.HWE)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "HWE:\nProp. Pops Sig", y= " \n ", x = " ") +
  scale_colour_viridis(option = "C") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.15, 0.6),
        legend.key.size = unit(0.2, "cm"))

Filt.HS.HT.HWE_Plot <- ggplot(FiltSummStats, aes(x = HT, y = HS, col = PropPopsSig.HWE)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(y = " \n ", x = " ") +
  scale_colour_viridis(option = "C") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=8,face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.position = "none")

Filt.scatter.HWE <- ggarrange(Filt.mHoHS.HWE_Plot, Filt.HS.HT.HWE_Plot, ncol = 2, labels = c("g)", "h)"), font.label = list(size = 8, face = "plain"))


Filt.mHoHS.FIS_Plot <- ggplot(FiltSummStats, aes(x = HS, y = meanHo, col = FIS)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "FIS", y= " \n ", x = "Mean He (HS)") +
  scale_colour_viridis(option = "B") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_line(size = 0.2),
        axis.text.x = element_text(colour="black",size=5,face="plain"),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.09, 0.6),
        legend.key.size = unit(0.2, "cm"))


Filt.HS.HT.FST_Plot <- ggplot(FiltSummStats, aes(x = HT, y = HS, col = FST)) +
  geom_jitter(size = 0.2, width = 0.003, height = 0.003, na.rm = TRUE) +
  labs(colour = "FST", y = " \n ", x = "HT") +
  scale_colour_viridis(option = "B") +
  scale_y_continuous(limits = c(0, 0.7)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  geom_abline(intercept = 0, slope = 1, size = 0.2, linetype = "dashed") + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = "transparent", size = 0.2),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),       
        axis.title.x = element_text(colour="black",size=10,face="bold"),
        axis.ticks.x = element_line(size = 0.2),
        axis.text.x = element_text(colour="black",size=5,face="plain"),
        axis.text.y = element_text(colour="black",size=5,face="plain"),
        axis.ticks.y = element_line(size = 0.2),
        legend.text = element_text(colour="black",size=4,face="plain"),
        legend.title = element_text(colour="black",size=5,face="bold"),
        legend.position = c(0.09, 0.6),
        legend.key.size = unit(0.2, "cm"))

Filt.scatter.Fstats <- ggarrange(Filt.mHoHS.FIS_Plot, Filt.HS.HT.FST_Plot, ncol = 2, labels = c("i)", "j)"), font.label = list(size = 8, face = "plain"))


## Correlation Plot
if(BLAST == "Y"){
  FiltDF.cor <- FiltSummStats[, c("CallRate", "RC_MeanTot", "RC_MeanRef", "RC_MeanSNP", "AB", "RepAvg", EVal.col, Cnt.col, "AvgPIC", "maf", "HS", "meanHo", "HT", "FIS", "FST", "PropPopsSig.HWE")]
  Filt.cor <- cor(FiltDF.cor, use="pairwise.complete", method = "spearman")
  
  colnames(Filt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "E Val.", "Aln. Count", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
  rownames(Filt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "E Val.", "Aln. Count", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
} else {
  FiltDF.cor <- FiltSummStats[, c("CallRate", "RC_MeanTot", "RC_MeanRef", "RC_MeanSNP", "AB", "RepAvg", "AvgPIC", "maf", "HS", "meanHo", "HT", "FIS", "FST", "PropPopsSig.HWE")]
  Filt.cor <- cor(FiltDF.cor, use="pairwise.complete", method = "spearman")
  
  colnames(Filt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
  rownames(Filt.cor) <- c("Call Rate", "Mean Total RC", "Mean Ref RC", "Mean SNP RC", "AB", "Mean Rep", "Mean PIC", "MAF", "HS", "Mean Ho", "HT", "FIS", "FST", "HWE")
}


# View:
Filt.hist.plot

Filt.mHoHSHT.plot <- ggarrange(Filt.scatter.CR, Filt.scatter.RC, Filt.scatter.Rep, Filt.scatter.HWE, Filt.scatter.Fstats, ncol = 1)
Filt.mHoHSHT.plot

corrplot.mixed(Filt.cor, upper = "circle", lower = "number", order = "FPC", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = 0.5, number.font = 1, cl.cex = 0.5, lower.col = "black", upper.col = c(viridis(n = 10, option = "D", begin = 0.1, end = 0.9), viridis(n = 10, option = "A", direction = -1, begin = 0.3)), outline = TRUE)


# Output figures
jpeg(filename = paste0(outpath, Outfile_Desc, "Corr.FinalFilt.jpg"), units =  "cm", width = 20, height = 20, res = 300)
corrplot.mixed(Filt.cor, upper = "circle", lower = "number", order = "FPC", tl.pos = "lt", tl.col = "black", tl.cex = 0.6, number.cex = 0.5, number.font = 1, cl.cex = 0.5, lower.col = "black", upper.col = c(viridis(n = 10, option = "D", begin = 0.1, end = 0.9), viridis(n = 10, option = "A", direction = -1, begin = 0.3)), outline = TRUE)
dev.off()

jpeg(filename = paste0(outpath, Outfile_Desc, "SummStat.FinalFilt.jpg"), units =  "cm", width = 28, height = 20, res = 300)
ggarrange(Filt.hist.plot, Filt.mHoHSHT.plot, ncol = 2)
dev.off()


# Create final PCoA and compare
gl.CR.RC.Rp.MAF.lD.Filt.pcoa <- gl.CR.RC.Rp.MAF.lD.SR.Filt
pop(gl.CR.RC.Rp.MAF.lD.Filt.pcoa) <- gl.CR.RC.Rp.MAF.lD.Filt.pcoa$other$ind.metrics[, GroupCol_PCoA]
pc <- gl.pcoa(gl.CR.RC.Rp.MAF.lD.Filt.pcoa)

gl.pcoa.plot(pc, gl.CR.RC.Rp.MAF.lD.Filt.pcoa, 
             scale = TRUE, 
             pt.size = 1,
             pop.labels = "pop", 
             save2tmp = TRUE)

# Figure out which temp file is the PCoA plot and save it
files_tempdir_plot <- list.files(tempdir())[which(str_match(list.files(tempdir()), "Plot") == "Plot")]
plot.last <- files_tempdir_plot[order(file.info(paste0(tempdir(), "/", files_tempdir_plot))$atime)][length(files_tempdir_plot[order(file.info(paste0(tempdir(), "/", files_tempdir_plot))$atime)])]
plot.no <- which(list.files(tempdir()) == plot.last)

Pcoa.filt <- gl.print.reports(print_report = plot.no)$plot

# Plot
Pcoa.filt <- Pcoa.filt +
  labs(title = "Filtered PCoA") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", 
                                    colour = "black", 
                                    fill = "transparent",
                                    size = 0.5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.title = element_text(colour="black",
                                  size=12, 
                                  face = "bold", 
                                  hjust = 0.5),
        axis.title = element_text(colour="black",
                                  size=10, 
                                  face = "bold"),
        axis.ticks = element_line(size = 0.2),
        axis.text.x = element_text(colour="black", 
                                   size=8, 
                                   face= "plain"),
        axis.text.y = element_text(colour="black", 
                                   size=8, 
                                   face= "plain"),
        legend.position = "none")

PCoA.Filt.plot <- ggarrange(ggarrange(Pcoa.raw, Pcoa.filt, ncol = 2), Samp.Map.pop, nrow = 2)

# View:
PCoA.Filt.plot

# Output figure:
jpeg(filename = paste0(outpath, Outfile_Desc, "Final.PCoA.jpg"), units =  "cm", width = 20, height = 28, res = 300)
PCoA.Filt.plot
dev.off()

# The PCoA, shows that the pattern found in the original PCoA hasn't really changed, although the PCs now explain slightly more variation (i.e. some of the noise has been removed).  

# This process will likely have to be repeated several times to get the thresholds right. Remember to filter to retain as many SNPs as possible, while also removing artifacts (driven by issues with data quality, sequencing error and non-independence) from the data. It's a trade off and as long as decisions are justifiable the end result will be a powerful, unbiased data set. In other words, don't choose thresholds based on the biological results expected, base them on quantifiable metrics that describe the quality of the data.   


################################
# Save final filtered data set #
################################

# Save genlight (so that we can load the genlight into another r session)
gl.FinalFilt <- gl.CR.RC.Rp.MAF.lD.SR.Filt
save(gl.FinalFilt, file=paste0(outpath, "gl.",Outfile_Desc, "_FinalFilt.rdata"))

# Save locus metrics
write.csv(gl.FinalFilt@other$loc.metrics, paste0(outpath, Outfile_Desc, ".FiltSumStats.csv"), row.names = F)


################
# Session Info #
################

# Output session info to file
writeLines(capture.output(sessionInfo()), paste0(outpath, Outfile_Desc, "_SessionInfo.txt"))
