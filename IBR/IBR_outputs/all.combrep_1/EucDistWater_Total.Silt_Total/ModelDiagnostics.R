library(ResistanceGA)
library(lme4)
library(DHARMa)
library(tidymodels)
library(ggplot2)
library(LMERConvenienceFunctions)
library(heavy)

setwd("/Users/robynshaw/Google Drive/Work/Pilbara_Small_Mammals/Manuscripts/Quoll_Paper/Full_Analysis/IBR/all.combrep_1_VRM_FinalForPaper/EucDistWater_Total.Silt_Total/")

# Transform surfaces
options(scipen = 999) # remove scientific notation so that axis is in 1000’s

# First plot the transformed surface for distance2water
# Get parameter values from Multisurface_Optim_Summary.txt in top performing model output

jpeg("Dist2Water_trans.jpg", width = 18, height = 18, units = "cm", res = 300)
Plot.trans(transformation = 1.50704, PARM = c(6.715469, 1621.32), Resistance = "../../../../Rasters_5km.Aligned_BufferCrop/EucDistWater_TotalBuffCrop.asc", Name = "Distance to water")
dev.off()

jpeg("Silt_trans.jpg", width = 18, height = 18, units = "cm", res = 300)
Plot.trans(transformation = 1.640336, PARM = c(0.592705, 2352.323), Resistance = "../../../../Rasters_5km.Aligned_BufferCrop/EucDistWater_TotalBuffCrop.asc", Name = "Distance to water")
dev.off()




###############################################
# REQUIRED HELPER FUNCTIONS FROM RESISTANCEGA #
###############################################

# The overdispersion function was found here:
# GLMM FAQ, Bolker et al. (2020):
# http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor.

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


# Get lower trioangle from dist matrix
lower <- function (matrix) {
  if (is.vector(matrix) == TRUE || dim(matrix)[1] != dim(matrix)[2]) {
    warning("Must provide square distance matrix with no column or row names")
  }
  lm <- matrix[lower.tri(matrix)]
  return(lm)
}


# Prepare random effects
ZZ.mat <- function (ID, drop = NULL) {
  ID.num <- ID
  if (any("corr_" %in% names(ID))) {
    if (any("pop1.pop" %in% names(ID))) {
      ID.num$pop1 <- as.numeric(ID.num$pop1.ind)
      ID.num$pop2 <- as.numeric(ID.num$pop2.ind)
      Zl.corr <- lapply(c("pop1.ind", "pop2.ind"), function(nm) Matrix::fac2sparse(ID[[nm]], 
                                                                                   "d", drop = FALSE))
    }
    else {
      ID.num$pop1 <- as.numeric(ID.num$pop1)
      ID.num$pop2 <- as.numeric(ID.num$pop2)
      Zl.corr <- lapply(c("pop1", "pop2"), function(nm) Matrix::fac2sparse(ID[[nm]], 
                                                                           "d", drop = FALSE))
    }
    for (i in 1:dim(Zl.corr[[1]])[1]) {
      Zl.corr[[1]][i, ] <- (ID.num$pop1 == i & ID.num$corr_ == 
                              1)
    }
    for (i in 1:dim(Zl.corr[[2]])[1]) {
      Zl.corr[[2]][i, ] <- (ID.num$pop2 == i & ID.num$corr_ == 
                              1)
    }
    ZZ.corr <- Reduce("+", Zl.corr[-1], Zl.corr[[1]])
    ZZ.corr <- Matrix::drop0(ZZ.corr)
  }
  if (any("pop1.pop" %in% names(ID))) {
    Zl <- lapply(c("pop1.ind", "pop2.ind"), function(nm) Matrix::fac2sparse(ID[[nm]], 
                                                                            "d", drop = FALSE))
    ZZ.ind <- Reduce("+", Zl[-1], Zl[[1]])
    Zl <- lapply(c("pop1.pop", "pop2.pop"), function(nm) Matrix::fac2sparse(ID[[nm]], 
                                                                            "d", drop = FALSE))
    ZZ.pop <- Reduce("+", Zl[-1], Zl[[1]])
    if (any("corr_" %in% names(ID))) {
      ZZ <- rbind(ZZ.ind, ZZ.pop, ZZ.corr)
    }
    else {
      ZZ <- rbind(ZZ.ind, ZZ.pop)
    }
  }
  else {
    Zl <- lapply(c("pop1", "pop2"), function(nm) Matrix::fac2sparse(ID[[nm]], 
                                                                    "d", drop = FALSE))
    ZZ <- ZZ.pop <- Reduce("+", Zl[-1], Zl[[1]])
    if (any("corr_" %in% names(ID))) {
      ZZ <- rbind(ZZ.pop, ZZ.corr)
    }
  }
  if (!is.null(drop)) {
    ZZ <- ZZ[, drop == 1]
  }
  return(ZZ)
}

############################################
# RUN MLPE AND PERFORM TESTS FOR MODEL FIT #
############################################

# Read in resistance distance (from best model output) and
# genetic distance (original input to resistanceGA)

resistance <- lower(read.csv("EucDistWater_Total.Silt_Total_commuteDistance_distMat.csv", header = FALSE))

response <- read.csv("../../Coord_Genetic_Data_Dh/GeneticMetrics_Dh.csv")$EucDist

# Generate population random effect
Samps <- as.matrix(read.csv("../../Coord_Genetic_Data_Dh/Coords_Dh.csv", stringsAsFactors = FALSE))
ID <- To.From.ID(nrow(Samps))
ZZ <- ZZ.mat(ID)
REML <- FALSE

mm <- resistance
m <- 0.5 * (sqrt((8 * length(mm)) + 1) + 1)
mm <- mm[which(mm != -1)]
cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
dat <- data.frame(ID, resistance = cs.matrix, response = response)
colnames(dat) <- c("pop1", "pop2", "resistance", "response")

# run model
MOD <- MLPE.lmm(resistance = resistance, 
                pairwise.genetic = response, 
                REML = F, 
                ID = ID, 
                ZZ = ZZ)

summary(MOD)
car::qqPlot(residuals(MOD))
overdisp_fun(MOD)
testDispersion(MOD)
qqnorm(residuals(MOD), main="MLPE Residuals")
qqline(residuals(MOD), col="blue")
confint.merMod(MOD, method = "boot")


MOD.rmOL <- romr.fnc(MOD, data = dat, trim = 2.5)
dat.new <- MOD.rmOL$data

unique(dat.new[, 1:2])
ID <- unique(dat.new[, 1:2])
ZZ <- ZZ.mat(ID)

MOD.rmOL.new <- MLPE.lmm(resistance = dat.new$resistance, 
                         pairwise.genetic = dat.new$response, 
                         REML = F, 
                         ID = ID, 
                         ZZ = ZZ)

summary(MOD.rmOL.new)
car::qqPlot(residuals(MOD.rmOL.new))
overdisp_fun(MOD.rmOL.new)
testDispersion(MOD.rmOL.new)
qqnorm(residuals(MOD.rmOL.new), main="MLPE Residuals")
qqline(residuals(MOD.rmOL.new), col="blue")
confint.merMod(MOD.rmOL.new)

# Test residuals for normality
set.seed(1234)
shapiro.test(residuals(MOD))# Strongly Rejected
ks.test(residuals(MOD), pnorm, 0,1)# Strongly Rejected
ks.test(residuals(MOD), pt, 2)# Strongly Rejected

shapiro.test(residuals(MOD.rmOL.new))# Strongly Rejected
ks.test(residuals(MOD.rmOL.new), pnorm, 0,1)# Strongly Rejected
ks.test(residuals(MOD.rmOL.new), pt, 2)# Strongly Rejected




library(heavy)
MOD.HT <- heavyLme(fixed = response ~ resistance,
              random = NULL,
              groups = ~ pop1,
              data = dat)
summary(MOD.HT)

car::qqPlot(MOD.HT$Resid$marginal)
car::qqPlot(MOD.HT$Resid$conditional)





library(ggpubr)
ggqqplot(MOD.HT$Resid$conditional)


# Check if model is overdispersed:
as.data.frame(overdisp_fun(MOD)) # Nope
as.data.frame(overdisp_fun(MOD.rmOL.new))

# A note about R2
# The marginal R2 values for all models relating to R. tunneyi capture rate are low. This reflects the fact that we had high trapping effort and low captures. This means there are many (true) zeros in our dataset, because R. tunneyi is patchily distributed across the landscape. A particular area may have suitable habitat with no captures, simply due to low densities. However, these vegetation types are significant predictors of capture rate when R. tunneyi are present. Furthermore, we are modeling a highly variable system and there is an element of individual R. tunneyi behaviour that is impossible to model (i.e. why do animals choose to move through the specific area where a trap was placed).

# In other words, our primary goal is to understand the relationships between the vegetation types (and area burnt later on) and R. tunneyi capture patterns in our model. For this objective, R2 is not particularly relevant. The regression coefficients define the relationship between each independent variable and the dependent variable. A small R2 doesn’t nullify or change the interpretation of the coefficient for an independent variable that is statistically significant.

# Instead, we compare the delta BIC of the top ranked model/s to the null model (intercept and [where relevant] the spatial autocovariate only), as a more informative measure of model fit. In this case, a smaller delta BIC indicates that the model is improved by the inclusion of fixed effects despite the low R2.


# Check goodness of fit
set.seed(123)
SimTop <- simulateResiduals(fittedModel = MOD, n = 1000)
plot(SimTop)

SimTop <- simulateResiduals(fittedModel = MOD.rmOL.new, n = 1000)
plot(SimTop)




# A note on significance of dispersion test from the authors of DHARMa
# Significance in hypothesis tests depends on the strength of the signal and the number of data points. Hence, the p-value alone is not a good indicator of the extent to which your residuals deviate from assumptions. Specifically, if you have a lot of data points, residual diagnostics will inevitably become significant, because having a perfectly fitting model is very unlikely. That, however, doesn’t necessarily mean that you need to change your model; it is at your discretion to decide whether this deviation is worth worrying about. e.g. A significant test with a dispersion parameter of 1.01 probably reflects a large data set, while a significant test with a dispersion value of 5 is clearly a reason to move to a model that accounts for overdispersion.

# In my case, the dispersion ratio is 1.16 so the significant test probably reflects the large number of data points. For interest, I also fit a model with an observation level random effect, however it was not a valid model.

