#eventually this could be a module run with Spades Experiment 
#These are Results figures created from simulation output
library(Require)
Require(data.table)
Require(ggplot2)
Require(arrow)
Require(terra)
Require(qs2)
Require(ggh4x)
source("R/results_functions.R")


gFolder <- googledrive::as_id("https://drive.google.com/drive/folders/1-jrQpHIyUPsveTH5gw1fsdyd3txj8TvP?usp=drive_link")
singleOutputPath <- "outputs/single"
focalOutputPath <- "outputs/focal"
ecoregionMap <- rast(file.path(focalOutputPath, "ecoregionMap_year2120.tif"))
overwriteFigures <- FALSE

#### set up #####
# fix names of ecoregionGroups
ecoTable <- terra::cats(ecoregionMap)
ecoVals <- data.table(ID = as.vector(ecoregionMap))
ecoVals <- ecoVals[!is.na(ID), .N, .(ID)]
ecoTable <- ecoVals[ecoTable, on = c("ID")]
ecoTable <- ecoTable[order(N, decreasing = TRUE)]
ecoD_names <- fread("ecodistricts.csv")
ecoD_names[, Ecodistrict := as.character(Ecodistrict)]

ecoTable <- ecoD_names[ecoTable, on = c("Ecodistrict" = "ecoregionName")]
lccs <- data.table(lcc = c("081", "210", "220", "230"), 
                   lccName = c("wetland", "conif.", "decid.", "mixed"))
ecoTable <- ecoTable[lccs, on = c("landcover" = "lcc")]
ecoTable[, name := paste0(Name, "/ ", lccName)]
rm(ecoVals, lccs)

#### figures #####


#no stochasticity in species traits, so use rep 1 for simplicity
singleSpp <- readRDS(file.path(singleOutputPath, "species_year2120.rds"))
focalSpp <- readRDS(file.path(focalOutputPath, "species_year2120.rds"))

singleSpp[, source := "single"]
focalSpp[, source := "focal"]
bothSpp <- rbind(singleSpp, focalSpp) 
speciesRep <- melt(bothSpp, 
                   id.vars = c("source", "growthCurveSource", "species"), 
                   measure.vars = c("longevity", "growthcurve", "mANPPproportion", 
                                    "mortalityshape", "inflationFactor"), 
                   variable.name = "growth param.",
                   value.name = "value"
)
#for reporting in-text
bothSpp[, .(varGC = var(growthcurve), meanGC = mean(growthcurve), 
            varANPP = var(mANPPproportion), meanANP = mean(mANPPproportion),
            varMS = var(mortalityshape), meanMS = mean(mortalityshape),
            varIF = var(inflationFactor), meanIF = mean(inflationFactor)), 
        .(source)]
#longevity has no perceptible difference 

speciesEco_single <- readRDS(file.path(singleOutputPath,  "speciesEcoregion_year2120.rds"))
speciesEco_focal <- readRDS(file.path(focalOutputPath,  "speciesEcoregion_year2120.rds"))
speciesEco_focal[, source := "focal"]
speciesEco_single[, source := "single"]
speciesEco <- rbind(speciesEco_focal, speciesEco_single)
speciesEco[, N := .N, .(speciesCode)]

rm(focalSpp, singleSpp, speciesEco_single,  speciesEco_focal)

#report this
setkey(bothSpp, species)
fwrite(bothSpp[, .(species, growthcurve, mANPPproportion, mortalityshape, inflationFactor, source)], "outputs/speciesStats.csv")
##### pixel-level results ####
# these are growth curves with all cohorts starting at age 0 - via Biomass_speciesYield
#first get the colors 
cols <- LandR::sppEquivalencies_CA[LandR %in% speciesEco$speciesCode & colorHex != "", .(colorHex, LandR)] |>
  unique()
colorHex <- cols$colorHex
names(colorHex) <- cols$LandR

meanMaxANPP <- ggplot(speciesEco, aes(y = maxANPP, x = speciesCode, col = source)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "species", y = "maxANPP (g/m2/yr)")

if (overwriteFigures) {
  ggsave("manuscript_figures/meanMaxANPP.png", meanMaxANPP, dpi = 300, width = 7, height = 4)
  googledrive::drive_upload("manuscript_figures/meanMaxANPP.png", path = gFolder,
                            name = "meanMaxANPP.png", overwrite = TRUE)
}

meanMaxB <- ggplot(speciesEco, aes(y = maxB, x = speciesCode, col = source)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "species", y = "MaxB (g/m2/yr)")

if (overwriteFigures) {
  ggsave("manuscript_figures/meanMaxB.png", meanMaxB, dpi = 300, width = 7, height = 4)
  #I don't know if this is a useful figure - it doesn't show much different than the results one 
  googledrive::drive_upload("manuscript_figures/meanMaxB.png", path = gFolder, 
                            name = "meanMaxB.png", overwrite = TRUE)
}

# get the original landscape
initialFocalCD <- readRDS(file.path(focalOutputPath, "cohortData_year2020.rds"))
initialFocalPG <- terra::rast(file.path(focalOutputPath, "pixelGroupMap_year2020.tif"))

#1. get Yield tables
yieldTablesAll <- lapply(list(focalOutputPath, singleOutputPath), 
                         FUN = getBiomassYieldCD)
yieldTablesAll[[1]][, source := "focal"]
yieldTablesAll[[2]][, source := "single"]
yieldTablesAll <- rbindlist(yieldTablesAll)

##pixels under age 10 have incorrect age due to a quirk in Biomass_core
##update it here
yieldTablesAll[age < 10, age := rank(B), .(speciesCode, pixelGroup, source)]
yieldTablesAll <- yieldTablesAll[ecoTable[, .(ecoregionGroup, name)], on = c("ecoregionGroup")]

#2. get subset of pixelGroups for figures
# in the event we want random ecoregionGroups for pixel-level results... 
# it is unclear if we want random or not.. alternatively, the most common pairing?
# but the most common pairings do not contain all species (e.g. no Pice_mariana)

### approach 1 - random samples ####
#randomly sample PGs, stratified by ecoregionGroup
setkey(102)
randomPGs <- yieldTablesAll[, .(pixelGroup = sample(pixelGroup, size = 1)), 
                            .(ecoregionGroup)]

#### approach 2 -  most common combinations by species ####
subCDLong <- LandR::addPixels2CohortData(initialFocalCD, initialFocalPG)
subCDLong <- subCDLong[, .(pixelIndex, ecoregionGroup, speciesCode)]
subCDLong[, nSpp := length(unique(speciesCode)), .(pixelIndex)]
#drop the single-species pixels, they are uninteresting
subCDLong <- subCDLong[nSpp > 1, .(combo = paste(unique(speciesCode), collapse = ", ")), 
                       .(pixelIndex, ecoregionGroup)]
subCDLong <- subCDLong[, .N, .(ecoregionGroup, combo)]
subCDLong[, mostCommon := max(N), .(ecoregionGroup)]
subCDLong <- subCDLong[N == mostCommon, .(ecoregionGroup, combo)]
#match this with the biomass_yieldTable pixelGorups
bytCD <- qs2::qs_read(file.path(focalOutputPath, 'cohortDataYield/cohortData_year000.qs2'))
bytCD <- bytCD[, .(combo = paste(unique(speciesCode), collapse = ", ")), .(pixelGroup, ecoregionGroup)]
mostCommonPGsubset <- bytCD[subCDLong, on = c("ecoregionGroup", "combo")]
# rm(bytCD, samplePix)


# this table contains the most common pairing of species by ERG, their associated PG,
# and the abundance of each ERG (N)
commonPGs <- ecoTable[mostCommonPGsubset, on = c("ecoregionGroup")]
commonPGs <- commonPGs[order(N, decreasing = TRUE)]


#use consistent colouring so species are the same between different plots
# this is the 6 most common ERG
mostCommonEG_pixelGG <- ggplot(yieldTablesAll[pixelGroup %in% commonPGs[1:6,]$pixelGroup], 
       aes(y = B/100, x = age, col = speciesCode)) + geom_line() + 
  theme_bw() + 
  scale_color_manual(name = "species", values = colorHex) + 
  labs(y = "Biomass (Mg/ha)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggh4x::facet_nested_wrap(vars(name, source), nrow = 3, ncol = 4)

if (overwriteFigures) {
  ggsave("manuscript_figures/mostCommonEG_pixelGG.png", mostCommonEG_pixelGG, 
         dpi = 300, height = 4, width = 7)
  googledrive::drive_upload("manuscript_figures/mostCommonEG_pixelGG.png", path = gFolder, 
                            name = "mostCommonEG_pixelGG.png", overwrite = TRUE)  
}
setkey(220)
six_randomPGs <- sample(randomPGs$pixelGroup, 6)
# [1] 3593  801 3431  681 2951 3192
#so - focal does not always produce a more mixed stand. 
random_pixelGG <- ggplot(yieldTablesAll[pixelGroup %in% six_randomPGs], 
                               aes(y = B/100, x = age, col = speciesCode)) + geom_line() + 
  theme_bw() + 
  scale_color_manual(name = "species", values = colorHex) + 
  labs(y = "Biomass (Mg/ha)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  ggh4x::facet_nested_wrap(vars(name, source), nrow = 3, ncol = 4)

if (overwriteFigures) {
  ggsave("manuscript_figures/random_pixelGG.png", random_pixelGG, 
         dpi = 300, height = 4, width = 7)
  googledrive::drive_upload("manuscript_figures/random_pixelGG.png", path = gFolder, 
                            name = "random_pixelGG.png", overwrite = TRUE)  
}

#### compare the mean of all combinations #### 
#unclear if this figure will be used 
# meanFromB <- yieldTablesAll[, .(B = mean(B), sdB = sd(B)), 
#                             .(speciesCode, age, source)]
# meanBgg <- ggplot(meanFromB, aes(x = age, y = B/100, col = speciesCode, linetype = source)) + 
#   geom_line(size = 1.2) + 
#   scale_color_manual(name = "species", values = colorHex) + 
#   labs(y = "B (Mg/ha)")


#these all have ample room for species to grow at differential rates
# limit landscape comparison to pixels below 0.25 inital B (ie with ample growing space)
# and min 2 species in a pixel (allowing competition) otherwise they grow identically between approaches
#53K are < 25%B, ~40K of those 

cdLong <-  LandR::addPixels2CohortData(initialFocalCD, pixelGroupMap = initialFocalPG)
nPixTotal <- length(unique(cdLong$pixelIndex))
maxB_eco <- speciesEco[source == "focal", .(ecoMaxB = max(maxB)), .(ecoregionGroup)]
cdLong[, sumB := sum(B), .(pixelIndex)]
cdLong <- cdLong[maxB_eco, on = c("ecoregionGroup")]
hist(cdLong$ecoMaxB - cdLong$sumB, xlab = "growing space (g/m2)")

cdLong[, current_sumB_pct := c(sumB/ecoMaxB) * 100]
cdLong[, nSpp := length(speciesCode), .(pixelIndex)]

#subset pixels from the landscape based on starting biomass and competition
subsetPGs <- getPixSubset(initialCD = initialFocalCD, initialPG = initialFocalPG,
                          speciesEcoregion = speciesEco[source == "focal"], prop = 0.25, minSpp = 2)

LandscapeB_overTime <- lapply(seq(2020, 2120, by = 10), 
                           pixelSummaryByThreshold, 
                           PixUnderThresh = subsetPGs,
                           focalOutputPath = focalOutputPath,
                           singleOutputPath = singleOutputPath
) |>
  rbindlist()

LandscapeB_overTime <- LandscapeB_overTime[, .(B = sum(B),meanB = mean(B)), 
                                               .(speciesCode, source, Year)]


# B is per g/m2 within a pixel, summed by pixel, and a pixel is 6.25 ha (by default)
BafterTime <- ggplot(LandscapeB_overTime, aes(y = B/100 * 6.25, x = Year, col = speciesCode, linetype = source)) + 
  geom_line(linewidth = 1.2) +
  labs(x = "Year", y = "Landscape biomass (Mg)") + 
  scale_color_manual(name = "species", values = colorHex) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

if (overwriteFigures) {
  ggsave("manuscript_figures/emptyPixels_after100yrs.png", BafterTime, dpi = 300, width = 7, height = 4)
  googledrive::drive_upload("manuscript_figures/emptyPixels_after100yrs.png", path = gFolder, 
                            name = "emptyPixels_after100yrs.png", overwrite = TRUE)
}


#find the percentage of initial PGs that yield mixed stands under each approach
mixedProp <- 0.8

yieldTablesAll[, sumB := sum(B), .(pixelGroup, source, age)]
yieldTablesAll[, nSpp := .N, .(pixelGroup, source, age)]
yieldTablesAll[, propB := B/sumB]
yieldTablesAll[, maxSpeciesProp := max(propB), .(pixelGroup, age, source)]
#take PGs with > 1 species, limit to the maximum
propB_summary <- yieldTablesAll[nSpp >1, 
                                .(isMixed = maxSpeciesProp < mixedProp,
                                  isLeading = propB >= mixedProp),
                                .(speciesCode, pixelGroup, source)]
propB_summary <- propB_summary[, .(isMixed = sum(isMixed), isNotMixed = sum(!isMixed),
                                   isLeading = sum(isLeading), isNotLeading = sum(!isLeading)), 
                                   .(source, speciesCode)]
propB_summary[, Leading := isLeading/(isLeading + isNotLeading) * 100]
propB_summary[, Mixed := isMixed/(isMixed + isNotMixed) * 100]
propB_summary[, Dominated := 100 - Leading - Mixed]


mixedStats <- melt.data.table(propB_summary, 
                                 id.vars = c("speciesCode", "source"), 
                                 measure.vars = c("Leading", "Dominated", "Mixed"), 
                              variable.name = "ValueInStand", value.name = "percent")

# label_map <- c(pctLeading = "species is leading", 
#                pctMixed = "species is in a mixed stand")
mixedGG <- ggplot(mixedStats, aes(y = percent, x = speciesCode, fill = ValueInStand)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  facet_wrap(~source, nrow = 1) + 
  labs(x = "species", y = "% across all pixelGroups") + 
  scale_fill_discrete(name = "Stand class")
  # scale_fill_discrete(name = "", labels = c(pctLeading = "Species is leading",
  #                                           pctMixed = "Stand is mixed"))
mixedGG



# how prevalent each species was in initial pg?
# initialAbundance <- unique(yieldTablesAll[, .(speciesCode, pixelGroup)])
# initialAbundance <- initialAbundance[, .N, .(speciesCode)]
# initialAbundance
# speciesCode     N
# <fctr> <int>
# 1:    Pice_gla   819
# 2:    Pice_mar   770
# 3:    Popu_bal   585
# 4:    Popu_tre  1264
# 5:    Pinu_con   607
# 6:    Pice_eng   386

# #I guess as a percentage change?
# finalBWide <- data.table::dcast(finalBwithCompetition, speciesCode + ecoregionGroup + N ~ source, value.var ="B")
# finalBWide[, pctChange := c(focal - single)/single * 100]
# finalBWide[, rawChange_sum := c(focal - single)]
# finalBWide[, rawChange_MgPerHa := rawChange_sum/N/100] #
# 
# pctBAfterTime <- ggplot(finalBWide, aes(x = pctChange, y = speciesCode)) + 
#   geom_bar(position = "dodge", stat = "identity") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   labs(x = "change in B by ecoregionGroup under Focal compared to Single (%)",
#        y = "species") + 
#   # scale_fill_manual(values = colorHex) + 
#   facet_wrap(~ecoregionGroup, ncol = 6)
# ggsave("manuscript_figures/Results_100yrPctByEcoregion.png", pctBAfter100, 
#        dpi = 300, height = 4, width = 7)
# googledrive::drive_upload("manuscript_figures/Results_100yrPctByEcoregion.png", path = gFolder, 
#                           name = "Results_100yrPctByEcoregion.png", overwrite = TRUE)  
# 
# rawBAfter100  <- ggplot(finalBWide, aes(x = rawChange, y = speciesCode)) + 
#   geom_bar(position = "dodge", stat = "identity") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   labs(x = "change in B under Focal compared to Single (Mg)",
#        y = "species") + 
#   # scale_fill_manual(values = colorHex) + 
#   facet_wrap(~ecoregionGroup, ncol = 6, scales = "free_x")
#   
# averagedB  <- ggplot(finalBWide[grep("_081", ecoregionGroup, invert = TRUE)], aes(x = rawChange_MgPerHa, y = speciesCode)) + 
#   geom_bar(position = "dodge", stat = "identity") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   labs(x = "change in B under Focal compared to Single (Mg)",
#        y = "species") + 
#   # scale_fill_manual(values = colorHex) + 
#   facet_wrap(~ecoregionGroup)
# 

