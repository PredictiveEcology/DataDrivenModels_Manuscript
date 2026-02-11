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

##### pixel-level results ####
# these are growth curves with all cohorts starting at age 0 - via Biomass_speciesYield
#first get the colors 
cols <- unique(LandR::sppEquivalencies_CA[LandR %in% speciesEco$speciesCode & colorHex != "", 
                                          .(colorHex, LandR)])
colorHex <- cols$colorHex
names(colorHex) <- cols$LandR

meanMaxANPP <- ggplot(speciesEco, aes(y = maxANPP, x = speciesCode, col = source)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "species", y = "maxANPP (g/m2/yr)")

ggsave("manuscript_figures/meanMaxANPP.png", meanMaxANPP, dpi = 300, width = 7, height = 4)
#I don't know if this is a useful figure - it doesn't show much different than the results one 
googledrive::drive_upload("manuscript_figures/meanMaxANPP.png", path = gFolder, 
                          name = "meanMaxANPP.png", overwrite = TRUE)

meanMaxB <- ggplot(speciesEco, aes(y = maxB, x = speciesCode, col = source)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "species", y = "MaxB (g/m2/yr)")

ggsave("manuscript_figures/meanMaxB.png", meanMaxB, dpi = 300, width = 7, height = 4)
#I don't know if this is a useful figure - it doesn't show much different than the results one 
googledrive::drive_upload("manuscript_figures/meanMaxB.png", path = gFolder, 
                          name = "meanMaxB.png", overwrite = TRUE)



#get the original landscape
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
       aes(y = B, x = age, col = speciesCode)) + geom_line() + 
  theme_bw() + 
  scale_color_manual(name = "species", values = colorHex) + 
  ggh4x::facet_nested_wrap(vars(name, source), nrow = 3, ncol = 4)
ggsave("manuscript_figures/mostCommonEG_pixelGG.png", mostCommonEG_pixelGG, 
       dpi = 300, height = 4, width = 7)
googledrive::drive_upload("manuscript_figures/mostCommonEG_pixelGG.png", path = gFolder, 
                          name = "mostCommonEG_pixelGG.png", overwrite = TRUE)  
 
setkey(220)
six_randomPGs <- sample(randomPGs$pixelGroup, 6)
#I stupidly ran setkey instead of set.seed, but for posterity
# 3533 1898 2013  630 1736 2946
#so - focal does not always produce a more mixed stand. 
random_pixelGG <- ggplot(yieldTablesAll[pixelGroup %in% six_randomPGs], 
                               aes(y = B, x = age, col = speciesCode)) + geom_line() + 
  theme_bw() + 
  scale_color_manual(name = "species", values = colorHex) + 
  ggh4x::facet_nested_wrap(vars(name, source), nrow = 3, ncol = 4)
ggsave("manuscript_figures/random_pixelGG.png", random_pixelGG, 
       dpi = 300, height = 4, width = 7)
googledrive::drive_upload("manuscript_figures/random_pixelGG.png", path = gFolder, 
                          name = "random_pixelGG.png", overwrite = TRUE)  

#### compare the mean of all combinations #### 
meanFromB <- yieldTablesAll[, .(B = mean(B), sdB = sd(B)), 
                            .(speciesCode, age, source)]
meanBgg <- ggplot(meanFromB, aes(x = age, y = B/100, col = speciesCode, linetype = source)) + 
  geom_line(size = 1.2) + 
  scale_color_manual(name = "species", values = colorHex) + 
  labs(y = "B (Mg/ha)")
meanFromB[speciesCode == "Pinu_ban" & age < 10]


#these all have ample room for species to grow at differential rates
# limit landscape comparison to pixels below 0.25 inital B (ie with ample growing space)
# and min 2 species in a pixel (allowing competition) otherwise they grow identically between approaches
subsetPGs <- getPixSubset(initialCD = initialFocalCD, initialPG = initialFocalPG,
                          speciesEcoregion = speciesEco, prop = 0.25, minSpp = 2)

finalBwithCompetition <- pixelSummaryByThreshold(subsetPGs,
                                                 focalOutputPath = focalOutputPath,
                                                 singleOutputPath = singleOutputPath,
                                                 year = 2070) |>
  rbindlist()
#kept ecoregions in the output - not sure if useful
finalBwithCompetitionSum <- finalBwithCompetition[, .(B = sum(B),meanB = mean(B)), 
                                               .(speciesCode, source)]

#area is wrong here but its wrong for both scenarios equally :) 
# B is per g/m2 within a pixel, summed by pixel, but we know a pixel is 6.25 ha
BafterTime <- ggplot(finalBwithCompetitionSum, aes(y = B/100 * 6.25, x = speciesCode, fill = source)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "species", y = "Biomass (Mg)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("manuscript_figures/Results_100yrWithCompetition.png", BafterTime, dpi = 300, width = 7, height = 4)
googledrive::drive_upload("manuscript_figures/Results_100yrWithCompetition.png", path = gFolder, 
                          name = "Results_100yrWithCompetition.png", overwrite = TRUE)

#I guess as a percentage change?
finalBWide <- data.table::dcast(finalBwithCompetition, speciesCode + ecoregionGroup + N ~ source, value.var ="B")
finalBWide[, pctChange := c(focal - single)/single * 100]
finalBWide[, rawChange_sum := c(focal - single)]
finalBWide[, rawChange_MgPerHa := rawChange_sum/N/100] #

pctBAfterTime <- ggplot(finalBWide, aes(x = pctChange, y = speciesCode)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "change in B by ecoregionGroup under Focal compared to Single (%)",
       y = "species") + 
  # scale_fill_manual(values = colorHex) + 
  facet_wrap(~ecoregionGroup, ncol = 6)
ggsave("manuscript_figures/Results_100yrPctByEcoregion.png", pctBAfter100, 
       dpi = 300, height = 4, width = 7)
googledrive::drive_upload("manuscript_figures/Results_100yrPctByEcoregion.png", path = gFolder, 
                          name = "Results_100yrPctByEcoregion.png", overwrite = TRUE)  

rawBAfter100  <- ggplot(finalBWide, aes(x = rawChange, y = speciesCode)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "change in B under Focal compared to Single (Mg)",
       y = "species") + 
  # scale_fill_manual(values = colorHex) + 
  facet_wrap(~ecoregionGroup, ncol = 6, scales = "free_x")
  
averagedB  <- ggplot(finalBWide[grep("_081", ecoregionGroup, invert = TRUE)], aes(x = rawChange_MgPerHa, y = speciesCode)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "change in B under Focal compared to Single (Mg)",
       y = "species") + 
  # scale_fill_manual(values = colorHex) + 
  facet_wrap(~ecoregionGroup)






####total relative biomass, mean aNPP by sim ####
#summarize cohort data by pixelGroup and species, join with pixelGroupMap, summarize the landscape

# focalLandscape <- landscapeSummary(focalOutputPath) |> rbindlist()
# focalLandscape[, source := "focal"]
# 
# singleLandscape <- landscapeSummary(singleOutputPath) |> rbindlist()
# singleLandscape[, source := "single"]
# landscapeDF <- rbind(singleLandscape, focalLandscape)
# 
# speciesBiomassPlot <- function(df, cols, y, species, ylab = "y",
#                                plotTitle = NULL, plotSubtitle = NULL) {
#   gg <- ggplot(data = df, aes(x = .data[["year"]], y = .data[[y]], fill = .data[[species]], group = .data[[species]])) +
#     geom_area(position = "stack") +
#     scale_fill_manual(values = cols) +
#     # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
#     scale_y_continuous() +
#     labs(x = "Year", y = ylab, fill = "Species", title = plotTitle, subtitle = plotSubtitle) +
#     theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
#     theme_bw(base_size = 16) + 
#     facet_wrap(~source, nrow = 1)
#   return(gg)
# }
# cols <- unique(LandR::sppEquivalencies_CA[LandR %in% landscapeDF$speciesCode & colorHex != "", 
#                                           .(colorHex, LandR)])
# plotCols <- cols$colorHex
# names(plotCols) <- cols$LandR
# speciesBiomassPlot(landscapeDF, cols = plotCols, y = "meanB", ylab = "mean landscape biomass (g/m2)", 
#                    species = "speciesCode")
# propB_all <- speciesBiomassPlot(landscapeDF, cols = plotCols, y = "propB", ylab = "proportion of landscape AGB", 
#                    species = "speciesCode")
# propB_100yr <- ggplot(landscapeDF[year == 2120], aes(y = propB, x = speciesCode, fill = source)) + 
#   geom_bar(position = 'dodge', stat = 'identity') + 
#   labs(y = "proportional landscape B after 100 yrs")
# ggsave(filename = "manuscript_figures/propB_100yr.png", 
#        propB_100yr, dpi = 300,height = 8, width = 8)
# 
# landscapeB_100yr <- ggplot(landscapeDF[year == 2120], aes(y = meanB/100, x = speciesCode, fill = source)) + 
#   geom_bar(position = 'dodge', stat = 'identity') + 
#   labs(y = "mean landscape B after 100 yrs (Mg/ha)")
# ggsave(filename = "manuscript_figures/landscapeB_100yr.png", 
#        landscapeB_100yr, dpi = 300, height = 8, width = 8)
# 


#compare some pixelGroups
# age is not great because without a spin up, and with poor initial age data, 
# many pixels have too much biomass for their age. Therefore, to capture differences, 
# examine pixels with available growing space 
#


#TODO: make this a panel of multiple pixelIndices - 
# report the ecoregion too
# ? join speciesCode with colour
#


#what did inflationFactor look like in the theoretical species?
#get the factorial used in the simulation:
#(they are identical between reps)


# not sure if we need this
# ggsave(PSPages, filename = "manuscript_figures/MontaneCordillera_PSP_standAge_histogram.png", height = 4, width = 8, dpi = 300)
# googledrive::drive_upload("manuscript_figures/MontaneCordillera_PSP_standAge_histogram.png", path = gFolder,
#                           name = "MontaneCordillera_PSP_standAge_histogram.png")

# age20CDs_sample <- age20CDs[pixelIndex %in% sample(x = pixelIndex, size = 4, replace = FALSE)]
# randomPixel <- ggplot(age20CDs_sample, aes(y = B, x = year, col = speciesCode, linetype = source)) + 
#   geom_line() +
#   theme_bw() + 
#   scale_fill_manual(colorHex) + 
#   labs(title = "random pixels with stand-age 20 in 2020", y = "B (g/m2)") + 
#   facet_wrap(~pixelIndex, scales = "free_y")
# randomPixel
# ggsave("manuscript_figures/randomPixel5135.png", randomPixel, width = 8, height = 4, dpi = 300)
