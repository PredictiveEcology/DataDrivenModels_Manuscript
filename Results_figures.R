#These are Results figures created from simulation output

library(Require)
Require(data.table)
Require(ggplot2)
Require(arrow)
source("R/results_functions.R")

gFolder <- googledrive::as_id("https://drive.google.com/drive/folders/1-jrQpHIyUPsveTH5gw1fsdyd3txj8TvP?usp=drive_link")


singleOutputPath <- "outputs/singleFitting_MC"
focalOutputPath <- "outputs/focalFitting_MC"
#no stochasticity in species traits, so use rep 1 for simplicity
singleRep1 <- readRDS(file.path(singleOutputPath, "rep1/species_year2120.rds"))
focalRep1 <- readRDS(file.path(focalOutputPath, "rep1/species_year2120.rds"))

singleRep1[, source := "single"]
focalRep1[, source := "focal"]
speciesRep1 <- rbind(singleRep1, focalRep1) 
speciesRep1[, maxB_multiplier := inflationFactor - 1] #TODO: maybe don't do this - or rethink how best to show differences
speciesRep <- melt(speciesRep1, 
                   id.vars = c("source", "growthCurveSource", "species"), 
                   measure.vars = c("longevity", "growthcurve", "mANPPproportion", 
                                    "mortalityshape", "maxB_multiplier"), 
                   variable.name = "growth param.",
                   value.name = "value"
)

#longevity has no perceptible difference 
sppPlot <- ggplot(speciesRep[`growth param.` != "longevity"], aes(y = value, x = species, fill = source)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_wrap(~`growth param.`, 
            scales = "free", ncol = 1) + theme_bw()
ggsave("manuscript_figures/sppPlot.png", sppPlot, dpi = 300, width = 8, height = 8)
googledrive::drive_upload("manuscript_figures/sppPlot.png", path = gFolder, name = "sppPlot.png")



#these will be identical between reps
speciesEco_single <- readRDS(file.path(singleOutputPath, "rep1", "speciesEcoregion_year2120.rds"))
speciesEco_focal <- readRDS(file.path(focalOutputPath, "rep1", "speciesEcoregion_year2120.rds"))

speciesEco_focal[, source := "focal"]
speciesEco_single[, source := "single"]
speciesEco <- rbind(speciesEco_focal, speciesEco_single)
speciesEco[, N := .N, .(speciesCode)]
meanMaxANPP <- ggplot(speciesEco, aes(y = maxANPP, x = speciesCode, col = source)) + 
  geom_boxplot() + 
  theme_bw() + 
  labs(x = "species", y = "maxANPP (g/m2/yr)")
ggsave("manuscript_figures/meanMaxANPP.png", meanMaxANPP, dpi = 300, width = 8, height = 4)
#I don't know if this is a useful figure - it doesn't show much different than the results one 
googledrive::drive_upload("manuscript_figures/meanMaxANPP.png", path = gFolder, name = "meanMaxANPP.png")

####total relative biomass, mean aNPP by sim ####
#summarize cohort data by pixelGroup and species, join with pixelGroupMap, summarize the landscape

focalLandscape <- lapply(list.files("outputs/focalFitting_MC/", full.names = TRUE), 
                         landscapeSummary) |>
  lapply(rbindlist) |>
  rbindlist()
singleLandscape <- lapply(list.files("outputs/singleFitting_MC/", full.names = TRUE), 
                         landscapeSummary) |>
  lapply(rbindlist) |>
  rbindlist()

#take mean of replicates
singleLandscape <- singleLandscape[, lapply(.SD, mean), 
                                   .SDcols =c ("meanANPP", "meanB", "meanMort", "propB"), 
                                   .(speciesCode, year)]
singleLandscape[, source := "single"]

focalLandscape <- focalLandscape[, lapply(.SD, mean), 
                                   .SDcols =c ("meanANPP", "meanB", "meanMort", "propB"), 
                                   .(speciesCode, year)]
focalLandscape[, source := "focal"]
landscapeDF <- rbind(singleLandscape, focalLandscape)

speciesBiomassPlot <- function(df, cols, y, species, ylab = "y",
                               plotTitle = NULL, plotSubtitle = NULL) {
  gg <- ggplot(data = df, aes(x = .data[["year"]], y = .data[[y]], fill = .data[[species]], group = .data[[species]])) +
    geom_area(position = "stack") +
    scale_fill_manual(values = cols) +
    # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    scale_y_continuous() +
    labs(x = "Year", y = ylab, fill = "Species", title = plotTitle, subtitle = plotSubtitle) +
    theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
    theme_bw(base_size = 16) + 
    facet_wrap(~source, nrow = 1)
  return(gg)
}
cols <- unique(LandR::sppEquivalencies_CA[LandR %in% landscapeDF$speciesCode & colorHex != "", 
                                          .(colorHex, LandR)])
plotCols <- cols$colorHex
names(plotCols) <- cols$LandR
speciesBiomassPlot(landscapeDF, cols = plotCols, y = "meanB", ylab = "mean landscape biomass (g/m2)", 
                   species = "speciesCode")
propB_all <- speciesBiomassPlot(landscapeDF, cols = plotCols, y = "propB", ylab = "proportion of landscape AGB", 
                   species = "speciesCode")
propB_100yr <- ggplot(landscapeDF[year == 2120], aes(y = propB, x = speciesCode, fill = source)) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  labs(y = "proportional landscape B after 100 yrs")
ggsave(filename = "manuscript_figures/propB_100yr.png", 
       propB_100yr, dpi = 300,height = 8, width = 8)

landscapeB_100yr <- ggplot(landscapeDF[year == 2120], aes(y = meanB/100, x = speciesCode, fill = source)) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  labs(y = "mean landscape B after 100 yrs (Mg/ha)")
ggsave(filename = "manuscript_figures/landscapeB_100yr.png", 
       landscapeB_100yr, dpi = 300, height = 8, width = 8)
#I don't know if this is interesting



#compare some pixelGroups
#with replicates this is done by looking at the trajectory of a pixel, which might change groups..
#TODO: wrap this in a function
#pick some pixels based on initial stand age
tempPG <- terra::rast("outputs/focalFitting_MC/rep1/pixelGroupMap_year2020.tif")
tempCD <- readRDS("outputs/focalFitting_MC/rep1/cohortData_year2020.rds")
nreps <- length(list.files("outputs/focalFitting_MC"))
ecoMap <- terra::rast("outputs/focalFitting_MC/rep1/ecoregionMap_year2120.tif")
age20PGs <- unique(tempCD[age  == c(20)]$pixelGroup)
age20PIs <- which(as.vector(tempPG) %in% age20PGs)
names(age20PIs) <- as.vector(tempPG)[age20PIs]
#randomly sample as some of the pixels have more than one pg
rm(tempPG, tempCD)
age20PIs <- sapply(age20PGs, function(pg){
  sample(age20PIs[names(age20PIs) == pg], size = 1)
})


singleCDs <- getCDLongFromOutput("outputs/singleFitting_MC", source = "single", pixelIndex = age20PIs) |>
  reproducible::Cache()
focalCDs <- getCDLongFromOutput("outputs/focalFitting_MC", source = "focal", pixelIndex = age20PIs) |>
  reproducible::Cache()
age20CDs <- rbind(focalCDs, singleCDs)

#this sums across reps, so divide by n reps
age20CDs <- age20CDs[, lapply(.SD, function(x){round(sum(x)/nreps)}), .SDcols = c("B", "mortality", "aNPPAct"), 
                     .(pixelIndex, year, speciesCode, ecoregionGroup, source)]
#TODO: confirm?
cols <- unique(LandR::sppEquivalencies_CA[LandR %in% age20CDs$speciesCode & colorHex != "", 
                                          .(colorHex, LandR)])
#TODO: make this a panel of multiple pixelIndices - 
# report the ecoregion too
# ? join speciesCode with colour
#
age20CDs <- age20CDs[cols, on = c("speciesCode" = "LandR")]
#5135 is figure prsented already

age20CDs_sample <- age20CDs[pixelIndex %in% sample(x = pixelIndex, size = 4, replace = FALSE)]
randomPixel <- ggplot(age20CDs_sample, aes(y = B, x = year, col = speciesCode, linetype = source)) + 
  geom_line() +
  theme_bw() + 
  scale_fill_manual(values = cols) + 
  labs(title = "random pixels with stand-age 20 in 2020", y = "B (g/m2)") + 
  facet_wrap(~pixelIndex, scales = "free_y")
randomPixel
ggsave("manuscript_figures/randomPixel5135.png", randomPixel, width = 8, height = 4, dpi = 300)


#what did inflationFactor look like in the theoretical species?
#get the factorial used in the simulation:
#(they are identical between reps)
factorial_spp <- list.files("outputs/focalFitting_MC/rep1", 
                            pattern = "speciesTableFactorial", full.names =TRUE)
factorial_spp <- read_ipc_file(list.files(factorial_spp, full.names = TRUE)[[1]]) |>
  as.data.table()
factorial_biomass <- list.files("outputs/focalFitting_MC/rep1", 
                                pattern = "cohortDataFactorial", full.names = TRUE)
factorial_biomass <- read_ipc_file(list.files(factorial_biomass, full.names = TRUE)[[1]]) |>
  as.data.table()
inflationFactorKey <- factorial_biomass[grep("A", speciesCode), .(maxB = max(B)), .(speciesCode)]
hist(inflationFactorKey$maxB)
inflationFactorKey <- factorial_spp[inflationFactorKey, on = c("species" = "speciesCode")]
temp <- melt.data.table(inflationFactorKey, id.vars = c("maxB", "species"), 
                        measure.vars =c ("mortalityshape", "growthcurve", "longevity", "mANPPproportion"),
                        variable.name = "trait", value.name = "value")

maxB_plot <- ggplot(temp, aes(y = maxB, x = value)) + geom_jitter() + 
  facet_wrap(~trait, scales = "free_x") + 
  theme_bw() + 
  labs(x = "maxB achieved in factorial")
ggsave(maxB_plot, file = "manuscript_figures/maxB_facetplot.png", dpi = 300, width = 8, height = 8)
googledrive::drive_upload("manuscript_figures/maxB_facetplot.png", path = gFolder, name = "maxB_facetplot.png")

# not sure if we need this
# ggsave(PSPages, filename = "manuscript_figures/MontaneCordillera_PSP_standAge_histogram.png", height = 4, width = 8, dpi = 300)
# googledrive::drive_upload("manuscript_figures/MontaneCordillera_PSP_standAge_histogram.png", path = gFolder,
#                           name = "MontaneCordillera_PSP_standAge_histogram.png")
