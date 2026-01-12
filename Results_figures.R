#These are Results figures created from simulation output


Require::Require(data.table)
Require::Require(ggplot2)

singleOutputPath <- "outputs/singleFitting_MC"
focalOutputPath <- "outputs/focalFitting_MC"
#no stochasticity in species traits, so use rep 1 for simplicity
singleRep1 <- readRDS(file.path(singleOutputPath, "rep1/species_year2120.rds"))
focalRep1 <- readRDS(file.path(focalOutputPath, "rep1/species_year2120.rds"))

singleRep1[, source := "single"]
focalRep1[, source := "focal"]
speciesRep1 <- rbind(singleRep1, focalRep1) 
speciesRep1[, maxB_multiplier := inflationFactor - 1]
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

#mANPPproportion is actually an integer (percentage) for memory reasons
#It has been set as a proportion here to reduce confusion
# ggplot(speciesRep1, aes(y = mANPPproportion/100, x = species, fill = source)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   theme_bw() +
#   labs(y = "mANPPproportion") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(speciesRep1, aes(y = growthcurve, x = species, fill = source)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   theme_bw() +
#   labs(y = "growthcurve") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# ggplot(speciesRep1, aes(y = mortalityshape, x = species, fill = source)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   theme_bw() +
#   labs(y = "mortalityshape") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# # #center inflation factor around 0 - instead of one (or set ymin to be 1)
# # ylims <- c(1, max(speciesRep1$inflationFactor) + 0.02)
# 
# speciesRep1
# 
# ggplot(speciesRep1, aes(y = c(inflationFactor - 1) * 100, x = species, fill = source)) + 
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_bw() +
#   labs(y = "maxB multiplier (%)") + 
#   # coord_cartesian(ylim = ylims) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



speciesEco_single <- readRDS(file.path(singleOutputPath, "rep1", "speciesEcoregion_year2120.rds"))
ecoMap <- terra::rast("outputs/focalFitting_MC/rep1/ecoregionMap_year2120.tif")

ecoMap
# ggplot(data = speciesEco_single, aes(y = maxB/100, x = speciesCode)) + 
#   geom_bar(stat = "identity", position = "dodge") + 
#   labs(x = "species", y = "maxB (") +
#   facet_wrap(~ecoregionGroup,scales = "free") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#comapre species within some BECs 

# plot(inSim$ecoregionLayer["ZONE_NAME"])
# speciesEco_single[, .N, .(ecoregionGroup)]
# 
# ####compare the distribution of BEC subzones
# #use `cats`
# ecozoneMap <- rast(file.path(singleOutputPath, "rep1", "ecoregionMap_year2100.tif"))
# terra::plot(ecozoneMap)
