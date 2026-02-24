#most, if not all, installed from running the sim
library(Require)
Require("data.table")
Require("reproducible")
Require("fpCompare")
Require("ggplot2")
Require("patchwork") #to stack some plots
Require("googledrive")
Require("arrow")
Require("scales")

#these are "Methods" figures that don't need to be created from simulation output 

gFolder <- googledrive::as_id("https://drive.google.com/drive/folders/1-jrQpHIyUPsveTH5gw1fsdyd3txj8TvP?usp=drive_link")
singleOutputPath <- "outputs/single"
focalOutputPath <- "outputs/focal"
#Manuscript figures

####fig 2####
#by default, the factorial does not retain the development and age-based mortality columns, only total mortality
#they were retained in this specific instance by modifying the module Biomass_core to 
# preserve these columns in the mortality_and_growth event
reproducible::checkPath("manuscript_figures/", create = TRUE)
factorial_w_mortality <- prepInputs(url = "https://drive.google.com/file/d/1qbXXkbjgVgf8vAEPWgvFh3oq_imyFgpt/view?usp=drive_link", 
                                    destinationPath = "inputs", 
                                    targetFile = "example_factorial_with_mort_subtypes.csv", 
                                    fun = "data.table::fread")
factorial_spp <- prepInputs(url = "https://drive.google.com/file/d/1qbXXkbjgVgf8vAEPWgvFh3oq_imyFgpt/view?usp=drive_link", 
                            destinationPath = "inputs", 
                            targetFile = "example_factorialSpp.csv", 
                            fun = "data.table::fread")
expSpp <- factorial_spp[longevity %==% 325 & mANPPproportion == 3.6 & growthcurve == 0.77 & mortalityshape == 19,]$species
expGC <- factorial_w_mortality[speciesCode %in% expSpp,]
ggData <- melt.data.table(expGC, id.vars = "age", measure.vars = c("mortality", "aNPPAct", "mAge", "mBio"), 
                          value.name = "B", variable.name = "param")

#colourblind-friendly 
okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#D55E00",  # vermillion
  "#000"
)

mortExampleFigA <- ggplot(data = ggData, aes(y = B, col = param, x = age)) +
  geom_line(linewidth = 1.2, linetype = "twodash") +
  theme_bw(base_size = 10) + 
  geom_line(aes(x = NA, y = NA, col = "biomass"), linewidth = 1.2) + 
  labs(y = "Biomass (g/m2)", colour = "Parameter ") + 
  scale_colour_manual(values = okabe_ito, 
                      breaks = c("aNPPAct", "mAge", "mBio", "mortality", "biomass"),
                      labels = c(aNPPAct = "ANPP",
                                 mAge = "age-mortality",
                                 mBio = "development-mortality", 
                                 mortality = "total mortality",
                                 biomass = "total biomass"))
  mortExampleFigB <- ggplot(data = expGC, aes(y = B, x = age)) +
  geom_line(linewidth = 1.2) + 
  theme_bw(base_size = 10) +
  scale_fill_discrete(name = "Biomass") + 
  labs(y = "Biomass (g/m2)")

mortExampleFig <- mortExampleFigA/mortExampleFigB
ggsave("manuscript_figures/mortExampleFig.png", mortExampleFig, dpi = 600, width = 6, height = 5)
googledrive::drive_upload("manuscript_figures/mortExampleFig.png", path = gFolder, name = "mortExampleFig.png", 
                          overwrite = TRUE)
rm(factorial_w_mortality, factorial_spp)

####fig 3####
#get 5 randoms and 1 non-random with excellent traits
#get the factorial used in the simulation:
# the name will vary so grep it
factorial_spp <- list.files(focalOutputPath, 
                            pattern = "speciesTableFactorial", full.names =TRUE) |>
  list.files(full.names = TRUE) |>
  read_ipc_file() |>
  as.data.table()

factorial_biomass <- list.files(focalOutputPath, 
                                pattern = "cohortDataFactorial", full.names =TRUE) |>
  list.files(full.names = TRUE) |>
  read_ipc_file() |>
  as.data.table()

set.seed(15)
#species beginning with A are the single-species cohorts
randomSp <- sample(factorial_spp[longevity < 400][grep("A", species),]$species, size = 6, replace = FALSE)
# randomSp <- c("A4037", "A675", "A2978", "A2598", "A177", "A261")
#because there isn't an instance of "zero" biomass that is preserved, manually add them here
# purely for illustrative purposes
randomGCFig_Dat <- factorial_biomass[speciesCode %in% randomSp] |> copy()
randomGCFig_Spp <- factorial_spp[species %in% randomSp][, .(species, longevity)] |>  copy()
randomGCFig_Spp[, c("age", "B", "speciesCode") := .(longevity, 0, species)]

randomGCFig_data <- rbind(randomGCFig_Spp, randomGCFig_Dat, fill = TRUE)
randomGCFig_data[, species := NULL]

#for figure legend #
temp <- factorial_spp[species %in% randomSp,] |> copy()
temp[, speciesLabel := paste0("g.c.: ", growthcurve, "; m.s.: ", mortalityshape, 
                              "; long.: ", longevity, "; mANPP: ", mANPPproportion)]
temp[, speciesLabel := paste0("Sp", 1:5, " - ", speciesLabel)]
temp <- temp[, .(species, speciesLabel)]
randomGCFig_data <- randomGCFig_data[temp, on = c("speciesCode" = "species")]

randomGCFig <- ggplot(data = randomGCFig_data, 
                      aes(y = B, x = age, col = speciesLabel)) + 
  geom_line(linewidth = 1.2, alpha = 0.7) + 
  theme_bw(base_size = 10) + 
  labs(y = "biomass (g/m2)", 
       col = "species") + 
  theme(legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 3, byrow = TRUE))
ggsave("manuscript_figures/randomGCFig.png", randomGCFig, dpi = 600, width = 7, height = 4)
googledrive::drive_upload("manuscript_figures/randomGCFig.png", path = gFolder, name = "randomGCFig.png", 
                          overwrite = TRUE)
rm(randomGCFigDat, randomGCFig_Spp)

##### Fig 4 ####
# fig 4 requires a separate factorial with all mortalityshapes, not just 21-25
# this was created by allowing fewer values of longevity and running with single
fT <- prepInputs(url = 'https://drive.google.com/file/d/1AeCMmh_rtzxXmqAU8od2b_w5g2cexDRt/view?usp=drive_link', 
                 targetFile = "factorialTraits_allMortalityShape.csv", fun = "data.table::fread", 
                 destinationPath = "outputs")
fB <- prepInputs(url = "https://drive.google.com/file/d/1AeCMmh_rtzxXmqAU8od2b_w5g2cexDRt/view?usp=drive_link", 
                 targetFile = "factorialBiomass_allMortalityShape.csv", 
                 destinationPath = "outputs", fun = "data.table::fread")
#keep longevity constant for figures
fT250 <- fT[longevity == 250,]
fB250 <- fB[speciesCode %in% fT250$species]
fT250 <- fT250[, .(species, growthcurve, mortalityshape, mANPPproportion, pixelGroup)]
fB250 <- fB250[fT250, on = c("speciesCode" = "species", "pixelGroup")]

#for each param of interest, hold the others constant
# these are the mean values in the factorial (mANPP ranged 2.4 to 6.4)
growthcurveSpp <- fT250[mortalityshape == 15 & mANPPproportion == 4.4,]$species
mortalityshapeSpp <- fT250[growthcurve == 0.5 & mANPPproportion == 4.4,]$species
mANPPpropSpp <- fT250[mortalityshape == 15 & growthcurve  == 0.5,]$species


leg_theme <- theme(
  legend.position = "bottom",
  legend.key.height = unit(2.5, "mm"),
  legend.key.width  = unit(3.0, "mm"),
  legend.spacing.y  = unit(1.5, "mm"),
  legend.margin     = margin(t = 0, r = 0, b = 0, l = 0),
  legend.box.margin = margin(t = -4, r = 0, b = 0, l = 0)
)

#it is much easier to let ggplot do the faceting



fig4A <- ggplot(fB250[speciesCode %in% growthcurveSpp, ],
                aes(age, B, group = speciesCode, colour = growthcurve)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  labs(y = "Biomass (g/m2)", x = "") +
  scale_colour_viridis_c(
    option = "D",
    breaks = pretty_breaks(3),
    guide = guide_colorbar(direction = "vertical")
  ) +
  theme_bw(base_size = 10) +  
  coord_cartesian(ylim = c(0,5000)) +
  leg_theme

fig4B <- ggplot(fB250[speciesCode %in% mortalityshapeSpp, ],
                aes(age, B, group = speciesCode, colour = mortalityshape)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  labs(x = "age", y = "") +
  scale_colour_viridis_c(
    option = "D",
    breaks = pretty_breaks(3),
    guide = guide_colorbar(direction = "vertical")
  ) +
  theme_bw(base_size = 10) + 
  coord_cartesian(ylim = c(0,5000)) +
  leg_theme + theme(axis.text.y = element_blank())

fig4C <- ggplot(fB250[speciesCode %in% mANPPpropSpp, ],
                aes(age, B, group = speciesCode, colour = mANPPproportion)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  labs(y = "", x = "") +
  scale_colour_viridis_c(
    option = "D",
    breaks = pretty_breaks(3),
    guide = guide_colorbar(direction = "vertical")
  ) +
  theme_bw(base_size = 10) + 
  coord_cartesian(ylim = c(0,5000)) +
  leg_theme + 
  theme(axis.text.y = element_blank())

fig4 <- fig4A + fig4B + fig4C

ggsave(
  filename = "manuscript_figures/fig4_growthcurveparam_fullRange.png",
  plot = fig4,
  width = 7, height = 5, dpi = 600
)

googledrive::drive_upload("manuscript_figures/fig4_growthcurveparam_fullRange.png", overwrite = TRUE,
                          path = gFolder, name = "fig4_growthcurveparam_fullRange.png")

#unsure if this is a methods or results figure
factorial_spp <- list.files(focalOutputPath, 
                            pattern = "speciesTableFactorial", full.names =TRUE)
factorial_spp <- read_ipc_file(list.files(factorial_spp, full.names = TRUE)[[1]]) |>
  as.data.table()
factorial_biomass <- list.files(focalOutputPath, 
                                pattern = "cohortDataFactorial", full.names = TRUE)
factorial_biomass <- read_ipc_file(list.files(factorial_biomass, full.names = TRUE)[[1]]) |>
  as.data.table()
inflationFactorKey <- factorial_biomass[grep("A", speciesCode), .(maxB = max(B)), .(speciesCode)]

inflationFactorKey <- factorial_spp[inflationFactorKey, on = c("species" = "speciesCode")]
temp <- melt.data.table(inflationFactorKey, id.vars = c("maxB", "species"), 
                        measure.vars =c ("mortalityshape", "growthcurve", "longevity", "mANPPproportion"),
                        variable.name = "trait", value.name = "value")

maxB_plot <- ggplot(temp, aes(y = maxB, x = value)) + geom_jitter() + 
  facet_wrap(~trait, scales = "free_x") + 
  theme_bw(base_size =10) + 
  labs(y = "maxB achieved in factorial")
ggsave(maxB_plot, file = "manuscript_figures/maxB_facetplot.png", dpi = 600, width = 7, height = 7)
googledrive::drive_upload("manuscript_figures/maxB_facetplot.png", path = gFolder, name = "maxB_facetplot.png", 
                          overwrite = TRUE)
