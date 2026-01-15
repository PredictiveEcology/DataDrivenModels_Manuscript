library(Require)
Require("data.table")
Require("reproducible")
Require("fpCompare")
Require("ggplot2")
Require("patchwork") #to stack some plots
Require("googledrive")
Require("arrow")

#these are "Methods" figures that don't need to be created from simulation output 

gFolder <- googledrive::as_id("https://drive.google.com/drive/folders/1-jrQpHIyUPsveTH5gw1fsdyd3txj8TvP?usp=drive_link")

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
  "#D55E00"  # vermillion
)

fig2a <- ggplot(data = ggData, aes(y = B, col = param, x = age)) + geom_line(linewidth = 1.2) + theme_bw() + 
  labs(y = "aboveground biomass (g/m2)", colour = "Parameter: ") + 
  scale_colour_manual(values = okabe_ito, 
                      breaks = c("aNPPAct", "mAge", "mBio", "mortality"),
                      labels = c(aNPPAct = "ANPP",
                                 mAge = "age-mortality",
                                 mBio = "development-mortality", 
                                 mortality = "total mortality")) + 
  theme(legend.position = "bottom")

fig2b <- ggplot(data = expGC, aes(y = B, x = age)) + geom_line(linewidth = 1.2) + theme_bw() + 
  labs(y = "aboveground biomass (g/m2)")

fig2 <- fig2a/fig2b
ggsave("manuscript_figures/fig2.png", fig2, dpi = 300, width = 8, height = 8)
googledrive::drive_upload("manuscript_figures/fig2.png", path = gFolder, name = "fig2.png")
rm(factorial_w_mortality, factorial_spp)
####fig 3####
#get 5 randoms and 1 non-random with excellent traits


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

# set.seed(14)
#species beginning with A are the single-species cohorts
# randomSp <- sample(factorial_spp[longevity < 400][grep("A", species),]$species, size = 6, replace = FALSE)
randomSp <- c("A265", " A1345", "A2299", "A5255", "A1256", "A688")

#because there isn't an instance of "zero" biomass that is preserved, manually add them here
# purely for illustrative purposes
fig3Dat <- factorial_biomass[speciesCode %in% randomSp] |> copy()
fig3Spp <- factorial_spp[species %in% randomSp][, .(species, longevity)] |>  copy()
fig3Spp[, c("age", "B", "speciesCode") := .(longevity, 0, species)]

fig3_data <- rbind(fig3Spp, fig3Dat, fill = TRUE)
fig3_data[, species := NULL]

#for figure legend #
temp <- factorial_spp[species %in% randomSp,] |> copy()

temp[, speciesLabel := paste0("g.c.: ", growthcurve, "; m.s.: ", mortalityshape, 
                              "; long.: ", longevity, "; mANPP: ", mANPPproportion)]
temp[, speciesLabel := paste0("Sp", 1:5, " - ", speciesLabel)]
temp <- temp[, .(species, speciesLabel)]
fig3_data <- fig3_data[temp, on = c("speciesCode" = "species")]

fig3 <- ggplot(data = fig3_data, 
               aes(y = B, x = age, col = speciesLabel)) + 
  geom_line(linewidth = 1.2, alpha = 0.7) + 
  theme_bw() + 
  labs(y = "aboveground biomass (g/m2)", 
       col = "species") + 
  theme(legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 3, byrow = TRUE))
ggsave("manuscript_figures/fig3.png", fig3, dpi = 300, width = 8, height = 4)
googledrive::drive_upload("manuscript_figures/fig3.png", path = gFolder, name = "fig3.png")
rm(fig3Dat, fig3Spp)

##### Fig 4 ####
#to redo figure 4, we need a separate factorial with all mortalityshapes, not just 21-25
# fig4_spp <- factorial_spp[longevity == 250]
# fig4_B <- factorial_biomass[speciesCode %in% fig4_spp$species]
# fig4B_mort <- fig4


