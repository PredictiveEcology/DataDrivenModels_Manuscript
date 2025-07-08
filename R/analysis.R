quickSpeciesComp

singleRep1 <- readRDS("outputs/singleFitting/rep1/species_year2091.rds")
multiRep1 <- readRDS("outputs/focalFitting/rep1/species_year2091.rds")

singleRep1[, source := "single"]
multiRep1[, source := "focal"]

speciesRep1 <- rbind(singleRep1, multiRep1) |>
  data.table::melt.data.table(id.vars = c("species", "source"),
                              measure.vars = c("mANPPproportion", "growthcurve", "inflationFactor", "mortalityshape"), 
                              variable.name = "trait", value.name = "value")
library(ggplot2)

gg <- ggplot(speciesRep1, aes(y = value, x = species, fill = source)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~trait, scales = "free") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg
