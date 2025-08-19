# This DF is a single species from the Biomass_factorial outputs. 
# When the module was run, both Biomass_speciseFactorial and Biomass_core 
# were edited to avoid removing mAge and mBio from cohortData
# This species has longevity 325, mANPP of 3.5, mortalityshape of 20, and growthcurve = 0.75. 
# Perhaps this will be unnecessary later - all other figures should come from the factorial output
mortDF <- prepInputs(destinationPath = "outputs", 
                     targetFile = "example_factorial_with_mort_subtypes.csv", 
                     fun = 'data.table::fread',
                     url = "https://drive.google.com/file/d/13kt6QYm0HQJdAX8ezPAS4-Sp3A9hPF9Q")
sppDF <- prepInputs(destinationPath = "outputs", 
                     targetFile = "example_factorialSpp.csv", 
                     fun = 'data.table::fread',
                     url = "https://drive.google.com/file/d/13kt6QYm0HQJdAX8ezPAS4-Sp3A9hPF9Q")
sppOfInterest <- sppDF[longevity == 300 & growthcurve == 0.71 & mortalityshape == 21 & mANPPproportion == 3.6,]
mortDF <- mortDF[speciesCode == sppOfInterest$species]
mortDF <- rbind(mortDF, data.table(age = 300, B = 0, aNPPAct = 0, mAge = 0, mBio = 0), fill = TRUE)
mortDF[, totalMort := sum(mBio, mAge, na.rm = TRUE), .(.I)]
#we can add a 0,0 end point for age = 300 for the plot


mortDF <- melt.data.table(mortDF, id.vars = c("age"), 
                          measure.vars = c("aNPPAct", "mAge", "mBio", "totalMort"),
                          variable.name = "param", value.name = "B")

fig1 <- ggplot(mortDF, aes(x = age, y = B/100, col = param)) + geom_line(size = 2, alpha = 0.8) + 
  labs(x = "age", y = "B (Mg/ha)", col = "Param.") + 
  scale_color_manual(values = c("#4daf4a", "#377eb8", "#984ea3", "#e41a1c"), 
                     labels = c("aNPP", "age mort.", "dev. mort.", "total mort.")) + 
  theme_bw(base_size = 14) 
ggsave(plot = fig1, filename = "outputs/mortality_example.png", 
      height = 4, width = 7)
fig1
