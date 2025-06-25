# This DF is a single species from the Biomass_factorial outputs. 
# When the module was run, both Biomass_speciseFactorial and Biomass_core 
#were edited to avoid removing mAge and mBio from cohortData
# This species has longevity 325, mANPP of 3.5, mortalityshape of 20, and growthcurve = 0.75. 
# Perhaps this will be unnecessary later - all other figures should come from the factorial output
mortDF <- fread("misc/mortCurve_325lng_3-5aNPP_20ms_75gc.csv")
mortDF[, totalMort := sum(mBio, mAge, na.rm = TRUE), .(.I)]
mortDF <- melt.data.table(mortDF, id.vars = c("age"), 
                          measure.vars = c("aNPPAct", "mAge", "mBio", "totalMort"),
                          variable.name = "param", value.name = "B")
ggplot(mortDF, aes(x = age, y = B/100, col = param)) + geom_line(size = 1.5, alpha = 0.9) + 
  labs(x = "age", y = "B (Mg/ha)") + 
  scale_color_manual(values = c("#4daf4a", "#377eb8", "#984ea3", "#e41a1c"), 
  labels = c("aNPP", "age mort.", "dev. mort.", "total mort.")) + 
  theme_bw()

