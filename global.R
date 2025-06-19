repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
# Need the latest version
if (tryCatch(packageVersion("SpaDES.project") < "0.1.1", error = function(x) TRUE)) {
  install.packages(c("SpaDES.project", "Require"), repos = repos)
}

projPath <- "~/git/DataDrivenModels_Manuscript"
customOpts <- list(gargle_oauth_email = "ianmseddy@gmail.com",
             gargle_oauth_cache = "~/google_drive_cache",
             gargle_oauth_client_type = "web",
             reproducible.inputPaths = "~/data"
             )
#usethis failed to install first time around...not sure why - 
inSim <- SpaDES.project::setupProject(
  packages = c("usethis", "googledrive", "httr2", "RCurl", "XML", "bcdata"),
  useGit=  TRUE,
  paths = list(projectPath = projPath,
               modulePath = file.path("modules"),
               cachePath = file.path("cache"),
               scratchPath = tempdir(),
               inputPath = file.path("inputs"),
               outputPath = file.path("outputs")
  ),
  modules = c("PredictiveEcology/Biomass_speciesFactorial",
              "PredictiveEcology/Biomass_speciesParameters",
              "PredictiveEcology/Biomass_core"),
  options = customOpts, 
  ),
  times = list(start = 2011, end = 2021),
  loadOrder = unlist(modules), #load order must be passed or BBDP will be sourced prior to fsDPF.
  #also canClimateData must come before fireSense
  climateVariablesForFire = list(ignition = c("CMDsm"),
                                 spread = c("CMDsm")),
  functions = "ianmseddy/NEBC@main/R/studyAreaFuns.R",
  # studyArea = {
  #   sa <- bcdata::
  #   sa <- sa[sa$FRU == FRU,]
  #   sa <- terra::vect(sa)
  #   sa <- terra::buffer(sa, 3000)
  #   sa <- terra::fillHoles(sa)
  # },
  # rasterToMatch = {
  #   rtm <- terra::rast(sa, res = c(250, 250))
  #   rtm[] <- 1
  #   rtm <- reproducible::postProcess(rtm, maskTo = sa)
  # },
  # rasterToMatch_biomassParam = rasterToMatch,
  # studyArea_biomassParam = studyArea,
  # sppEquiv = {
  #   spp <- LandR::speciesInStudyArea(studyArea = studyArea,
  #                                    dPath = "inputs",)
  #   sppEquiv <- LandR::sppEquivalencies_CA[KNN %in% spp$speciesList,]
  #   sppEquiv <- sppEquiv[LANDIS_traits != "",]
  #   sppEquiv <- sppEquiv[grep(pattern = "Spp", x = sppEquiv$KNN, invert = TRUE),]
  # },
  # objectSynonyms = list(c("rstLCC","rstLCC2011"),
  #                       # c("flammableRTM", "flammableRTM2011"),
  #                       # c("nonForest_timeSinceDisturbance", "nonForest_timeSinceDisturbance2011"),
  #                       # c("landcoverDT", "landcoverDT2011"),
  #                       c("standAgeMap", "standAgeMap2011")
  # ),
  params = list(
    .globals = list(.studyAreaName = paste0("FRU", FRU),
                    dataYear = 2011,
                    .plots = "png",
                    .useCache = FALSE,
                    sppEquivCol = "LandR"),
    Biomass_borealDataPrep = list(
      overrideAgeInFires = FALSE,
      overrideBiomassInFires = FALSE,
      .useCache = c(".inputObjects")
    ),
    fireSense_dataPrepFit = list(
      igAggFactor = 4),
    canClimateData = list(
      ".useCache" = c(".inputObjects"))
    # fireSense_IgnitionFit = list(".plots" = )
  )
)
